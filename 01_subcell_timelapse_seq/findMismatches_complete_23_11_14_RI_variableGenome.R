#!/n/app/R/4.3.1-gcc-9.2.0/bin/Rscript

# Interpret cigar strings of sam file to find mismatches in reads and genomic location when mapped to genome
# For Brendan, written on 7/10/18 by K. Lachance

# NOTE: Modified in 9/18/18 by K. Lachance to also look-up genomic nt, as well as calculate total number of mismatches and T>C mismatches specifically in each read

# NOTE: Modified on 2/4/19 by K. Lachance to also calculate total number of bases aligning in each read (sum of all values in cigar string, ignoring softclipped S)

# NOTE: Modified on 3/21/19 by K. Lachance to fix strand-flipping

# NOTE: Modified on 6/14/19 by Brendan to include "N" in CIGAR which indicates a spliced read aligned by STAR (treat as a deletion), new genome fasta (combined mouse and fly), and new SAM flags to indicate strand

# NOTE: Modified on 7/26/19 by K. Lachance to re-fix strand-flipping for new library prep protocol

# NOTE: Modified on 10/22/19 by Brendan to work with pseudo-bed files instead of raw sam alignment

# NOTE: Modified on 10/23/19 by Brendan to fix cases where bedtools getfasta fails because genome coordinates are in scientific notation

# NOTE: Modified on 3/15/21 by M. Couvillion to add fasta path as command line input

# NOTE: Modified on 11/14/2023 by R. Ietswaart to make compatible with R/4.3.1 and batch submissions, removed dependency on gtools.

# Load bash modules for later analysis
system("module load gcc/9.2.0")
system("module load bedtools/2.30.0")

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
path <- paste0(args[2],'/')
fastapath <- args[3]

# Create output file name
filename = paste0(path, dat_file, "_MM.bed") # For running in parallel

# Load data
cat("Loading data into R...\n")
dat <- read.table(paste0(path, dat_file), fill = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL, header = FALSE)

# Remove rows that are all NAs but first column (artifact of import)
dat = dat[dat[,9] != "",]

# Create output table to hold new information
m = nrow(dat)
out <- cbind(dat, matrix(NA, nrow = m, ncol = 11)) # 11 extra columns will hold read nt quality, read nt position, read nt, and genome pos, genome nt, all genome positions covered by the read, #MM, #T>C, and a binary indication if the read contains a T>C, and total aligned bases, respectively

# For each read, parse cigar string
cat("Running through each read to find mismatches...\n")
for (i in 1:m) {
    
    # Track progress
    if (i %% 1000 == 0) {
       cat(sprintf("Done with %i reads...\n", i))
    }

    # Create empty lists to store valuable information for each read
    read_pos = c()
    read_let = c()
    read_quality = c()
    gen_pos = c()
    gen_let = c()
    gen_pos_covered = c()
    all_covered_quality = c()

    # Get cigar string (in column 6)
    cig = out[i,8]

    # Get the read contents (in column 9, will be useful for getting mismatched letters later)
    read = out[i, 9]
    
    # Get the quality scores
    quality_scores_string = out[i,10]

    # Parse cigar string by delimiters (=, D, I, X, N, and S) but keep delimiters attached to value
    C = unlist(strsplit(cig, "(?<=[=DIXNS])", perl = TRUE))

    # Have counters for position in read and genome
    p_read = 1 # Counter for current position in read, start at the beginning for each read
    p_gen = as.numeric(out[i, 2]) # Counter for current position in the genome, start at genomic position at beginning of each read (in column 4)
    total_align = 0 # Total number of bases that align from read

    # Step through each part of cigar string for each read
    for (j in 1:length(C)) {
    	
	# Separate each step into number and delimiter
	tmp = substr(C[j], 0, nchar(C[j])-1)
	tmp = c(tmp, substr(C[j], nchar(C[j]), nchar(C[j])))

	# If the delimiter is =, then just add that to the position counter (both)
	if (tmp[2] == "=") {
	  gen_pos_covered = c(gen_pos_covered, seq(p_gen, as.numeric((p_gen+as.numeric(tmp[1])-1)), by=1))
	  all_covered_quality = c(all_covered_quality, as.character(substr(quality_scores_string, p_read, (p_read + as.numeric(tmp[1])-1))))
	   p_read = p_read + as.numeric(tmp[1])
	   p_gen = p_gen + as.numeric(tmp[1])
	   total_align = total_align + as.numeric(tmp[1])
	}

	# If the delimiter is D, then just add that to the genome position counter (not the read position counter)
	if (tmp[2] == "D") {
	  #gen_pos_covered = c(gen_pos_covered, seq(p_gen, as.numeric((p_gen+as.numeric(tmp[1])-1)), by=1))
	   p_gen = p_gen + as.numeric(tmp[1])
	   # total_align = total_align + as.numeric(tmp[1]) # Remove if don't want to count deletions toward total read aligned count
	}
	
	# If the delimiter is N, then just add that to the genome position counter (not the read position counter)
	if (tmp[2] == "N") {
	   p_gen = p_gen + as.numeric(tmp[1])
	   # total_align = total_align + as.numeric(tmp[1]) 
	}

	# If the delimiter is I, then just add that to the read position counter (not the genome position counter)
	if (tmp[2] == "I") {
	   p_read = p_read + as.numeric(tmp[1])
	   total_align = total_align + as.numeric(tmp[1]) # Remove if don't want to count insertions toward total read aligned count
	}
	
	# If the delimiter is S, these have been soft-clipped, which bowtie2 ignores when reporting the genomic coordinates where the read
	# aligns to in the sam file, so then just add that to the read position counter as well
	if (tmp[2] == "S") {
	   p_read = p_read + as.numeric(tmp[1])
	}

	# If the delimiter is X, then that is a mismatch!
	if (tmp[2] == "X") {
	
	   # Now, the number associated with the X indicates the number of consecutive mismatches, so we need to step through each one
	   for (k in 1:as.numeric(tmp[1])) {
	   
		   # Store the current p_read and p_gen values as positions of mismatches
	   	   read_pos = c(read_pos, p_read)
	   	   gen_pos = c(gen_pos, p_gen)

	   	   # Find the character of the sequence quality of the read
	   	   q_score_char = as.character(substr(quality_scores_string, p_read, p_read))
	   	   
	   	   # Get the ASCII integer that corresponds to the quality characters
	   	   q_score <- utf8ToInt(q_score_char)
	   	   
	   	   # Get the actual quality score from the ASCII integer
	   	   q_score_real <- q_score-33
	   	   
	   	   # Record quality score
	   	   read_quality = c(read_quality, q_score_real)
	   	   
	   	   all_covered_quality = c(all_covered_quality, q_score_char)
	   	   
	   	   # Next, find the letter associated with the current read position
		    mm.old = substr(read, p_read, p_read)
		   
		   #  Now, look up the corressponding genomic nucleotide at the mismatch position
                   if (out[i, 7] == 147 | out[i, 7] == 99 | out[i, 7] == 72 | out[i, 7] == 152 | out[i, 7] == 0) { # Strand information is held in the second column of the sam file (163 and 99 for +, 147 and 83 for -)
                      strand = "+"
                   } else {
                      strand = "-"
                   }
                   
          #   If on the - strand, flip the letter
		   if (strand == "-") {
		      if (toupper(mm.old) == "A") {
		      	 mm.new = "T"
		      }
		      if (toupper(mm.old) == "C") {
   		        mm.new	= "G"
		      }
		      if (toupper(mm.old) == "G") {
		      	mm.new = "C"
		      }
		      if (toupper(mm.old) == "T") {
   		        mm.new	= "A"
	              }
       		   }
		   if (strand == "+") {
		      mm.new = toupper(mm.old)
		   }

			read_let = c(read_let, mm.new)

		   # Create a data frame with just one row with the information required to look up the genomic nucleotide corresponding to the mismatch
		   t = data.frame(Read = out[i, 1], Start = gen_pos[length(gen_pos)], End = gen_pos[length(gen_pos)] + 1, Space1 = ".", Space2 = "0", Strand = strand, stringsAsFactors = FALSE)
		   t$Start<-format(t$Start, scientific=F)
		   t$End<-format(t$End, scientific=F)
		   write.table(t, paste0(path, dat_file,"_analysis/temp1.bed"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

		   # Run bedtools to get the corresponding genomic nucleotide
		   system(paste0("bedtools getfasta -tab -s -fi ", fastapath," -bed ", path, dat_file,"_analysis/temp1.bed > ", path, dat_file,"_analysis/temp2.txt"))

		   # Read in genomic nucleotide
		   nt.old = read.table(paste0(path, dat_file,"_analysis/temp2.txt"), colClasses = c("character", "character"))[1,2]

		   nt.new = toupper(nt.old)

		   # Record genomic nucleotide
		   gen_let = c(gen_let, nt.new)

		   # Clean-up temporary files
		   system(paste0("rm -fr ", path, dat_file,"_analysis/temp1.bed ", path, dat_file,"_analysis/temp2.txt"))
		   
		   # add one to both position counters (once for each consecutive mismatch)
		   gen_pos_covered = c(gen_pos_covered, p_gen) 	   
		   p_read = p_read + 1
		   p_gen = p_gen + 1
		   total_align = total_align + 1

            } # End of stepping through each consecutive mismatch

	} # End of if statement for whether cigar sub-string is indicating a mismatch

    } # End of for loop stepping through a cigar string

    # Calculate total read statistics: number of mismatches, number of T>C mismatches, and the binary switch if the read contains T>C mismatches

    nMM = length(read_pos) # Number of mismatches

    nTC = 0 # Number of T>C mismatches
    TC = 0 # Flag for if read contains T>C mismatches
    
    # Convert all covered quality string into score
    
    all_covered_quality <- paste(all_covered_quality, collapse='')
    
    all_covered_quality_scores <- utf8ToInt(all_covered_quality)-33

    # Step through each mismatch and check if it is a T>C mismatch
        if (nMM > 0) {
            for (l in 1:nMM) {
      		  if (read_let[l] == "C" & (gen_let[l] == "T")) {
                       nTC = nTC + 1
                       TC = 1
                  }
            }
        } 

    # When entire cigar string has been processed, check to make sure that read_pos has reached the length of the read
    #if (p_read != nchar(read)) {
       #cat(sprintf("WARNING: Read #%i does not have a matching cigar string and read length!\n", i))
    #}

    # If mismatch lists have been populated, save lists in dataframe
    if (nMM > 0) {
        # Save lists of mismatches to last three columns
    	out[i, ncol(out)-10] = paste(read_quality, collapse=', ')
      	out[i, ncol(out)-9] = paste(read_pos, collapse=', ')
    	out[i, ncol(out)-8] = paste(read_let, collapse=', ')
    	out[i, ncol(out)-7] = paste(gen_pos, collapse=', ')
	out[i, ncol(out)-6] = paste(gen_let, collapse=', ')
	out[i, ncol(out)-5] = paste(gen_pos_covered, collapse=', ')
	out[i, ncol(out)-4] = paste(all_covered_quality_scores, collapse=', ')
	out[i, ncol(out)-3] = nMM
	out[i, ncol(out)-2] = nTC
	out[i, ncol(out)-1] = TC
	out[i, ncol(out)] = total_align

    } else {
      	# Save NAs
        out[i, ncol(out)-10] = NA
        out[i, ncol(out)-9] = NA
        out[i, ncol(out)-8] = NA
        out[i, ncol(out)-7] = NA
        out[i, ncol(out)-6] = NA
        out[i, ncol(out)-5] = paste(gen_pos_covered, collapse=', ')
        out[i, ncol(out)-4] = paste(all_covered_quality_scores, collapse=', ')
        out[i, ncol(out)-3] = 0
        out[i, ncol(out)-2] = 0
        out[i, ncol(out)-1] = 0
        out[i, ncol(out)] = total_align

    } # End of if statement to save lists to data frame

} # End of for loop stepping through each read

# Write file to output file
cat("Writing results to output file...\n")
write.table(out, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

