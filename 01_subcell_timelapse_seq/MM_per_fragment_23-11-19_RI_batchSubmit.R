#!/n/app/R/4.3.1-gcc-9.2.0/bin/Rscript

# Written by Brendan and Kate 2/2020
# Modified 10/2020 by Mary to write fragment as sequence instead of read1 ** uses read1 sequence by default in overlapping region even if they differ there
# Modified 11/2023 by Robert Ietswaart for batch running with R/4.3.1: automate necessary package installation

args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
path <- paste0(args[2], '/')
fname <- args[3]

# install and load packages
lib_folder <- "~/R/rna_flow_packages/"
# Create the folder if it does not exist
if (!dir.exists(lib_folder)) {
  dir.create(lib_folder, recursive = TRUE)
}

# Function to check if a package is installed
is_package_installed <- function(pkg_name, lf) {
  return(pkg_name %in% rownames(installed.packages(lib.loc = lf)))
}

# Install package if not installed in the specified library folder
if (!is_package_installed("tidyr", lib_folder)) {
   options(repos = c(CRAN = "https://cloud.r-project.org"))
   install.packages("tidyr", lib=lib_folder)
}

if (!is_package_installed("splitstackshape", lib_folder)) {
   options(repos = c(CRAN = "https://cloud.r-project.org"))
   install.packages("splitstackshape", lib=lib_folder)
}

library(tidyr, lib.loc=lib_folder)
library(splitstackshape, lib.loc=lib_folder)
library(data.table, lib.loc=lib_folder)

# output file name
new_filename = paste0(path, dat_file, "_analysis/", fname, "_fragments.bed")

# load current mismatch data, which must be sorted by read name
table <- read.csv(paste0(path, dat_file), sep="\t", header=FALSE, stringsAsFactors = FALSE)

table_new <- table[order(table$V4),]

table <-table_new

# change classes
table$V11 <- as.character(table$V11)
table$V12 <- as.character(table$V12)
table$V13 <- as.character(table$V13)
table$V14 <- as.character(table$V14)
table$V15 <- as.character(table$V15)
table$V16 <- as.character(table$V16)
table$V17 <- as.character(table$V17)

# get number of reads
number_of_reads <- nrow(table)

# just run this to set up new table, this will be overwritten later. 
new_table<-table[1,]

# set this counter 
i<-1

# for every pair of reads (row number 1, 3, 5, etc), do the following:
for (x in seq(from=1, to=number_of_reads, by=2)) {
  
  # for each read, make a table containing the genome positions each cover and the corresponding read nt quality
    read_1 <- data.frame(strsplit(table[x,16], ", "), strsplit(table[x,17], ", "))
    read_2 <- data.frame(strsplit(table[x+1,16], ", "), strsplit(table[x+1,17], ", "))
    
    # save how many genomic nts are covered with this pair
    all_positions_covered=c(as.numeric(as.character(read_1[,1])), as.numeric(as.character(read_2[,1])))
    sum_positions_covered=length(unique(all_positions_covered))
  
    # Get fragment sequence
    read1 <- as.numeric(unlist(strsplit(table[x,16], ", ")))
    read2 <- as.numeric(unlist(strsplit(table[x+1,16], ", ")))
    seq1 <- as.character(table[x,9])
    seq2 <- as.character(table[x+1,9])
    if (table[x,6] == '+') {
    addR2toEnd = length(setdiff(read2, read1))
    R2toAdd = substring(seq2, nchar(seq2)-addR2toEnd+1, nchar(seq2))
    frag = paste0(seq1, R2toAdd)
    start = table[x,2]
    end = table[x+1,3]
    }
    if (table[x,6] == '-') {
    addR2toFront = length(setdiff(read2, read1))
    R2toAdd = substring(seq2, 1, addR2toFront)
    frag = paste0(R2toAdd, seq1)
    start = table[x+1,2]
    end = table[x,3]
    }

  # if the start coordinate of read 1 matches the start coordinate of read 2, AND they have the same CIGAR, then they are already identical
  # in this case, just set the whole sam line for read 1 as the fragment info
  if  (table[x,2] == table[(x+1),2] && table[x,8] == table[(x+1),8]) {
    new_table[i,]=table[x,]
    # But change seq to fragment instead of read1
    new_table[i,9]=frag
    # And change coordinates
    new_table[i,2]=start
	new_table[i,3]=end
    # column 12 is now number of genomic nts that are covered by this pair
    new_table[i,12]=sum_positions_covered
    
    # otherwise, the reads are different
  } else {
    
    #first calculate total number of mismatches between reads 1 and 2
    num_MM<- table[x,18]+table[(x+1),18]
    
    # if there are no mismatches in either read, then it doesn't matter, just keep read 1 information
    if (num_MM == 0) {
      new_table[i,]=table[x,]
      # But change seq to fragment instead of read1
      new_table[i,9]=frag
      # And change coordinates
      new_table[i,2]=start
	  new_table[i,3]=end

      # column 12 is now number of genomic nts that are covered by this pair
      new_table[i,12]=sum_positions_covered
      
    } else {
      
      # if there are mismatches in at least one of the reads, 
      # make data frame with mismatches from read 1, remove NAs, and separate each mismatch into its own row
      mismatch_table_r1<-data.frame(table[x,c(11,13,14,15)], stringsAsFactors = FALSE)
      mismatch_table_r1<-cSplit(mismatch_table_r1, splitCols=c("V11","V13","V14","V15"), sep=", ", direction="long", makeEqual = FALSE, type.convert=FALSE)
      label<-rep("r1",nrow(mismatch_table_r1))
      mismatch_table_r1<-cbind(label, mismatch_table_r1)
      
      # make the same table for read 2
      mismatch_table_r2<-data.frame(table[x+1,c(11,13,14,15)], stringsAsFactors = FALSE)
      mismatch_table_r2<-cSplit(mismatch_table_r2, splitCols=c("V11","V13","V14","V15"), sep=", ", direction="long", makeEqual = FALSE, type.convert=FALSE)
      label<-rep("r2",nrow(mismatch_table_r2))
      mismatch_table_r2<-cbind(label, mismatch_table_r2)
      
      # combine mismatches from both reads, and get rid of NAs (if one read contains no MM)
      mismatch_table_both<-rbind(mismatch_table_r1, mismatch_table_r2)
      mismatch_table_both<-mismatch_table_both[complete.cases(mismatch_table_both), ]
      
      # identify mismatches that are found on both reads and move to new table
      mismatch_keep<-unique(mismatch_table_both[duplicated(mismatch_table_both[,3:5]),])
      
      # move mismatches that are not found on both reads to new table
      mismatch_examine<-mismatch_table_both[!(duplicated(mismatch_table_both[,3:5]) | duplicated(mismatch_table_both[,3:5], fromLast = TRUE)), ]
      
      # if there are any...
      if (nrow(mismatch_examine)>0) {
      
      # examine each of these
      for (k in seq(1, nrow(mismatch_examine), by=1)) {
        
        # if the mismatch falls on read 1...
        if (mismatch_examine[k,1]=="r1") {
          
          # check to see if it is covered by read 2
          if (mismatch_examine[k,4] %in% read_2[,1]) {
            
            # if it is on read 2, get the row
            read_2_row <- which(read_2[,1]==as.numeric(mismatch_examine[k,4]))
            
            # now get the quality
            read_2_qual <- read_2[read_2_row,2]
            read_2_qual <- as.numeric(as.character(read_2_qual))
            
            # if read 1 has a higher quality...
            if (as.numeric(mismatch_examine[k,2])>read_2_qual) {
              
              # then keep the mismatch
              mismatch_keep = rbind(mismatch_examine[k,], mismatch_keep)
              
            # otherwise read 2 has an equal or higher quality, so don't keep that mismatch. 
            } else {
            } # end comparing read quality
            
            # otherwise it isn't covered by read 2, so keep the mismatch from read 1  
             } else {
              mismatch_keep = rbind(mismatch_examine[k,], mismatch_keep)
              
          } # end check if the mismatch on read 1 is covered by read 2
          
          # otherwise the mismatch falls on read 2, so do the opposite (check read 1)
        } else {
          
          # check to see if it is covered by read 1
          if (mismatch_examine[k,4] %in% read_1[,1]) {
            
            # if it is on read 1, get the row
            read_1_row <- which(read_1[,1]==as.numeric(mismatch_examine[k,4]))
            
            # now get the quality
            read_1_qual <- read_1[read_1_row,2]
            read_1_qual <- as.numeric(as.character(read_1_qual))
            
            # if read 2 has a higher quality...
            if (as.numeric(mismatch_examine[k,2])>read_1_qual) {
              
              # then keep the mismatch
              mismatch_keep = rbind(mismatch_examine[k,], mismatch_keep)
              
              # otherwise read 1 has an equal or higher quality, so don't keep that mismatch. 
            } else {
            } # end comparing read quality
            
            # otherwise it isn't covered by read 2, so keep the mismatch from read 1  
          } else {
            mismatch_keep = rbind(mismatch_examine[k,], mismatch_keep)
            
          } # end check if the mismatch on read 1 is covered by read 2 
          
        }
        
      } # end loop for examine each mismatch not in each read
        
        # if there are no remaining mismatches to examine, do nothing
      } else {
        
      } # end loop for mismatch_examine table
      
      # create a mismatch table for all mismatches that will be kept
      mismatch_table <- mismatch_keep[,2:5]
      
      # if there are mismatches to be kept
      if (nrow(mismatch_table)>0) {
        
        # if there any cases where the two reads have different mismatches at the same genomic position with equal quality scores (likely very rare),
        # then just keep the first one (aka, randomly pick one)
       mismatch_table<-mismatch_table[!duplicated(mismatch_table[,3:4]),]
        
      # which of these mismatches are T>C mismatches?
      num_TC<-0
      for (c in seq(from=1, to=nrow(mismatch_table))){
        if (mismatch_table[c,2]=="C" & mismatch_table[c,4]=="T") {
          num_TC<-num_TC+1
        }
      }
      
      # does the fragment contain a T>C mismatch?
      read_TC<-0
      if (num_TC>0) {
      read_TC<-read_TC+1
      }
      
      # transpose so you can use unite function
      mismatch_table<-transpose(mismatch_table)
      
      # now concat back together
      cat_mismatch_table<- unite(mismatch_table, col="new", sep=", ")
      
      # now store this new mismatch info as V11, V12, and V13 in new_table
      # first copy original sam entry
      new_table[i,]<-table[x,]
      # But change seq to fragment instead of read1
      new_table[i,9]=frag
      # And change coordinates
      new_table[i,2]=start
      new_table[i,3]=end

      # column 11 is quality of mismatches
      new_table[i,11]=cat_mismatch_table[1,1]
      # column 12 is now number of genomic nts that are covered by this pair
      new_table[i,12]=sum_positions_covered
      new_table[i,13]=cat_mismatch_table[2,1]
      new_table[i,14]=cat_mismatch_table[3,1]
      new_table[i,15]=cat_mismatch_table[4,1]
      # fix total number of mismatches
      new_table[i,18]=ncol(mismatch_table)
      # report new number of T>C mismatches
      new_table[i,19]=num_TC
      # report of the fragment contains a T>C mismatch (binary)
      new_table[i,20]=read_TC
      
      # otherwise there are no mismatches to be kept
      } else {
        
        # first copy original sam entry
        new_table[i,]<-table[x,]
        # But change seq to fragment instead of read1
        new_table[i,9]=frag
        # And change coordinates
        new_table[i,2]=start
	    new_table[i,3]=end

        # column 11 is quality
        new_table[i,11]=NA
        # now replace with new mismatches
        new_table[i,12]=sum_positions_covered
        new_table[i,13]=NA
        new_table[i,14]=NA
        new_table[i,15]=NA
        # fix total number of mismatches
        new_table[i,18]=0
        # report new number of T>C mismatches
        new_table[i,19]=0
        # report of the fragment contains a T>C mismatch (binary)
        new_table[i,20]=0
        
        
      } # end loop for are there any mismatches to be kept
    }
  
  } 
    i<-i+1
}

# save ouput
write.table(new_table, file=new_filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
