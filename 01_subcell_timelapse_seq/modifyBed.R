#!/n/app/R/3.5.1/bin/Rscript

# use: ./modifyBed.R .bed .sam

# set bed and sam files as arguments 1 and 2, respectively
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
sam_file <- args[2]

# set new file name
new_file_name <- paste(bed_file, "withC.bed", sep="")

# count the max number of columns in sam file (since this varies based on single or paired end data)
nF <- max(count.fields(sam_file, sep = "\t"))

# load in current bed and sam files
bed_original <- read.table(bed_file, fill = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL, header = FALSE)
sam_original <- read.table(sam_file, fill = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL, header = FALSE, col.names = paste0("V",seq_len(nF)))

new_bed<-cbind(bed_original$V1, bed_original$V2, bed_original$V3, bed_original$V4, bed_original$V5, bed_original$V6, sam_original$V2, sam_original$V6, sam_original$V10, sam_original$V11)

write.table(new_bed, file=new_file_name, row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")