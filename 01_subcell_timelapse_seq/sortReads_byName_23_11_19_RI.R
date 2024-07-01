#!/n/app/R/4.3.1-gcc-9.2.0/bin/Rscript

# set file as first argument
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
path <- args[2]

# new file to be saved
new_filename = paste0(path, '/', dat_file, "_sort.bed")

# load current data
table <- read.csv(paste0(path, '/', dat_file), sep="\t", header=FALSE)

# sort
table_new <- table[order(table$V4),]

# save
write.table(table_new, file=new_filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
