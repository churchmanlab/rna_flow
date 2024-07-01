#!/n/app/R/3.5.1/bin/Rscript

# set file as first argument
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]

# new file to be saved
new_filename = paste0(dat_file, "_sort_for_noSingles.bed")

# load current data
table <- read.csv(dat_file, sep="\t", header=FALSE)

# sort
table_new <- table[order(table$V4),]

# save
write.table(table_new, file=new_filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
