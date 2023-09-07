library(data.table)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)


inBED <- args[1]
outBED <- args[2]




DT <- data.table(read.table(inBED, sep = '\t', header=FALSE))

# Get min coord for each transcript
minsDT <- setDT(DT)[ , .SD[which.min(V2)], by = V4]   # Min of groups
maxsDT <- setDT(DT)[ , .SD[which.max(V3)], by = V4]   # Min of groups

outDT <- data.table(minsDT$V1, minsDT$V2, maxsDT$V3, minsDT$V4, minsDT$V5, minsDT$V6)

# Write files
write.table(outDT, file=outBED, row.names=FALSE, col.names=FALSE, sep=("\t"), quote=FALSE)


