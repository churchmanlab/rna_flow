#!/n/app/R/3.5.1/bin/Rscript

# use: ./top1000genes.R *.tsv *ptc.genes.txt timepoint

# load libs
library(dplyr, lib.loc="/n/groups/churchman/bms36/programs/3.5/")

# set GS output and list of protein coding genes and gtf and sample name as arguments 1 and 2 and 3 and 4, respectively
args <- commandArgs(trailingOnly = TRUE)
GS_out <- args[1]
ptc_list <- args[2]
gtf <- args[3]
sample <- args[4]
timepoint <- args[5]

# get correct column number from timepoint
if (timepoint==1) { 
  column <- 10
} else if (timepoint==2) {
  column <- 14
} else if  (timepoint==3) {
  column <- 18
} else if  (timepoint==4) {
  column <- 22
} else if  (timepoint==5) {
  column <- 26
} else {
  column <- 0
}

# set new file name
new_file_name1 <- paste(sample, "_top1000genes_turnover.bed", sep="")
new_file_name2 <- paste(sample, "_top1000genes.turnover.MAPs.txt", sep="")

# load GS output
complete <- read.delim(GS_out, header=TRUE, sep = "\t")

# load list of ptc genes
ptc <- read.delim(ptc_list, header=TRUE, sep = "\t")

# keep only protein coding genes
complete_ptc <- inner_join(complete, ptc, by='Gene')

# identify genes that are most labeled:

# First remove genes with less than 100 reads
background_filtered <- complete_ptc %>% filter(complete_ptc[,3] >= 100 & complete_ptc[,4] >= 100 & complete_ptc[,5] >= 100 & complete_ptc[,6] >= 100 & complete_ptc[,7] >= 100)

# Next remove genes with background labeling >1%
background_filtered <- background_filtered %>% filter(background_filtered[,10] < 0.01)

# Sort genes according to timepoint MAPs and keep top 100
out <- background_filtered %>% arrange(desc(background_filtered[,column]))
out <- out[1:1000,]

# Retain only gene and timepoint MAP values
out <- out[,c(1,2,column)]

# now get bed file containing all exons 
gtf_file <- read.delim(gtf, header=FALSE, sep = "\t")

# only keep genes that are in top 100
gtf_genes <- inner_join(out, gtf_file, by=c('Gene'='V1'))
gtf_genes <- gtf_genes %>% filter(V4=='exon')

# make fake score column
gtf_genes$Score <- 0

# make bed file
bed_out <- gtf_genes %>% select(V2, V5, V6, Gene, Score, V8)

# save outputs
write.table(bed_out, new_file_name1, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(out, new_file_name2, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")


