#!/n/app/R/3.5.1/bin/Rscript

###
# After running mismatch scripts to get reads_MM.fragments.sort.bed ** full fragment sequences **, make shortened file with run_AWKforTCperT_fragments

### USE:    1. Run mismatch scripts to get reads_MM.fragments.bed *with full fragment 		
###			   sequences*
###			2. Use run_AWKforTCperT_fragments.sh to make files with select column and only 
###			   reads of interest (e.g. MTall, MTnorRNA)
###			3. Make directory 'MMfrequency' on same level as directories that contain 
###			   reads_MM.fragments.bed
###			4. sbatch -p short -t 0-01:00 
###				--wrap="../Scripts/MismatchFrequencyAveTsAndNumReadsWithMM.R"


library(data.table, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5/")
library(stringr, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5/")
library(purrr, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5/")
library(rlist, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5/")


########### Modify this section for each experiment ##########
##############################################################
Exp = 'T_cyto' # 'TL3' mouse_H TL5
args <- commandArgs(trailingOnly = TRUE)
MapMethod <- 'NA' 
reads <- 'top1000genes_turnover' # 'MTall' All MTnorRNA
path = '/n/groups/churchman/bms36/2021-05-16_T_U/'  # 
samples=c('T1_cyto', 'T2_cyto', 'T3_cyto', 'T4_cyto', 'T5_cyto')
# samples=c('H1', 'H2', 'H3', 'H4', 'H5')
#samples=paste0(Exp '_', c('0m', '7m', '15m', '30m','45m', '60m', '90m', '120m', '240m'))

# e.g. files are in path, samples[i],'_', MapMethod,'/frag_MMfrequency_',reads,'.txt')

##############################################################
##############################################################




NumSamps = length(samples)


# shortnames = c('AriTC0','AriTC1', 'AriTC2','AriTC3','AriTC4','AriTC5')

# names = c()
# for (i in c(1:NumSamps)){
# names = c(names, paste0('name',i))
# }
# samples = c()
# for (i in c(1:NumSamps)){
# samples = c(samples, get(names[i]))
# }

shortnames = samples

# Get data, make into DTs
for (i in c(1:NumSamps)) {
  # get reads
  assign(paste0('mDT',i), data.table(read.table(paste0(path, 'STAR_2023/', samples[i],'/',samples[i],'_MM_temp_turnover/',samples[i],'_frag_MMfrequency_All.txt'), header=TRUE, sep='\t', quote='',stringsAsFactors = FALSE)))
  
  
  # Make new columns with number of Ts and the conversion rate
  get(paste0('mDT',i))[, Tcount := ifelse(strand == '+', as.numeric(lapply(str_count(seq, pattern='T'),'[[',1)) + TC_mismatches, as.numeric(lapply(str_count(seq, pattern='A'),'[[',1)) + TC_mismatches)]
  get(paste0('mDT',i))[, ConvRate := round(TC_mismatches/Tcount, digits = 2)]
  
  
  
  # Get average and std dev of T count
  assign(paste0('meanTcount',i), format(mean(get(paste0('mDT',i))$Tcount, na.rm=TRUE), digits=3))
  assign(paste0('stddev',i), format(sd(get(paste0('mDT',i))$Tcount, na.rm=TRUE), digits=3))
  
  # Get number of reads with x TC conv (use 30 as arbitrary max)
  for (j in c(0:30)){
    assign(paste0(paste0('TC',j,'_'),i),nrow(get(paste0('mDT',i))[get(paste0('mDT',i))$TC_mismatches == j]))
  }
  
  # Get total number of reads 
  assign(paste0('tot',i),nrow(get(paste0('mDT',i))))
  
  # Get fragment lengths and T counts
  assign(paste0('lengths', i),nchar(get(paste0('mDT',i))$seq))
  assign(paste0('TcountsDist', i),get(paste0('mDT',i))$Tcount)
}



# Get 2d joint frequency distribution for T counts and TC conversions [n,k]
for (i in c(1:NumSamps)){
  DT <- get(paste0('mDT',i))
  JointFreq <- table(DT$Tcount,DT$TC_mismatches)
  write.table(JointFreq, file=paste0(path, 'STAR_2023/', samples[i],'/',samples[i],'_MM_temp_turnover/',samples[i],'_',reads,'_TcountANDTCconv.txt'), sep=("\t"), quote=FALSE, col.names=NA)
}



# Make plots for fragment length and T counts
samples = samples # shortnames
pdf(paste0(path, 'n_k_data_2023/', Exp, "_", reads,'_LengthAndTcountDist.pdf'), width = 4, height = 8)
par(mfrow=c(2,1), cex.lab = 1)

for (i in c(1:NumSamps)){
  hist(get(paste0('lengths',i)), main= paste0(samples[i], ' fragments length distribution'), xlab = paste0(samples[i], ' lengths'), breaks = 50)
  assign(paste0('Thist',i), hist(get(paste0('TcountsDist',i)), main= paste0(samples[i], ' T counts distribution'), xlab = paste0('Number of Ts'), breaks = 135))
}
# measure skewness and kurtosis
# library(moments)
# assign(paste0('skew',i), skewness(get(paste0('mDT',i))$Tcount))
# assign(paste0('kurt',i), kurtosis(get(paste0('mDT',i))$Tcount))
dev.off()



# Make table with counts for distribution of T counts and 
for (i in c(1:NumSamps)) { # Make the new table names 
  assign(paste0('TcountDT',i), data.table(table(get(paste0('mDT',i))$Tcount)))
}
for (i in c(1:NumSamps)) { # Name columns in each table
  d=get(paste0('TcountDT',i))
  colnames(d)[1]='Tcount'
  colnames(d)[2]=samples[i]
  d[, Tcount := as.numeric(Tcount)] # Make Tcount column numeric
  assign(paste0('TcountDT',i),d)
}
# Find max T count
Tcounts = c()
for (i in c(1:NumSamps)) {
  Tcounts = c(Tcounts, get(paste0('TcountDT',i))$Tcount)
}
maxTcount = max(Tcounts)
FillTable = data.table(Tcount=seq(0,maxTcount))

# Fix tables to include all values of T counts
for (i in c(1:NumSamps)) { # Name columns in each table
  d=get(paste0('TcountDT',i))
  d=merge(d,FillTable, all.y=TRUE)
  # d[,x := NULL]
  # d[is.na(d)] <- 0
  assign(paste0('TcountDT',i),d)
}

mylist=list()
for (i in c(1:NumSamps)) { # Make list of data frames
  mylist[[i]] <- get(paste0('TcountDT',i))
}
TcountDT <- Reduce(merge, mylist)
setorder(TcountDT)

write.table(TcountDT, file=paste0(path, 'n_k_data_2023/', Exp, '_frag_', reads, '_TcountDistributions.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)





# Plot average # Ts per read with std dev and num with TC MM

# Get max number of TC conversions in data
TC_mm = c()
for (i in c(1:NumSamps)) {
  TC_mm = c(TC_mm, get(paste0('mDT',i))$TC_mismatches)
}
maxTC = max(TC_mm)

pdf(paste0(path, 'n_k_data_2023/', Exp, "_", reads,'_AveTs_NumReadsWithMM.pdf'), width = 9, height = 10)

par(mfrow=c(3,1), cex.lab = 1)

# For barplot T counts rates

Tcountslist=c()
sdevlist=c()
totslist=c()
for (i in c(1:NumSamps)){
  Tcountslist=c(Tcountslist, get(paste0('meanTcount',i)))
  sdevlist=c(sdevlist, get(paste0('stddev',i)))
  totslist=c(totslist, get(paste0('tot',i)))
}
Tcounts <- as.numeric(Tcountslist) 
sdevs <- as.numeric(sdevlist) 
tots <-  as.numeric(totslist)
minReadCounts = min(totslist[1:(length(totslist)-1)])

ylimitsTcount = c(0,55)

xx = barplot(Tcounts, beside=TRUE, names.arg=samples, col=c('orange'), cex.names = .8, ylab = 'Mean number T in fragment (std dev)', main = 'T counts per fragment', ylim=ylimitsTcount) # 

# Add std dev as text to top of bars
text(x = xx, y = Tcounts, label = sdevs, pos = 3, cex = 0.8, col = "black")
# Add total number of reads as text to bottom of bars
text(x = xx, y = 1, label = tots, pos = 3, cex = 0.8, col = "black")
# Add mean value near top of bar
text(x = xx, y = Tcounts-5, label = Tcounts, pos = 3, cex = 0.8, col = "black")


for (j in c(0:maxTC)){
  assign(paste0('TClist',j),c())
}
for (j in c(0:maxTC)){
  for (i in c(1:NumSamps)){
    assign(paste0('TClist',j),c(get(paste0('TClist',j)), get(paste0(paste0('TC',j,'_'),i))))
  }
}
# For frequency of each count
for (j in c(0:maxTC)){
  assign(paste0('TC',j,'s'),as.numeric(get(paste0('TClist',j))))
}

vec = c()
for (j in c(0:maxTC)){
  vec = c(vec,get(paste0('TC',j,'s')))
}
TCs <- matrix(vec, nrow=NumSamps, ncol=maxTC+1)

# Colors
if (NumSamps > 7) {
  cols = c('dodgerblue','darkblue', 'blue','darkorchid','violet', 'forestgreen','green', 'yellowgreen', 'yellow','gold1', 'orange', 'red', 'red3','brown')} # skyblue
if (NumSamps == 6 | NumSamps == 7) {
  cols = c('dodgerblue','violet', 'forestgreen', 'gold1', 'orange', 'red', 'brown')}
if (NumSamps == 9) {
  cols = c('dodgerblue','darkblue','violet', 'forestgreen','green', 'yellow', 'orange', 'red', 'brown')} # skyblue
if (NumSamps <8) {
  cols = c('dodgerblue','violet', 'forestgreen', 'gold1', 'orange', 'red', 'brown')}
if (NumSamps == 5) {
  cols = c('dodgerblue', 'forestgreen', 'gold1', 'orange', 'red')}
if (NumSamps == 4) {
  cols = c('dodgerblue', 'green', 'gold1', 'red')}
if (Exp == 'TL5' & NumSamps == 8) {
  cols = c('dodgerblue', 'green', 'gold1', 'red','dodgerblue', 'green', 'gold1', 'red')} # skyblue

cols = cols[1:NumSamps]


ylimitsBar1 = c(0, 1.2*max(TC1s))


# Plot number of fragments with n TC conversions
names = as.character(c(0:maxTC))
xxx = barplot(TCs, beside=TRUE, names.arg=names, col=cols, cex.names = 1, ylab = 'Frequency', main = 'Number of fragments with n TC conversions', legend.text=samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8),ylim=ylimitsBar1) 

text(x = xxx+.3, y = TCs+(max(ylimitsBar1/50)), label = TCs, pos = 3, cex = 0.6, srt=90,col = 'black')


# Now plot same as above but normalized to read counts in sample
# Make matrix with readcounts to divide TCs matrix
readcountMat = matrix(tots, nrow=NumSamps, ncol=maxTC+1)
NormTCs = round(TCs/readcountMat*minReadCounts, 0)

ylimitsBarNorm = c(0, 1.1*max(NormTCs[,2]))

xxx = barplot(NormTCs, beside=TRUE, names.arg=names, col=cols, cex.names = 1, ylab = 'Frequency', main = paste0('Number of fragments with n TC conversions\n normalized per ',minReadCounts,' reads'), legend.text=samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8),ylim=ylimitsBarNorm) 

text(x = xxx+.3, y = NormTCs+(max(ylimitsBarNorm)/50), label = NormTCs, pos = 3, cex = 0.6, srt=90,col = 'black')

dev.off()

