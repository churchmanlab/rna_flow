#! /usr/bin/Rscript
# 20240211
# Robert Ietswaart
# Source: 2022-12-27_kd_ii_analysis_v_2024.R by Brendan Smalec

# KDii unspliced junction analysis

# set output directory
dir_out <- "/n/groups/churchman/ri23/bseq/GS20240122_KD_ii/figures/"

# load libs
library(dplyr) # dplyr_1.1.3
library(ggplot2) # ggplot2_3.4.4
library(tidyr) # tidyr_1.3.0
# library(ggpubr) # ggpubr_0.6.0

base_size <- 6

GS_path_exon <- "/n/groups/churchman/ri23/bseq/GS20240112_KD_ii/" #GS20230110_RBPii #GS20240112_KD_ii
GS_path_nascentRNA <- "/n/groups/churchman/ri23/bseq/GS20240122_KD_ii/"

wt_GS_path <- "/n/groups/churchman/ri23/bseq/GS20231201_K562/"

# load data
# human_PTC_genes <- read.delim("/Users/bsmalec/Dropbox (HMS)/PhD/Sequencing_data/genomes/edited_genomes/hg38/human_ptc_list.txt", header=TRUE, sep="\t")
rep_prefixes <- c("rep1", "rep2")
compartments <- c("tot", "nuc")
sample_prefix <- c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1","scramble")
read_types <- c('retained_introns','retained_introns_with_exons', 'no_retained_introns', 'protein_coding', 'exons_RI', 'exons_BMS')  #c('unspliced_junctions','introns', 'not_introns','exons')  

# load files
for (sample in sample_prefix) {
    for (rep_prefix in rep_prefixes) {
        for (comp in compartments) {
            for (rt in read_types){
                # if (rt == 'exons'){
                #   filename = paste0(GS_path_exon, "GS20240112_", sample, "_", comp, "_", rep_prefix, ".tsv")
                # }
                if (rt == 'exons_BMS'){
                  filename = paste0(GS_path_exon, "GS20240112_", sample, "_", comp, "_", rep_prefix, ".tsv")
                }
                else if (rt == 'exons_RI'){
                  filename = paste0(GS_path_nascentRNA, "GS20240122_", sample, "_", comp, "_", rep_prefix, "_default.tsv")
                }
                else{
                  filename = paste0(GS_path_nascentRNA, "GS20240122_", sample, "_", comp, "_", rep_prefix, '_',rt, ".tsv")
                }
                assign( paste0(sample, "_", comp, "_", rep_prefix, '_', rt), read.delim(filename, header=TRUE, sep = "\t"))
            }#for read_type
        }#for compartments
    }#for reps
}#for samples


# wt
# wt_tot_rep1 <- read.delim(paste0(wt_GS_path, "GS20231201_T_tot.tsv"), header=TRUE, sep = "\t")
# wt_tot_rep2 <- read.delim(paste0(wt_GS_path, "GS20231201_U_tot.tsv"), header=TRUE, sep = "\t")
# wt_nuc_rep1 <- read.delim(paste0(wt_GS_path, "GS20231201_T_nuc.tsv"), header=TRUE, sep = "\t")
# wt_nuc_rep2 <- read.delim(paste0(wt_GS_path, "GS20231201_U_nuc.tsv"), header=TRUE, sep = "\t")

# load bayesian model data - may need to update until finalized
human_rates <- read.delim("/n/groups/churchman/ri23/bseq/Bayes20240120_K562/Bayes_Rates_20240120_human_final.tsv", header=T)
# mouse_rates <- read.delim("/n/groups/churchman/ri23/bseq/Bayes20240120_3T3/Bayes_Rates_20240120_mouse_final.tsv", header=T)

# load bayes factor data
bayes_factor_human <- read.delim("/n/groups/churchman/ri23/bseq/BayesFactor20240112/Bayes_factor_20240112_human_final.tsv", header=T, sep='\t')
# bayes_factor_mouse <- read.delim("/n/groups/churchman/ri23/bseq/BayesFactor20240112/Bayes_factor_20240112_mouse_final.tsv", header=T, sep='\t')

# get human genes where bayes factor >100 in both replicates
human_bayes_list <- bayes_factor_human %>% select(Gene, Symbol, PUND) 
human_bayes_list_nd <- human_bayes_list %>% filter(PUND=='True') %>% mutate(nuc_deg='Yes')
human_bayes_list_no_nd <- human_bayes_list %>% filter(PUND=='False') %>% mutate(nuc_deg='No')
# human_bayes_list_no_nd <- anti_join(human_bayes_list, human_bayes_list_nd, by='Gene') %>% mutate(nuc_deg='No') #PUND==NA gets labeled no here
human_bayes_list <- rbind(human_bayes_list_nd, human_bayes_list_no_nd)
human_rates_i <- inner_join(human_rates, human_bayes_list, by='Gene') %>% select(-Symbol.y) %>% rename(Symbol=Symbol.x)

# #TEMP Detained intron list instead of PUNDs
# DI33 <- read.delim("/n/groups/churchman/ri23/bseq/Comparison20240122/GenesDev2015Sharp_DIgenes_33.csv", header=F, sep=',',col.names=c('Symbol'))
# DI64 <- read.delim("/n/groups/churchman/ri23/bseq/Comparison20240122/GenesDev2015Sharp_DIgenes_64.csv", header=F, sep=',',col.names=c('Symbol'))
# DI33 <- inner_join(DI33, bayes_factor_human %>% select(Gene, Symbol),by='Symbol')
# DI64 <- inner_join(DI64, bayes_factor_human %>% select(Gene, Symbol),by='Symbol')
# DI33$nuc_deg <- 'Yes'
# DI64$nuc_deg <- 'Yes'
# '%!in%' <- function(x,y)!('%in%'(x,y))
# human_symbol_list <- bayes_factor_human %>% select(Symbol,Gene)
# # human_symbol_list_no_DI <- human_symbol_list[human_symbol_list$Gene%!in% DI33$Gene,]
# human_symbol_list_no_DI <- human_symbol_list[human_symbol_list$Gene%!in% DI64$Gene,]
# human_symbol_list_no_DI$nuc_deg <- "No"
# # human_bayes_list <- rbind(DI33, human_symbol_list_no_DI)
# human_bayes_list <- rbind(DI64, human_symbol_list_no_DI)
# human_rates_i <- inner_join(human_rates, human_bayes_list, by='Gene') %>% select(-Symbol.y) %>% rename(Symbol=Symbol.x)

# blank table
gg_all <- data.frame(matrix(ncol=8))
colnames(gg_all) <- c("Gene","Symbol","MAP","Sample","Time","RNA","Rep","Read_type") #c("Gene","Symbol","MAP","Sample","Time","RNA","Rep")

for (rt in read_types){

  # total RNAs - everything except wt
  rep_prefixes <- c("rep1", "rep2")
  sample_prefix <- c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1","scramble")
  
  for (sample in sample_prefix) {
    for (rep_prefix in rep_prefixes) {
        temp <- get(paste0(sample, "_tot_", rep_prefix, '_', rt))
        temp2 <- temp[,c(1,2,5,9)]
        temp3 <- temp2 %>% mutate(Sample=sample, Time1=0, Time2=60, RNA='tot', Rep=rep_prefix, Read_type=rt)
        colnames(temp3) <- c("Gene", "Symbol", "MAP", "MAP", "Sample", "Time", "Time", "RNA", "Rep", "Read_type")
        temp4 <- rbind(temp3[,c(1,2,3,5,6,8,9,10)], temp3[,c(1,2,4,5,7,8,9,10)])
        gg_all <- rbind(gg_all, temp4)
      }
  }
  
  # nuc RNAs - everything except wt
  for (sample in sample_prefix) {
    for (rep_prefix in rep_prefixes) {
      temp <- get(paste0(sample, "_nuc_", rep_prefix, '_', rt))
      temp2 <- temp[,c(1,2,5)]
      temp3 <- temp2 %>% mutate(Sample=sample, Time=60, RNA='nuc', Rep=rep_prefix, Read_type=rt)
      colnames(temp3) <- c("Gene", "Symbol", "MAP", "Sample", "Time", "RNA", "Rep", "Read_type")
      gg_all <- rbind(gg_all, temp3)
    }
  }

}#for read_types

# # total RNAs - wt
# rep_prefixes <- c("rep1", "rep2")
# sample_prefix <- c("wt")
# 
# for (sample in sample_prefix) {
#   for (rep_prefix in rep_prefixes) {
#     temp <- get(paste0(sample, "_tot_", rep_prefix))
#     temp2 <- temp[,c(1,2,5,17)]
#     temp3 <- temp2 %>% mutate(Sample=sample, Time1=0, Time2=60, RNA='tot', Rep=rep_prefix)
#     colnames(temp3) <- c("Gene", "Symbol", "MAP", "MAP", "Sample", "Time", "Time", "RNA", "Rep")
#     temp4 <- rbind(temp3[,c(1,2,3,5,6,8,9)], temp3[,c(1,2,4,5,7,8,9)])
#     gg_all <- rbind(gg_all, temp4)
#   }
# }
# 
# # nuc RNAs - wt
# for (sample in sample_prefix) {
#   for (rep_prefix in rep_prefixes) {
#     temp <- get(paste0(sample, "_nuc_", rep_prefix))
#     temp2 <- temp[,c(1,2,17)]
#     temp3 <- temp2 %>% mutate(Sample=sample, Time=60, RNA='nuc', Rep=rep_prefix)
#     colnames(temp3) <- c("Gene", "Symbol", "MAP", "Sample", "Time", "RNA", "Rep")
#     gg_all <- rbind(gg_all, temp3)
#   }
# }

gg_all <- gg_all %>% filter(Time>-1)
# write.table(gg_all, file=paste0(dir_out,"kd_ii_gs_MAP_nascentRNA.txt"), quote=F, sep="\t", row.names = F)


ggplot(gg_all, aes(x=Rep, y=MAP)) +
  geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
  geom_boxplot(aes(fill=Read_type), outlier.shape = NA, lwd=0.1) +
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylab("MAP fraction new RNA") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Replicate") +
  ylim(c(0,0.5))+
  guides(fill=guide_legend(title="Read_type"))
ggsave(paste0(dir_out,"MAPs_read_type_RI.pdf"), plot=last_plot(), width=7, height=5, units="in", dpi=300)


# compare kd MAP / scramble MAP 
pseudo_count <- 0.000001

# blank table
gg_compare1 <- data.frame(matrix(ncol=10))
colnames(gg_compare1) <- c("Gene","Symbol.x","kd_over_scr","kd_minus_scr","kd_percent_change","robert","Sample","RNA","Rep","Read_type")

for (rt in read_types){
  # compare kd MAP / scramble MAP - tot
  rep_prefixes <- c("rep1", "rep2")
  sample_prefix <- c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1")
  for (sample in sample_prefix) {
    for (rep_prefix in rep_prefixes) {
      temp1 <- get(paste0(sample, "_tot_", rep_prefix,'_', rt))
      temp2 <- get(paste0("scramble_tot_", rep_prefix,'_', rt))
      temp3 <- inner_join(temp1, temp2, by='Gene')
      temp3$kd_over_scr <- (temp3[,9] + pseudo_count)/(temp3[,18] + pseudo_count) #default
      temp3$kd_minus_scr <- temp3[,9] - temp3[,18]
      temp3$kd_percent_change <- 100*(temp3[,9] - temp3[,18])/(temp3[,18] + pseudo_count)
      temp3$robert <- (1 - temp3[,9] + pseudo_count)/(1 - temp3[,18] + pseudo_count)
      temp3 <- temp3 %>% mutate(Sample=sample, RNA="tot",  Rep=rep_prefix, Read_type=rt)
      temp4 <- temp3 %>% select(Gene, Symbol.x, kd_over_scr, kd_minus_scr, kd_percent_change, robert, Sample, RNA, Rep, Read_type)
      assign( paste0(sample, "_", rep_prefix, "_compare1"), temp3)
      gg_compare1 <- rbind(gg_compare1, temp4)
    }
  }
  
  # # wt
  # sample <- "wt"
  # for (rep_prefix in rep_prefixes) {
  #   temp1 <- get(paste0("wt_tot_", rep_prefix))
  #   temp2 <- get(paste0("scramble_tot_", rep_prefix))
  #   temp3 <- inner_join(temp1, temp2, by='Gene')
  #   temp3$kd_over_scr <- temp3[,17]/temp3[,30]
  #   temp3$kd_minus_scr <- temp3[,17]-temp3[,30]
  #   temp3$kd_percent_change <- abs(temp3[,30]-temp3[,17])/temp3[,30]
  #   temp3$robert <- (1-temp3[,17])/(1-temp3[,30])
  #   temp3 <- temp3 %>% mutate(Sample=sample, RNA="tot",  Rep=rep_prefix)
  #   temp4 <- temp3 %>% select(Gene, Symbol.x, kd_over_scr, kd_minus_scr, kd_percent_change, robert, Sample, RNA, Rep)
  #   assign( paste0(sample, "_", rep_prefix, "_compare1"), temp3)
  #   gg_compare1 <- rbind(gg_compare1, temp4)
  # }
  
  # nuc
  for (sample in sample_prefix) {
    for (rep_prefix in rep_prefixes) {
      temp1 <- get(paste0(sample, "_nuc_", rep_prefix, '_', rt))
      temp2 <- get(paste0("scramble_nuc_", rep_prefix, '_', rt))
      temp3 <- inner_join(temp1, temp2, by='Gene')  ##default
      temp3$kd_over_scr <- (temp3[,5] + pseudo_count)/(temp3[,10] + pseudo_count)
      temp3$kd_minus_scr <- temp3[,5]-temp3[,10]
      temp3$kd_percent_change <- 100*(temp3[,5]-temp3[,10])/(temp3[,10] + pseudo_count)
      temp3$robert <- (1 - temp3[,5] + pseudo_count)/(1 - temp3[,10] + pseudo_count)
      temp3 <- temp3 %>% mutate(Sample=sample, RNA="nuc", Rep=rep_prefix, Read_type=rt)
      temp4 <- temp3 %>% select(Gene, Symbol.x, kd_over_scr, kd_minus_scr, kd_percent_change, robert, Sample, RNA, Rep, Read_type)
      assign( paste0(sample, "_", rep_prefix, "_compare1"), temp3)
      gg_compare1 <- rbind(gg_compare1, temp4)
    }
  }
  
  # #wt
  # sample <- "wt"
  # for (rep_prefix in rep_prefixes) {
  #   temp1 <- get(paste0("wt_nuc_", rep_prefix))
  #   temp2 <- get(paste0("scramble_nuc_", rep_prefix))
  #   temp3 <- inner_join(temp1, temp2, by='Gene')
  #   temp3$kd_over_scr <- temp3[,17]/temp3[,26]
  #   temp3$kd_minus_scr <- temp3[,17]-temp3[,26]
  #   temp3$kd_percent_change <- abs(temp3[,26]-temp3[,17])/temp3[,26]
  #   temp3$robert <- (1-temp3[,17])/(1-temp3[,26])
  #   temp3 <- temp3 %>% mutate(Sample=sample, RNA="nuc", Rep=rep_prefix)
  #   temp4 <- temp3 %>% select(Gene, Symbol.x, kd_over_scr, kd_minus_scr, kd_percent_change, robert, Sample, RNA, Rep)
  #   assign( paste0(sample, "_", rep_prefix, "_compare1"), temp3)
  #   gg_compare1 <- rbind(gg_compare1, temp4)
  # }
}#for read_types

gg_compare1 <- gg_compare1 %>% filter(!is.na(Gene))

gg_compare1_bayes <- inner_join(gg_compare1, human_bayes_list, by='Gene')

gg_compare1_bayes$Sample <- factor(gg_compare1_bayes$Sample, levels=c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1","wt"))
gg_compare1_bayes$Read_type <- factor(gg_compare1_bayes$Read_type, levels=read_types)

# write.table(gg_compare1_bayes, file=paste0(dir_out,"kd_ii_calcs_vs_scramble_nascentRNA.txt"), quote=F, sep="\t", row.names = F)

# rename the compartments (currently in shorthand)
comp_real_names <- data.frame(matrix(c("tot", "Total RNA", 
                                       "nuc", "Nuclear RNA"), ncol=2, byrow=T))

# rename the replicates (currently in shorthand)
rep_real_names <- data.frame(matrix(c("rep1", "Rep 1", 
                                      "rep2", "Rep 2"), ncol=2, byrow=T))

gg_compare1_bayes_new <- full_join(gg_compare1_bayes, comp_real_names, by=c('RNA'='X1')) %>% select(-RNA) %>% rename(RNA=X2) %>% full_join(., rep_real_names, by=c('Rep'='X1')) %>% select(-Rep) %>% rename(Rep=X2)


# compare_means(kd_minus_scr ~ nuc_deg, data=gg_compare1_bayes_new, method='wilcox.test', group.by=c("Rep","Sample","RNA"), p.adjust.method = "fdr")

# ggplot(gg_compare1_bayes_new, aes(x=Rep, y=kd_over_scr)) +
#   geom_boxplot(aes(fill=Read_type), outlier.shape = NA)+#outlier.alpha = 1) +
#   facet_grid(cols=vars(Sample), rows=vars(RNA)) +
#   theme_bw(base_size = base_size) +
#   ylim(c(0,3)) #+
# # ggsave(paste0(dir_out,"MAPs_fc_all_samples.pdf"), plot=last_plot(), width=3.5, height=2.5, units="in", dpi=300)

# ggplot(gg_compare1_bayes_new, aes(x=Rep, y=kd_over_scr)) +
#   geom_boxplot(aes(fill=interaction(Read_type, nuc_deg), group=interaction(Read_type, nuc_deg)), outlier.shape = NA) +
#   facet_grid(cols=vars(Sample), rows=vars(RNA)) +
#   theme_bw(base_size = base_size) +
#   ylim(c(0,3)) #+

ggplot(gg_compare1_bayes_new, aes(x=Rep, y=kd_minus_scr)) +
  geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
  geom_boxplot(aes(fill=interaction(nuc_deg,Read_type), group=interaction(nuc_deg,Read_type, Rep)), outlier.shape = NA, lwd=0.1) +
  # geom_boxplot(aes(fill=Read_type), outlier.shape = NA, lwd=0.1) + #outlier.alpha = 0
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylab("Change in fraction new RNA (KD - scr)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Replicate") +
  ylim(c(-0.3,0.5))+#c(-0.5,1)) +
  guides(fill=guide_legend(title="PUND.Read_type"))
# ggsave(paste0(dir_out,"MAPs_delta_read_type.nuc_deg.pdf"), plot=last_plot(), width=7, height=5, units="in", dpi=300)
ggsave(paste0(dir_out,"MAPs_delta_read_type.nuc_deg_RI.pdf"), plot=last_plot(), width=7, height=5, units="in", dpi=300)

ggplot(gg_compare1_bayes_new, aes(x=Rep, y=kd_minus_scr)) +
  geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
  geom_boxplot(aes(fill=Read_type), outlier.shape = NA, lwd=0.1) +
  # geom_boxplot(aes(fill=Read_type), outlier.shape = NA, lwd=0.1) + #outlier.alpha = 0
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylab("Change in fraction new RNA (KD - scr)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Replicate") +
  ylim(c(-0.25,0.35))+# ylim(c(-0.5,1)) +
  guides(fill=guide_legend(title="Read_type"))
# ggsave(paste0(dir_out,"MAPs_delta_read_type.pdf"), plot=last_plot(), width=7, height=5, units="in", dpi=300)
ggsave(paste0(dir_out,"MAPs_delta_read_type_RI.pdf"), plot=last_plot(), width=7, height=5, units="in", dpi=300)

ggplot(gg_compare1_bayes_new, aes(x=Rep, y=log(kd_over_scr))) +
  geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
  geom_boxplot(aes(fill=Read_type), outlier.shape = NA, lwd=0.1) + #outlier.alpha = 0
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylab("Log Fold change in fraction new RNA (kd / scr)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Replicate") +
  ylim(c(-2,3)) #+
  # guides(fill=guide_legend(title="PUND\nGene"))
# ggsave(paste0(dir_out,"MAPs_logfc_all_samples.pdf"), plot=last_plot(), width=3.5, height=2.5, units="in", dpi=300)

#old
# ggplot(gg_compare1_bayes_new, aes(x=Rep, y=kd_minus_scr)) +
#   geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
#   geom_boxplot(aes(fill=nuc_deg), outlier.alpha = 0, lwd=0.2) +
#   facet_grid(cols=vars(Sample), rows=vars(RNA)) +
#   theme_bw(base_size = base_size) +
#   ylim(c(-0.75,0.25)) +
#   ylab("Change in fraction new RNA") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         strip.background=element_rect(fill="white"),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   xlab("Replicate") +
#   guides(fill=guide_legend(title="PUND\nGene"))
# ggsave(paste0(dir_out,"MAPs_delta_all_samples.pdf"), plot=last_plot(), width=4, height=3, units="in", dpi=300)

####################################################################################### For fig 2
temp_nuc <- gg_compare1_bayes_new %>% filter(RNA=='Nuclear RNA')
temp_nuc <- temp_nuc %>% unite(new_row, c("Sample","Rep"))
temp_nuc$nuc_deg <- factor(temp_nuc$nuc_deg, levels=c("Yes","No"))
temp_nuc$new_row <- factor(temp_nuc$new_row, levels=c("wt_Rep 2", "wt_Rep 1", "ZFC3H1_Rep 2", "ZFC3H1_Rep 1", "PABPN1_Rep 2", "PABPN1_Rep 1", "EXOSC10_Rep 2", "EXOSC10_Rep 1", "DIS3_Rep 2","DIS3_Rep 1"))

ggplot(temp_nuc, aes(x=new_row, y=kd_minus_scr)) +
  geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
  geom_boxplot(aes(fill=Read_type), outlier.alpha = 0, lwd=0.2, ) +
  theme_bw(base_size = base_size) +
  ylab("Change in fraction new nuclear RNA\n(KD - scr)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=0),
        legend.position = "none") +
  xlab("Replicate") +
  guides(fill=guide_legend(title="PUND\nGene")) +
  scale_y_continuous(position = "right") +
  # scale_fill_manual(values=c('red3','grey50')) + 
  coord_flip() 
# ggsave(paste0(dir_out,"MAPs_delta_all_samples_fig.pdf"), plot=last_plot(), width=2.25, height=2.5, units="in", dpi=300)


# just rep 1
temp_nuc <- gg_compare1_bayes_new %>% filter(RNA=='Nuclear RNA')
temp_nuc2 <- temp_nuc %>% filter(Rep=="Rep 1")
temp_nuc2$nuc_deg <- factor(temp_nuc2$nuc_deg, levels=c("Yes","No"))
temp_nuc2$Sample <- factor(temp_nuc2$Sample, levels=c("wt", "ZFC3H1", "PABPN1", "EXOSC10", "DIS3"))

#MS figure, NOT UPDATED YET
# ggplot(temp_nuc2, aes(x=Sample, y=kd_minus_scr)) +
#   geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
#   geom_boxplot(aes(fill=nuc_deg), outlier.alpha = 0, lwd=0.2, ) +
#   theme_bw(base_size = base_size) +
#   ylab("Change in fraction new RNA") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.position = "none") +
#   xlab("Replicate") +
#   guides(fill=guide_legend(title="PUND\nGene")) +
#   scale_y_continuous(position = "right", limits=c(-0.75,0.25)) +
#   scale_fill_manual(values=c('red3','grey50')) + 
#   coord_flip() 
# # ggsave(paste0(dir_out,"MAPs_delta_all_samples_fig_rep1.pdf"), plot=last_plot(), width=2.25, height=2, units="in", dpi=300)

# just rep 2
temp_nuc <- gg_compare1_bayes_new %>% filter(RNA=='Nuclear RNA')
temp_nuc2 <- temp_nuc %>% filter(Rep=="Rep 2")
temp_nuc2$nuc_deg <- factor(temp_nuc2$nuc_deg, levels=c("Yes","No"))
temp_nuc2$Sample <- factor(temp_nuc2$Sample, levels=c("wt", "ZFC3H1", "PABPN1", "EXOSC10", "DIS3"))

#MS figure, NOT UPDATED YET
# ggplot(temp_nuc2, aes(x=Sample, y=kd_minus_scr)) +
#   geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
#   geom_boxplot(aes(fill=nuc_deg), outlier.alpha = 0, lwd=0.2, ) +
#   theme_bw(base_size = base_size) +
#   ylab("Change in fraction new RNA") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.position = "none") +
#   xlab("Replicate") +
#   guides(fill=guide_legend(title="PUND\nGene")) +
#   scale_y_continuous(position = "right", limits=c(-0.75,0.25)) +
#   scale_fill_manual(values=c('red3','grey50')) + 
#   coord_flip() 
# # ggsave(paste0(dir_out,"MAPs_delta_all_samples_fig_rep2.pdf"), plot=last_plot(), width=2.25, height=2, units="in", dpi=300)


# old
# ggplot(temp_nuc, aes(x=new_row, y=log2(robert))) +
#   geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
#   geom_boxplot(aes(fill=nuc_deg), outlier.alpha = 0, lwd=0.2, ) +
#   theme_bw(base_size = base_size) +
#   ylab("log2 (1-MAP(kd))/(1-MAP(scr))") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.position = "none") +
#   xlab("Replicate") +
#   guides(fill=guide_legend(title="PUND\nGene")) +
#   scale_fill_manual(values=c('red3','grey50')) + 
#   scale_y_continuous(position = "right", limits=c(-2,4)) +
#   coord_flip() 
# ggsave(paste0(dir_out,"MAPs_delta_all_samples_robert.pdf"), plot=last_plot(), width=2.25, height=2.5, units="in", dpi=300)

# old
# ggplot(temp_nuc, aes(x=new_row, y=log2(kd_over_scr))) +
#   geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
#   geom_boxplot(aes(fill=nuc_deg), outlier.alpha = 0, lwd=0.2, ) +
#   theme_bw(base_size = base_size) +
#   ylab("log2 (MAP(KD)/MAP(scr))") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.position = "none") +
#   xlab("Replicate") +
#   guides(fill=guide_legend(title="PUND\nGene")) +
#   scale_fill_manual(values=c('red3','grey50')) + 
#   scale_y_continuous(position = "right", limits=c(-4,2)) +
#   coord_flip() 
# ggsave(paste0(dir_out,"MAPs_delta_all_samples_fc.pdf"), plot=last_plot(), width=2.25, height=2.5, units="in", dpi=300)


# old
# temp_tot <- gg_compare1_bayes_new %>% filter(RNA=='Total RNA')
# temp_tot <- temp_tot %>% unite(new_row, c("Sample","Rep"))
# temp_tot$nuc_deg <- factor(temp_tot$nuc_deg, levels=c("Yes","No"))
# temp_tot$new_row <- factor(temp_tot$new_row, levels=c("wt_Rep 2", "wt_Rep 1", "ZFC3H1_Rep 2", "ZFC3H1_Rep 1", "PABPN1_Rep 2", "PABPN1_Rep 1", "EXOSC10_Rep 2", "EXOSC10_Rep 1", "DIS3_Rep 2","DIS3_Rep 1"))

# old
# ggplot(temp_tot, aes(x=new_row, y=kd_minus_scr)) +
#   geom_hline(yintercept = 0, lwd=0.1, color="grey60") +
#   geom_boxplot(aes(fill=nuc_deg), outlier.alpha = 0, lwd=0.2, ) +
#   theme_bw(base_size = base_size) +
#   ylab("Change in fraction new RNA") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.position = "none") +
#   xlab("Replicate") +
#   guides(fill=guide_legend(title="PUND\nGene")) +
#   scale_y_continuous(position = "right", limits=c(-0.75,0.25)) +
#   scale_fill_manual(values=c('red3','grey50')) + 
#   coord_flip() 
# ggsave(paste0(dir_out,"MAPs_delta_all_samples_fig_tot.pdf"), plot=last_plot(), width=2.25, height=2.5, units="in", dpi=300)


ggplot(gg_compare1_bayes, aes(x=Rep, y=kd_percent_change)) +
  geom_boxplot(aes(fill=Read_type), outlier.alpha = 0) +
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylim(c(-100,100))

ggplot(gg_compare1_bayes_new, aes(x=Rep, y=kd_percent_change)) +
  # geom_boxplot(aes(fill=Read_type), outlier.alpha = 0) +
  geom_boxplot(aes(fill=interaction(nuc_deg,Read_type), group=interaction(nuc_deg,Read_type, Rep)), outlier.shape = NA, lwd=0.1) +
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylim(c(-100,100))+ 
  guides(fill=guide_legend(title="PUND.Read_type"))

############################################################################################# delta MAP vs K nuc deg

# look not just at punds

# MAP this chunk 
human_nuc_deg_rep1 <- human_rates_i %>%
  select(Gene, Symbol,
         T.half_life_nucdeg.MAP, T.half_life_nucdeg.0.025.quantile, T.half_life_nucdeg.0.975.quantile, 
         T.k_nucdeg.MAP, T.k_nucdeg.0.025.quantile, T.k_nucdeg.0.975.quantile,
         T.k_nucdeg.underflow, nuc_deg) %>%
  rename(half_life_MAP=T.half_life_nucdeg.MAP, half_life_lower=T.half_life_nucdeg.0.025.quantile, half_life_upper=T.half_life_nucdeg.0.975.quantile, 
         rate_MAP=T.k_nucdeg.MAP, rate_lower= T.k_nucdeg.0.025.quantile, rate_upper=T.k_nucdeg.0.975.quantile,
         fit_chi2=T.k_nucdeg.underflow) %>%
  mutate(Species='Human', Rate='Nuclear_degradation', Time='Nuclear_degradation', Rep='Rep_1')
human_nuc_deg_rep1 <- human_nuc_deg_rep1 %>% filter(half_life_MAP>-1)

human_nuc_deg_rep2 <- human_rates_i %>% 
  select(Gene, Symbol,
         U.half_life_nucdeg.MAP, U.half_life_nucdeg.0.025.quantile, U.half_life_nucdeg.0.975.quantile, 
         U.k_nucdeg.MAP, U.k_nucdeg.0.025.quantile, U.k_nucdeg.0.975.quantile,
         U.k_nucdeg.underflow, nuc_deg) %>%
  rename(half_life_MAP=U.half_life_nucdeg.MAP, half_life_lower=U.half_life_nucdeg.0.025.quantile, half_life_upper=U.half_life_nucdeg.0.975.quantile, 
         rate_MAP=U.k_nucdeg.MAP, rate_lower= U.k_nucdeg.0.025.quantile, rate_upper=U.k_nucdeg.0.975.quantile,
         fit_chi2=U.k_nucdeg.underflow) %>%
  mutate(Species='Human', Rate='Nuclear_degradation', Time='Nuclear_degradation', Rep='Rep_2')
human_nuc_deg_rep2 <- human_nuc_deg_rep2 %>% filter(half_life_MAP>-1)



nuc_deg_both <- inner_join(human_nuc_deg_rep1, human_nuc_deg_rep2, by='Gene')
nuc_deg_both$mean <- sqrt(nuc_deg_both$half_life_MAP.x * nuc_deg_both$half_life_MAP.y) #(nuc_deg_both$half_life_MAP.x + nuc_deg_both$half_life_MAP.y)/2
nuc_deg_both$mean_rate <- sqrt(nuc_deg_both$rate_MAP.x * nuc_deg_both$rate_MAP.y)#(nuc_deg_both$rate_MAP.x + nuc_deg_both$rate_MAP.y)/2

punds_temp <- inner_join(nuc_deg_both, gg_compare1_bayes_new, by='Gene')
punds_temp_rep1 <- punds_temp %>% filter(Rep == 'Rep 1')

for (rt in read_types){
# rt <- 'not_introns'#'introns'#'unspliced_junctions'#'exons'
  punds_temp2_rep1 <- punds_temp_rep1 %>% filter(RNA=='Nuclear RNA' & Sample!='wt' & nuc_deg.x=='No' & Read_type==rt)
  punds_temp3_rep1 <- punds_temp_rep1 %>% filter(RNA=='Nuclear RNA' & Sample!='wt' & nuc_deg.x=='Yes' & Read_type==rt)

  punds_temp2_rep1$line <- 2^(60/punds_temp2_rep1$mean)

  ggplot(punds_temp2_rep1, aes(x=log2(mean), y=log2(robert))) + #  log2(half_life_MAP.x)
    geom_line(aes(x=log2(mean), y=log2(line)), lwd=0.4) +
    geom_point(alpha = 0.05, col='gray50', cex=0.3) +
    facet_wrap(vars(Sample), nrow=2) +
    theme_bw(base_size = base_size) +
    scale_x_continuous(limits=c(0,11), breaks=c(0, 2, 4, 6, 8, 10), labels=c("1", "4", "16", "64", "256", "1,024")) +
    ylab("(1-MAPkd)/(1-MAPscr)") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background=element_rect(fill="white"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("Nuclear degradation half-life (minutes)") +
    scale_y_continuous(limits=c(-2,6), breaks=c(-2, 0, 2, 4, 6), labels=c("0.25", "1", "4", "16", "64")) +
    geom_point(data=punds_temp3_rep1, aes(x=log2(mean), y=log2(robert)), alpha = 0.1, col='red4', cex=0.3) +
    ggtitle(rt)
  ggsave(paste0(dir_out,"MAPs_delta_v_nuc_deg_half_life_rep1_log_nuc_all_genes_",rt,".pdf"), plot=last_plot(), width=2.5, height=2.5, units="in", dpi=300)

}#for rt


# cors
# blank table
cor_table <- data.frame(matrix(ncol=2))
colnames(cor_table) <- c("Sample","Spearman_cor")

samples <- c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1")
reps <- c("Rep 1", "Rep 2")
RNA <- c("Nuclear RNA")

i <- 1

for (RNA_type in RNA) {
  for (sample in samples) {
    for (rep in reps) {
      temp <- punds_temp %>% filter(Sample==sample & Rep==rep & RNA==RNA_type)
      c <- cor(log2(temp$mean), log2(temp$robert), method='spearman', use='complete.obs')
      c <- formatC( round( c, 2 ), format='f', digits=2 )
      c_text <- paste0("r == ", c)
      n_text <- paste0("n = ", nrow(temp))
      cor_table[i,1] <- paste0(RNA_type, " ", sample, " ", rep)
      cor_table[i,2] <- c_text
      i <- i+1
    }
  }
}



# load nuc deg rates (MAPs) from all_figures_load_rates.R

# look at punds only

nuc_deg_both <- inner_join(human_nuc_deg_rep1, human_nuc_deg_rep2, by='Gene')

# MAP rates here
nuc_deg_both$mean <- (nuc_deg_both$half_life_MAP.x + nuc_deg_both$half_life_MAP.y)/2
nuc_deg_both$mean_rate <- (nuc_deg_both$rate_MAP.x + nuc_deg_both$rate_MAP.y)/2

punds_temp <- inner_join(nuc_deg_both, gg_compare1_bayes_new, by='Gene')
punds_temp_rep1 <- punds_temp %>% filter(Rep == 'Rep 1')

ggplot(punds_temp_rep1, aes(x=log2(mean), y=kd_minus_scr)) +
  geom_point(alpha = 0.1) +
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  xlim(c(0,8)) +
  ylab("Change in fraction new RNA") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white")) +
  xlab("log2(Nuclear degradation half-life)") 
# ggsave(paste0(dir_out,"MAPs_delta_v_nuc_deg_half_life_rep1.pdf"), plot=last_plot(), width=5, height=3, units="in", dpi=300)

# log[ (1-MAP(KD)) / (1-MAP(scramble)) ]
ggplot(punds_temp_rep1, aes(x=log2(mean), y=log2(robert))) +
  geom_point(alpha = 0.1) +
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  scale_x_continuous(limits=c(0,11), breaks=c(0, 2, 4, 6, 8, 10), labels=c("1", "4", "16", "64", "256", "1,024")) +
  ylab("(1-MAPkd)/(1-MAPscr)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white")) +
  xlab("Nuclear degradation half-life") +
  scale_y_continuous(limits=c(-2,6), breaks=c(-2, 0, 2, 4, 6), labels=c("0.25", "1", "4", "16", "64"))
# ggsave(paste0(dir_out,"MAPs_delta_v_nuc_deg_half_life_rep1_log.pdf"), plot=last_plot(), width=5, height=3, units="in", dpi=300)

# same but just nuc
punds_temp2_rep1 <- punds_temp_rep1 %>% filter(RNA=='Nuclear RNA')

punds_temp2_rep1$line <- 2^(60/punds_temp2_rep1$mean)

plot <- punds_temp2_rep1 %>% filter(Sample!='wt')

ggplot(plot, aes(x=log2(mean), y=log2(robert))) +
  geom_line(aes(x=log2(mean), y=log2(line)), color='gray60', lwd=0.4) +
  geom_point(alpha = 0.1, col='red4', cex=0.3) +
  facet_wrap(vars(Sample), nrow=2) +
  theme_bw(base_size = base_size) +
  scale_x_continuous(limits=c(0,11), breaks=c(0, 2, 4, 6, 8, 10), labels=c("1", "4", "16", "64", "256", "1,024")) +
  ylab("(1-MAPkd)/(1-MAPscr)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Nuclear degradation half-life (minutes)") +
  scale_y_continuous(limits=c(-2,6), breaks=c(-2, 0, 2, 4, 6), labels=c("0.25", "1", "4", "16", "64"))
# ggsave(paste0(dir_out,"MAPs_delta_v_nuc_deg_half_life_rep1_log_nuc.pdf"), plot=last_plot(), width=2.5, height=2.5, units="in", dpi=300)


# log[ (1-MAP(KD)) / (1-MAP(scramble)) ]
plot <- punds_temp_rep1 %>% filter(Sample!='wt')

ggplot(plot, aes(x=mean_rate, y=log(robert))) +
  geom_point(alpha = 0.1) +
  facet_grid(cols=vars(Sample), rows=vars(RNA)) +
  theme_bw(base_size = base_size) +
  ylab("ln((1-MAPkd)/(1-MAPscr))") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background=element_rect(fill="white")) +
  xlab("k_nd") 
# ggsave(paste0(dir_out,"MAPs_delta_v_nuc_deg_rate.pdf"), plot=last_plot(), width=5, height=3, units="in", dpi=300)


# cors
# blank table
cor_table <- data.frame(matrix(ncol=2))
colnames(cor_table) <- c("Sample","Spearman_cor")

samples <- c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1")
reps <- c("Rep 1", "Rep 2")
RNA <- c("Total RNA", "Nuclear RNA")

i <- 1

for (RNA_type in RNA) {
  for (sample in samples) {
    for (rep in reps) {
      temp <- punds_temp %>% filter(Sample==sample & Rep==rep & RNA==RNA_type)
      c <- cor(log2(temp$mean), log2(temp$robert), method='spearman')
      c <- formatC( round( c, 2 ), format='f', digits=2 )
      c_text <- paste0("r == ", c)
      n_text <- paste0("n = ", nrow(temp))
      cor_table[i,1] <- paste0(RNA_type, " ", sample, " ", rep)
      cor_table[i,2] <- c_text
      i <- i+1
    }
  }
}


########### OLD TEMP TEST with UL background subtraction
#tot
for (rt in read_types){
  # compare kd MAP / scramble MAP - tot
  rep_prefixes <- c("rep1", "rep2")
  sample_prefix <- c("DIS3", "EXOSC10", "PABPN1", "ZFC3H1")
  for (sample in sample_prefix) {
    for (rep_prefix in rep_prefixes) {
      temp1 <- get(paste0(sample, "_tot_", rep_prefix,'_', rt))
      temp2 <- get(paste0("scramble_tot_", rep_prefix,'_', rt))
      temp3 <- inner_join(temp1, temp2, by='Gene')
      # temp3$kd_over_scr <- (temp3[,9] + pseudo_count)/(temp3[,18] + pseudo_count) #default
      # temp3$kd_minus_scr <- temp3[,9] - temp3[,18]
      # temp3$kd_percent_change <- 100*(temp3[,9] - temp3[,18])/(temp3[,18] + pseudo_count)
      # temp3$robert <- (1 - temp3[,9] + pseudo_count)/(1 - temp3[,18] + pseudo_count)
      temp3$kd_over_scr <- (pmax(temp3[,9] - temp3[,5], pseudo_count)) / (pmax(temp3[,18] - temp3[,14], pseudo_count)) ###TEMP: with total UL background substraction as correction
      temp3$kd_minus_scr <- pmax(temp3[,9] - temp3[,5], pseudo_count) - pmax(temp3[,18] - temp3[,14], pseudo_count)
      temp3$kd_percent_change <- 100*(pmax(temp3[,9] - temp3[,5], pseudo_count) - pmax(temp3[,18] - temp3[,14], pseudo_count)) / (pmax(temp3[,18] - temp3[,14], pseudo_count))
      temp3$robert <- (1 - pmax(temp3[,9] - temp3[,5], pseudo_count))/(1 - pmax(temp3[,18] - temp3[,14], pseudo_count))
      temp3 <- temp3 %>% mutate(Sample=sample, RNA="tot",  Rep=rep_prefix, Read_type=rt)
      temp4 <- temp3 %>% select(Gene, Symbol.x, kd_over_scr, kd_minus_scr, kd_percent_change, robert, Sample, RNA, Rep, Read_type)
      assign( paste0(sample, "_", rep_prefix, "_compare1"), temp3)
      gg_compare1 <- rbind(gg_compare1, temp4)
    }
  }
  

# nuc
for (sample in sample_prefix) {
  for (rep_prefix in rep_prefixes) {
    temp1 <- get(paste0(sample, "_nuc_", rep_prefix, '_', rt))
    temp2 <- get(paste0("scramble_nuc_", rep_prefix, '_', rt))
    # temp3 <- inner_join(temp1, temp2, by='Gene')  ##default
    # temp3$kd_over_scr <- (temp3[,5] + pseudo_count)/(temp3[,10] + pseudo_count)
    # temp3$kd_minus_scr <- temp3[,5]-temp3[,10]
    # temp3$kd_percent_change <- 100*(temp3[,5]-temp3[,10])/(temp3[,10] + pseudo_count)
    # temp3$robert <- (1 - temp3[,5] + pseudo_count)/(1 - temp3[,10] + pseudo_count)
    temp1_tot <- get(paste0(sample, "_tot_", rep_prefix, '_', rt))
    temp2_tot <- get(paste0("scramble_tot_", rep_prefix, '_', rt))
    temp3 <- inner_join(temp1, temp2, by='Gene')
    temp3 <- inner_join(temp3, temp1_tot, by='Gene')
    temp3 <- inner_join(temp3, temp2_tot, by='Gene')
    temp3$kd_over_scr <- (pmax(temp3[,5] - temp3[,15], pseudo_count)) / (pmax(temp3[,10] - temp3[,24], pseudo_count)) ###TEMP TEST: with total UL background substraction as correction
    temp3$kd_minus_scr <- pmax(temp3[,5] - temp3[,15], pseudo_count) - pmax(temp3[,10] - temp3[,24], pseudo_count)
    temp3$kd_percent_change <- 100*(pmax(temp3[,5] - temp3[,15], pseudo_count) - pmax(temp3[,10] - temp3[,24], pseudo_count)) / (pmax(temp3[,10] - temp3[,24], pseudo_count))
    temp3$robert <- (1 - pmax(temp3[,5] - temp3[,15], pseudo_count))/(1 - pmax(temp3[,10] - temp3[,24], pseudo_count))
    temp3 <- temp3 %>% mutate(Sample=sample, RNA="nuc", Rep=rep_prefix, Read_type=rt)
    temp4 <- temp3 %>% select(Gene, Symbol.x, kd_over_scr, kd_minus_scr, kd_percent_change, robert, Sample, RNA, Rep, Read_type)
    assign( paste0(sample, "_", rep_prefix, "_compare1"), temp3)
    gg_compare1 <- rbind(gg_compare1, temp4)
  }
}