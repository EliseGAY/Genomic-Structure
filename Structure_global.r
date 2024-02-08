#!/usr/bin/env Rscript

#========================================================================================#
#========================================================================================#
# Pop structure : 

  # - ACP Global on all chromosome 
  # - FST
  # - sNMF
  # - SFS
  # - DIV

# Aim : Run general analysis of population structure on gVCF

# Authors : Stefano Mona , Elise Gay, Romuald Laso-Jadart 
#           2023
#	Please inform the authors before sharing
#========================================================================================#
#========================================================================================#

#============================#
#============================#
# ------ Load libraries ----
#============================#
#============================#
library(LEA)
library(vcfR)
library(rlist)
library(PopGenome)
library(qvalue)
library(pegas)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(withr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(pcadapt)

#--------------#
# Metadata pop
#--------------#

# the pop table has to be ordered in the same way as all the VCF header
#'''
#samples	pop	
#sample_1	pop1
#sample_2	pop1
#sample_3	pop1
#sample_4	pop2
#sample_5	pop2
#'''

table_pop=read.table("metadata/Samples_table.txt", header = TRUE, row.names = 1)
table_pop

samples=row.names(table_pop)
samples

pop=table_pop$pop
pop

pop1=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP1"),]))
pop1
pop2=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP2"),]))
pop2
pop3=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP3"),]))
pop3
pop4=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP4"),]))
pop4

lista_pop<-list(pop1,pop2,pop3,pop4)
names(lista_pop)=c("pop1","pop2","pop3","pop4")

# for div computing 
lista_pop_all=list(pop1,pop2,pop3,pop4,c(pop1,pop2,pop3,pop4))
names(lista_pop_all) = c("pop1","pop2","pop3","pop4","pop1+pop2+pop3+pop4")

#-----------------------------#
# Load Stefano Mona functions
#-----------------------------#
source("R_functions/libreria_filtri_VCF_e_SFS_unfolded.r")

#===============================================#
#===============================================#
# ------ Pairwise FST on whole chr ----
#===============================================#
#===============================================#

#----------#
# Infos  : 
#----------#

# INPUT 
#-------#
# VCF with chosen filters : 20% Na and MAF = 0.05
# Position vector : vector of total position called in the chr

# FUNCTION ARGS
#---------------------#
# res_FST<-calcola_fst_pairwise_bootstrap(data_FST, lista_pop = lista_pop, MAF, bootstrap)

# data_FST : input VCF
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part
# MAF : Minor allele frequency chosen (usually 0.05)
# bootstrap : nb of repetition. Bootstrap are made by computing FST on random pairs of population

# OUTPUT
#---------------------#
# table of pairwise FST values in each windows 
# "mid_wind","FST","low_bound", "upper_bound","Nb_SNP_FST"

# Load data / samples set up
#-----------------------------#
# Load input VCF depending on analysis
data_FST=read.vcfR("data/global_FST/super_1/WGS_21_SUPER_1_TAG_Flowqual_Noindels_Norepeat_SNP_DP10_50_10_200_Na20_MAF005_order.vcf")
colnames(data_FST@gt)

# Compute FST pairwise
#----------------------#
# Usually : MAF = 0.01 or 0.05 / Boot = 1000
res_FST=calcola_fst_pairwise_bootstrap(data_FST, 
                                       lista_pop = lista_pop, 
                                       0.05, 
                                       10)

# res_FST[1] : pop1-pop2 pop1-pop3 pop1-pop4 pop2-pop3 pop2-pop4 pop3-pop4
# [1] 0.111517482 0.096383436 0.154419450 0.006691872 0.067096883 0.064647038

# res_FST[2]
# [,1]         [,2]          [,3]          [,4]          [,5]          [,6]
# pairwise_fst_temp  0.021129757 -0.004930209 -0.0103571582  1.395912e-03  1.182520e-02  0.0113430827
# pairwise_fst_temp  0.016506281  0.008023368  0.0001254892 -2.019709e-03  2.076501e-02  0.0154640942
# pairwise_fst_temp -0.005050721 -0.016065455  0.0001254892 -1.100178e-03  1.739254e-02  0.0174395518
# pairwise_fst_temp  0.025454615 -0.007660121  0.0001254892 -2.820894e-03  4.464543e-02  0.0126687982
# pairwise_fst_temp  0.017705668  0.015680395  0.0001344521  4.903638e-05 -6.079505e-03  0.0076506819

# compute significance fst
#---------------------------#
colnames(res_FST[[2]]) = c("POP1-POP2","POP1-POP3","POP1-POP4","POP2-POP3","POP2-POP4","POP3-POP4")

# plot FST distribution among bootstrap and compare with real fst value
pop_list=seq(1,6) # for each pairewise pop : 1 to 6

for (i in pop_list){
  p=ggplot()+
    geom_density(data=as.data.frame(res_FST[[2]][,i]), 
                 aes(x= `res_FST[[2]][, i]`,
                     stat = "density")) +
    geom_vline(xintercept = res_FST[[1]][i])
  
  plot_name=paste("FST_plot_POP_comp",i, sep = "_")
  assign(plot_name, p, envir=parent.frame())
  
  # get p-value = mean fst bootstraps value which are >= 'real' fst
  pval_com_name=paste("pval_comp", i, sep = "_" )
  pval_com= sum(res_FST[[2]][,i]>=res_FST[[1]][i])/10
  assign(pval_com_name, pval_com, envir=parent.frame())
  
  # IF needed : get the absolute number of fst boot, which are >= to 'real fst'
  # table(res_FST[[2]][,i]>=res_FST[[1]][i])
}

FST_plot_POP_comp_1
FST_plot_POP_comp_2
FST_plot_POP_comp_3
FST_plot_POP_comp_4
FST_plot_POP_comp_5
FST_plot_POP_comp_6

# ------ FST global ------ 
#=========================================#
#=========================================#
# FST global in sliding windows on all pop
#=========================================#
#=========================================#

#----------#
# Infos  : 
#----------#

# INPUT 
#-------#
# VCF with chosen filters : here 20% Na + MAF 0.05
# Because VCF will be read by popgenome : VCF HAS TO BE GUNZIPED AND PUT IN ONE FILE by FOLDER

# FUNCTION ARGS
#---------------------#

# boot_popgenome<-function(dati, lista_pop, bootstrap)
# dati = vcf
# bootstrap (int)
# lista_pop = list of population 

# OUTPUT
#---------------------#
# - global FST with bootstrap distribution 
# - pairwise FST with bootstrap distrib

# STEP 1 
#---------#
# Gunzip each vcf and put each VCF in individuals folder name 'Chr_vcf_popgenome'

#------------------------#
# STEP 2 : Run the loop
#------------------------#
bootstrap=10
# read data
vcf_data=readData("data/global_FST/super_1/",format = "VCF")
res=boot_popgenome(vcf_data, lista_pop, bootstrap)

#-------------#
# FST summary
#-------------#

# results with "calcola_fst_pairwise_bootstrap" : reynolds FST
# [1] 0.106380006 0.092142425 0.154419450 0.006612063 0.070189834 0.068777718

# FST global with popgenome
# Fst                  2.5%                   50%                 97.5%            nloci/snps 
# "0.13558030765765" "-0.0240334320084205" "0.00237811967706347"  "0.0263811548990177"           "273518738" 

# FST pairwise with popgenome
# result with popgenome : Fixation Index based on minor.allele frequencies (Hudson)
# pop_1/pop_2 0.14861414 -0.04030557  0.009117634 0.07403241
# pop_1/pop_3 0.14053908 -0.02745922  0.007356764 0.03605225
# pop_1/pop_4 0.22905353 -0.04480268 -0.019159681 0.22905353
# pop_2/pop_3 0.01895009 -0.02282777 -0.001567434 0.02103520
# pop_2/pop_4 0.12904505 -0.04178652  0.002703320 0.06282386
# pop_3/pop_4 0.12324716 -0.02498501  0.004716508 0.04277268

#======================# 
#======================# 
# ------ ACP All chr----
#======================#
#======================#

# Compute PCA on the entire VCF given

# INPUT 
#-------#
# VCF with chosen filters : 20% Na and MAF filtered (MAF can also be put directly in the PCadapt function)
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part

# FUNCTION ARGS
#------------------#
# respca_10 = pcadapt(input = read.pcadapt(path_vcf, type = "vcf"),  min.maf=0.01, K = 10)
  # VCF : abslotue path to your VCF file (already filtered)
  # min.maf = minimum maf threshold
  # K = number of components to be computed

# OUTPUT
#---------#
# PCA plot according to  and population assignation (see metadata)

#----------#
# run ACP
#----------#
# Run PCA
respca_10 = pcadapt(input = read.pcadapt("data/PCA_FST/WGS_21_SUPER_1_TAG_Flowqual_Noindels_Norepeat_SNP_DP10_50_10_200_Na20_MAF005_order.vcf", 
                                         type = "vcf"),
                    min.maf=0.05, K = 10)

# add samples name in table result (same order as VCF)
scores = data.frame(respca_10$scores)
rownames(scores) = samples

# Compute percent of variance explained by all axis (here 10 axis)
PEV = paste(round((respca_10$singular.values^2)*100,digits=1),"%",sep="")
PEV
# "11.5%" "10.6%" "6.6%"  "6.2%"  "6%"    "5.8%"  "5.6%"  "5.5%"  "5.4%"  "5.4%" 

# Plot PCA
ggplot() +
    geom_point(data=scores,
               aes(x=scores$X1, y=scores$X2,
                   colour=pop),
               size=2)+

    labs(x=PEV[1], y = PEV[2]) +
    scale_colour_manual(name = "Sampling sites",
                         values=rainbow(length(unique(pop)))) +

    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5))

#=====================#
#=====================#
# ---- sNMF ----
#=====================#
#=====================#
# INPUT 
#-------#
  # VCFs filtered for Na and DP on R (for NA and DP)
  #	==> MAF at 0.005 added with bcftools
  #	==> Samples were ordered by population (formated file for snmf anlaysis)
  #	==> On this file : plink was run to made the ".ped" file needed to run snmf.

# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part

# METHODS
#------------------#
# snmf : "estimates admixture coefficients using sparse Non-Negative Matrix Factorization algorithms"

# OUTPUT
#---------#
# table_all : contains all cross entropy for 10 rep in 100 snmf runs in each K.
# A folder is created  containing the snmf project
# plot of admixture for each individuals
 
# Run sNMF
#----------#
# go in the directory. The absolute path has to be short, ortherwise snmf won't work
setwd("data/snmf/")
# create project with nb K + nb repetition chosen
project = snmf("WGS_21_SUPER_1_TAG_Flowqual_Noindels_Norepeat_SNP_DP10_50_10_200_Na20_MAF005_order_Plink.ped",
               K=1:8,
               entropy=T,
               repetitions = 20,
               project = "new")

# To re-load the project already created, use:
project = load.snmfProject("VCF.Plink.snmfProject")

# plot cross-entropy criterion of all runs of the project
plot(project, cex = 1.2, col = "lightblue", 
     pch = 19, 
     xaxp=c(0,35,35))

# Run sNMF with repetition
#--------------------------#
# One hundred Runs of snmf on K = 8 and with 10 repetition

# initiate table for K = 8
table_all=data.frame("K=2"=numeric(), 
                      "K=3"=numeric(), 
                      "K=4"=numeric(),
                     "K=5"=numeric(), 
                     "K=6"=numeric(), 
                     "K=7"=numeric(),
                     "K=8"=numeric())
table_all

# To choose the best K, with repetition :
# Ten Runs of snmf on K = 8 and with 10 repetition each : 100 run of snmf
for (i in seq(1:10)){
  project = snmf("WGS_21_SUPER_1_TAG_Flowqual_Noindels_Norepeat_SNP_DP10_50_10_200_Na20_MAF005_order_Plink.ped",
                 K=1:8,
                 entropy=T,
                 repetitions = 10,
                 project = "new")
  table_i=as.data.frame(cbind(cross.entropy(project, K = 2),
                cross.entropy(project, K = 3),
                cross.entropy(project, K = 4),
                cross.entropy(project, K = 5),
                cross.entropy(project, K = 6),
                cross.entropy(project, K = 7),
                cross.entropy(project, K = 8)))
  table_all=rbind(table_all, table_i)
}

# format table to plot
summary(table_all)
table_all_melt=melt(table_all)

# plot boxplot of cross entropy
ggplot()+
  geom_boxplot(aes(x=table_all_melt$variable, 
              y=table_all_melt$value))

# plot admixture
#----------------#
for (i in c(1,2,3)) {

  ce = cross.entropy(project, K = i)
  best = which.min(ce)
  qmatrix = Q(project, K = i,run = best)
  rownames(qmatrix)= samples # carefull of the samples order in list
  meltedqmatrix = melt(qmatrix)
  
  #One color for one K
  my.colors=c("orange","violet","lightgreen","grey40","dodgerblue","goldenrod","firebrick2","forestgreen", rainbow(7))

  
  p=plot(ggplot() +
         geom_bar(data=meltedqmatrix,aes(x=Var1,y=value,fill=Var2),stat="identity",
                  show.legend=F) +
         theme(text = element_blank(), 
               #axis.text = element_text(size = 7, angle = 90),
               axis.text.x=element_blank(),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                               colour = "grey"),
               panel.background = element_rect(fill = "white",
                                               colour = "black",
                                               size = 0.5,
                                               linetype = 'solid'),
               axis.ticks = element_blank()) +
         xlab("Individuals") +
         ylab("Ancestry proportion") +
         ggtitle(paste("K = ",i,", Minimal cross-entropy = ",ce,sep="")) +
         scale_fill_manual(values=my.colors))+
    geom_vline(xintercept = c(47, 71, 105), size=1,linetype="dashed")
	
	name=paste(i, "plot.pdf", sep="")
	ggsave(p, filename = name, width = 10, height = 10, device = "pdf")

}

#==============#
#==============#
# ---- SFS ----
#==============#
#==============#

# INPUT 
#-------#
# VCF with required filters (generally : Na = 0%, Heterozygous rate max = 80%)

# METHODS
#------------------#
# PEGAS function : site.spectrum. Create SFS on DNAbin sequences
# calcola_normalized_foldedSFS : used to fold the SFS

# OUTPUT
#---------#
# sfs table

# change directory if needed:
setwd("~/../OneDrive - MNHN/museum/Shared_tools/Pop_structure/examples_structure_analyses/")
# read VCF :
VCF=read.vcfR("data/SFS/WGS_21_SUPER_1_TAG_Flowqual_Noindels_Norepeat_SNP_DP10_50_10_200_Na0_het80_order.vcf")
# Create DNA sequences
DNAbin <-vcfR2DNAbin(VCF,  extract.indels = TRUE,
                     consensus = FALSE,
                     extract.haps = F,
                     unphased_as_NA = F,
                     asterisk_as_del = FALSE,
                     ref.seq = NULL,
                     start.pos = NULL,
                     verbose = TRUE)


# Compute folded site frequency spectrum with the function site.spectrum (library pegas required)
sfs_folded_GWS<-site.spectrum(DNAbin, folded=T)
sfs_folded_GWS

# Normalized by the spectrum by the number of sites 
# Stefano's function: calcola_normalized_foldedSFS
# ARGUMENT : Folded SFS
sfs_folded_GWS_norm<-calcola_normalized_foldedSFS(sfs_folded_GWS)
sfs_folded_GWS_norm_table=as.data.frame(sfs_folded_GWS_norm)

# Plot 
ggplot() +
  geom_point(aes(y=sfs_folded_GWS_norm_table$eta_2, x=seq(1,21)), color="blue")+
  labs(x="", y = "")+
  scale_x_continuous(breaks = seq(0, 21, by = 1), limits = c(0,21))+
  scale_y_continuous(breaks = seq(0,1, by = 0.1), limits = c(0,1))


#==============================#
#==============================#
# ---- diversity indexes ----
#==============================#
#==============================#

# INPUT 
#-------#
# VCF with required filters (generally : Na = 0%, Heterozygous rate max = 80%)

# METHODS
#------------------#
# PEGAS function : site.spectrum. Create SFS on DNAbin sequences
# calcola_normalized_foldedSFS : used to fold the SFS

# OUTPUT
#---------#
# 4 diversity indexes : 
#   - Nb de sites ségregeants
#   - Mean pairwise difference
#   - S
#   - tajima D

div_list=calcola_TD_folded(sfs_folded_GWS_norm)

div_list

# > div_list
# [[1]]
# [1] 2.147975
# 
# [[2]]
# [1] 3.66808
# 
# [[3]]
# [1] 2.01861
# 
# [[4]]
# [1] 10.09723