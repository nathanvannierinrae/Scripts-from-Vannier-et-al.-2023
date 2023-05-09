library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(vegan)
library(tidyverse)
####Matrix vs Roots####
#####recovering normalized counts#######
cts=read.csv("strain.OG.qf", sep="\t", h=T)
for (f in 1:length(levels(cts$Strain))) {
  tryCatch({
    counts=cts[which(cts$Strain==levels(cts$Strain)[f]),]
    cts2=sapply(counts[,c(3:8)], as.integer)
    row.names(cts2)=counts[,2]
    cond=factor(c("M","M","R", "M","R","R"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    h=counts(dds, normalized=TRUE)
    h=as.data.frame(h)
    h$Strain<-paste("",levels(cts$Strain)[f], sep="")
    write.csv(h, paste("normcounts/", levels(cts$Strain)[f], "_normcounts_MvsR.csv", sep = ""))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} 


wd <-setwd("normcounts")
ls()
list <- list.files("normcounts")
list
dataFiles <- lapply(Sys.glob("*.csv"), read.csv)

combinedstrains=bind_rows(dataFiles, .id = "column_label")
strains=c("Root61","Root322","Root22","Root102","Root172","Root149","Root491","Root685","Root123D2","Root670","Root154","Root1277","Root569","Root68","Root1280","Root70","Root179","Root559","Root667","Root181")

combinedstrains=combinedstrains[which(combinedstrains$Strain%in%strains),]
tab=combinedstrains[,c(2,9,3:8)]
colnames(tab)[1]="OG" 

###################
###reshape table###
###################


tabSB1=tab[,c(1,2,3)]
tabSB2=tab[,c(1,2,4)]
tabSB3=tab[,c(1,2,6)]
tabSM1=tab[,c(1,2,5)]
tabSM2=tab[,c(1,2,7)]
tabSM3=tab[,c(1,2,8)]

resSB1=reshape(tabSB1, idvar=c("OG"), timevar="Strain", direction="wide")
resSB2=reshape(tabSB2, idvar=c("OG"), timevar="Strain", direction="wide")
resSB3=reshape(tabSB3, idvar=c("OG"), timevar="Strain", direction="wide")
resSM1=reshape(tabSM1, idvar=c("OG"), timevar="Strain", direction="wide")
resSM2=reshape(tabSM2, idvar=c("OG"), timevar="Strain", direction="wide")
resSM3=reshape(tabSM3, idvar=c("OG"), timevar="Strain", direction="wide")

resSB1=resSB1[,-1]
resSB2=resSB2[,-1]
resSB3=resSB3[,-1]
resSM1=resSM1[,-1]
resSM2=resSM2[,-1]
resSM3=resSM3[,-1]



resSB1[is.na(resSB1)] <- 0
resSB2[is.na(resSB2)] <- 0
resSB3[is.na(resSB3)] <- 0
resSM1[is.na(resSM1)] <- 0
resSM2[is.na(resSM2)] <- 0
resSM3[is.na(resSM3)] <- 0





#########################################################################
#############Calculate mean expression in roots/matrix samples###########
#########################################################################


SB_SM=cbind(resSB1,resSB2,resSB3,resSM1,resSM2,resSM3)
SB_SM=as.data.frame(SB_SM)

####

tax_list<-levels(as.factor(combinedstrains$Strain))
tax_list=as.data.frame(tax_list)
colnames(tax_list)[1]<-"Strain"
colnames(SB_SM)=tax_list$Strain
taxonomy=read.csv("tax_table_DNA_bact.txt", sep="", h=T)
tax_list=inner_join(tax_list,taxonomy, by="Strain")

tax_list$sample<-"Matrix"


SB_tax_list=tax_list

tax_list$sample<-"Roots"
SM_tax_list=tax_list

sample_data=rbind(SM_tax_list,SM_tax_list,SM_tax_list,SB_tax_list,SB_tax_list,SB_tax_list)



###PCOA roots+soil###
library(phyloseq)

sample_data$couleur<-"brown"
sample_data[which(sample_data$sample=="Roots"),10]<-"green"
samples = sample_data(sample_data)
colnames(SB_SM)=rownames(samples)
OTU = otu_table(SB_SM, taxa_are_rows = TRUE)
GP1 <- phyloseq(OTU,samples)
GP1_norm  = transform_sample_counts(GP1, function(x) x / sum(x) )
GP1_filter=filter_taxa(GP1_norm, function(x) mean(x) > 0.0001, TRUE)


GP.ord <- ordinate(GP1_filter, "MDS", "bray")
p1 = plot_ordination(GP1_filter,  GP.ord, axes=c(3,4),type="samples", shape="sample", color="sample")
p1  + geom_point(size=3)  + scale_colour_manual(values=c("Roots"="forestgreen", "Matrix"="chocolate2"))+  ggtitle("Expression profile of OG in Matrix and Roots samples") 

p1 = plot_ordination(GP1_filter,  GP.ord, axes=c(3,4),type="samples", shape="sample", color="Order")
p1  + geom_point(size=3)  + scale_colour_manual(values=c("Rhizobiales"="#009933", "Xanthomonadales"="gold","Actinomycetales"="#995c00", "Caulobacterales"="darkslategray2", "Pseudomonadales"="#e68a00","Sphingomonadales"="#009999", "Burkholderiales"="#80ffaa"))+  ggtitle("Expression profile of OG in Matrix and Roots samples") 
