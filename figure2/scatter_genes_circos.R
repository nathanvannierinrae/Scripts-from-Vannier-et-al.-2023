######################################################
###############PreparingKaryotypefile###############
######################################################
library(tidyverse)
library(dplyr)
library(stringr)

alldata=read.table("file_for_itol_bacteria.txt")
#######list 03-2020#######
strains=c("Root61","Root322","Root22","Root102","Root172","Root149","Root491","Root685","Root123D2","Root670","Root154","Root1277","Root569","Root68","Root1280","Root70","Root179","Root559","Root667","Root181")
######list 04-2020#######
strains=c("Root61","Root322","Root22","Root172","Root149","Root491","Root685","Root123D2","Root670","Root154","Root1277","Root569","Root68","Root1280","Root70","Root179","Root559","Root667","Root53", "Root101")




filtered_datas=alldata[which(alldata$Strain%in%strains),]
wd<-setwd("E:/data_14-03-20/filesForDESeq202004/filesForDESeq_bacteria_nodup/strain_genes")
ls()
list<-list.files("E:/data_14-03-20/filesForDESeq202004/filesForDESeq_bacteria_nodup/strain_genes")
list
dataFiles<-lapply(Sys.glob("*.csv"),read.csv)

combinedstrains=bind_rows(dataFiles, .id = "column_label")



reads_genes=read.table("E:/data_14-03-20/filesForDESeq202004/filesForDESeq_bacteria_nodup/strain.genes.qf",h=T)

head(combinedstrains)

colnames(combinedstrains)[2]<-"Name"

abundance=read.csv("file_for_itol_bacteria.txt", sep="")

root_abundance=abundance[,c(1,41)]
soil_abundance=abundance[,c(1,42)]
root_abundance=as.data.frame(root_abundance)
soil_abundance=as.data.frame(soil_abundance)
colnames(soil_abundance)[1]<-"Strain"
colnames(root_abundance)[1]<-"Strain"
abundances=inner_join(root_abundance,soil_abundance,by='Strain')
filtered_datas=left_join(filtered_datas,abundances,by="Strain")
filtered_datas=filtered_datas[order(-filtered_datas$mean_abundance_roots ),]

joined_genes=left_join(combinedstrains,reads_genes,by='Name')

joined_data_genes=left_join(filtered_datas,joined_genes,by='Strain')

joined_data_genes=joined_data_genes[which(is.na(joined_data_genes$log2FoldChange)==FALSE),]

#########PourgarderlesNA########################################################
########joined_data_genes[,48:49][is.na(joined_data_genes[,48:49])]<-1###########
##############joined_data_genes[,46:47][is.na(joined_data_genes[,46:47])]<-0###########
joined_data_genes
names_1=filter(joined_data_genes,joined_data_genes$log2FoldChange>0.25)
names_2=filter(joined_data_genes,joined_data_genes$log2FoldChange<(-0.25))


names=rbind(names_1,names_2)
joined_data_genes=joined_data_genes[which(joined_data_genes$Name%in%names$Name),]


nb_genes_strain=joined_data_genes%>%
group_by(Strain)%>%
summarise(no_rows=length(Strain))

joined_data_genes=left_join(joined_data_genes,nb_genes_strain,by='Strain')

filtered_datas=left_join(filtered_datas,nb_genes_strain,by='Strain')

karyotype=as.data.frame(filtered_datas$Strain)
karyotype[,1]<-"chr"
karyotype[,2]<-"-"
karyotype[,3]<-filtered_datas$Strain
karyotype[,4]<-filtered_datas$Strain
karyotype[,5]<-"0"
karyotype[,6]<-filtered_datas$no_rows
karyotype[,7]<-filtered_datas$Strain



write.table(karyotype,"karyotype_genes_21strains_09-04.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

######generating the file for scatter circos##########

#########Creating the plotting areas for each gene####################
budubuh=joined_data_genes

list_terms<-as.factor(filtered_datas$Strain)
output_col=matrix(ncol=3,nrow=nrow(budubuh)+1)
output_rest=matrix(ncol=57,nrow=1)
colnames(output_rest)=colnames(budubuh)
a=0

for(i in 1:length(list_terms)){
dat=budubuh[which(budubuh$Strain%in%list_terms[i]),]
lengthdat=nrow(dat)
output=matrix(ncol=2,nrow=lengthdat)
b=1

for(j in 1:lengthdat){
output[j,1]<-b
output[j,2]<-b
b=b+1
}

a=a+lengthdat

output_col[(a+1-lengthdat):a+1,1:2]<-output

output_rest=rbind(output_rest,dat)

}
output_data=cbind(output_rest,output_col)
output_data=output_data[-1,1:59]
write.table(output_data,"table_for_circos_genes.txt")

#########Creating the scatter file for circos####################

table_circos=read.table("table_for_circos_genes.txt",h=T)
table_circos[,60]<-""
head(table_circos)
list_names_pvalue=table_circos[which(table_circos$padj<0.05),c(44,46,49)]
up_reg=list_names_pvalue[which(list_names_pvalue$log2FoldChange>0),]
dow_reg=list_names_pvalue[which(list_names_pvalue$log2FoldChange<0),]
table_circos[which(table_circos$Name%in%up_reg$Name),60]<-"color=238,50,50"
table_circos[which(table_circos$Name%in%dow_reg$Name),60]<-"color=88,77,232"
regulated=rbind(up_reg,dow_reg)
table_circos[which(!(table_circos$Name%in%regulated$Name)),60]<-"color=161,162,175"

table_circos[which(table_circos$Name%in%up_reg$Name),60]<-"color=red"
table_circos[which(table_circos$Name%in%dow_reg$Name),60]<-"color=blue"
regulated=rbind(up_reg,dow_reg)
table_circos[which(!(table_circos$Name%in%regulated$Name)),60]<-"color=grey"





scatter=table_circos[,c(1,58,59,46,60)]
scatter$log2FoldChange<-round(scatter$log2FoldChange,digits=5)
write.table(scatter,"scatter_genes_21strains_08-04.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

summary(scatter$log2FoldChange)

#########Creating the highlight file for the taxonomy in circos####################
karyotype=read.table("karyotype_genes_20strains_09-04.txt")

highlight=as.data.frame(karyotype$V3)
highlight[,2]<-0
highlight[,3]<-karyotype$V6
highlight[,4]<-c("fill_color=255,255,0","fill_color=180,95,4","fill_color=180,95,4","fill_color=0,153,51","fill_color=0,153,51","fill_color=230,138,0","fill_color=255,255,0","fill_color=128,255,170","fill_color=0,153,51","fill_color=0,153,153","fill_color=230,138,0","fill_color=0,153,51","fill_color=230,138,0","fill_color=180,95,4","fill_color=255,255,0","fill_color=230,138,0","fill_color=0,153,51","fill_color=179,255,255","fill_color=0,153,51","fill_color=180,95,4")

"fill_color=255,255,0"
"stroke_thickness=2,stroke_color=white"

c("255,255,0","180,95,4","180,95,4","0,153,51","0,153,51","230,138,0","255,255,0","128,255,170","0,153,51","0,153,153","230,138,0","0,153,51","230,138,0","180,95,4","255,255,0","230,138,0","0,153,51","179,255,255","0,153,51","180,95,4")
write.table(highlight,"highlight_genes_20strains_08-04.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")



#########Creating the barplots file for the abundances in circos####################


table_circos_abundances=left_join(table_circos,abundances,by='Strain')

histogram_soil=table_circos_abundances[,c(1,85,86,90)]
write.table(histogram_soil,"histograme_soil_20strains_08-04.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


histogram_roots=table_circos_abundances[,c(1,58,59,61)]
write.table(histogram_roots,"histograme_roots_21strains_08.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

summary(scatter$log2FoldChange)


