library(reshape2)
library(funrar)
library(tidyverse)
library(dplyr)
taxonomy=read.csv("tax_table_DNA_bact.txt", h=T, sep="\t")
abundances=read.table("abundance_bacteria.csv", h=T, row.names = 1, sep=",")


colnames(abundances)<-c("Strain", "Matrix 1", "Matrix 2", "Roots 1","Matrix 3", "Roots 2", "Roots 3")
abundances=as.data.frame(abundances)
abundances=abundances[,c(1,4,6,7,2,3,5)]
relative_abundances=make_relative(as.matrix(t(abundances[,2:7])))
relative_abundances[is.na(relative_abundances)]<-0
relative_abundances=t(relative_abundances)
relative_abundances=as.data.frame(relative_abundances)
relative_abundances$Strain<-abundances$Strain

#####Select pathways to work on ######
list_OG=read.csv("top200_annotation.csv", h=T, sep="\t")
list_OG_selected=list_OG[which(list_OG$cog=="J"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="C"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="M"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="L"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="P"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="U"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="K"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="F"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="E"),1]
#list_OG_selected=list_OG[which(list_OG$cog=="G"),1]
setwd("LFC/")
list <- list.files("LFC/")
list=list[1:31]

result=matrix(nrow=length(list), ncol=4)
colnames(result)=c("mean_FC","sumFC","Strain", "root_colonization")
result=as.data.frame(result)
for (f in 1:length(list)) {
  tryCatch({
    temporary_data=read.csv(paste("LFC/",list[f], sep=""), h=T, sep=",")
    temporary_data=temporary_data[which(temporary_data$X%in%list_OG_selected),]
    temporary_data$Strain<-paste("",list[f], sep="")
    temporary_data=temporary_data %>% separate(Strain, c("Strain"), sep="O")
    temporary_data=temporary_data[which(is.na(temporary_data$log2FoldChange)=="FALSE"),]
    result$mean_FC[f]<-mean(temporary_data$log2FoldChange)
    result$sumFC[f]<-sum(temporary_data$log2FoldChange)
    
    result$Strain[f]<-temporary_data$Strain[1]
    relative_abundances_temp=relative_abundances[which(relative_abundances$Strain%in%temporary_data$Strain),]
    result$root_colonization[f]<-mean(t(relative_abundances_temp[1,-(4:7)]))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} 

result=result[which(is.na(result$mean_FC)=="FALSE"),]

strains=c("Root61","Root322","Root22","Root102","Root172","Root149","Root491","Root685","Root123D2","Root670","Root154","Root1277","Root569","Root68","Root1280","Root70","Root179","Root559","Root667","Root181")
result_20=result[which(result$Strain%in%strains),]


result_20_J=result_20
#result_20_C=result_20
#result_20_M=result_20
#result_20_L=result_20
#result_20_P=result_20
#result_20_U=result_20
#result_20_K=result_20
#result_20_F=result_20
#result_20_E=result_20
#result_20_G=result_20
colnames(result_20_J)[1]<-"mean_FC_J"
#colnames(result_20_C)[1]<-"mean_FC_C"
#colnames(result_20_M)[1]<-"mean_FC_M"
#colnames(result_20_L)[1]<-"mean_FC_L"
#colnames(result_20_P)[1]<-"mean_FC_P"
#colnames(result_20_U)[1]<-"mean_FC_U"
#colnames(result_20_K)[1]<-"mean_FC_K"
#colnames(result_20_F)[1]<-"mean_FC_F"
#colnames(result_20_E)[1]<-"mean_FC_E"
#colnames(result_20_G)[1]<-"mean_FC_G"

#############These data are used to make the regression plots as described below #######
library(car)
mod=lm(log(result_20_J$root_colonization)~result_20_J$mean_FC_J)
anova(mod)
summary(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))


p3 <- ggplot(result_20, aes(x=result_20_J$mean_FC_J, y=log(result_20_J$root_colonization)), fill=Strain) +
  geom_point() +
  geom_smooth(method=lm , color="Black", fill="#9aaba7", se=TRUE, fullrange = TRUE)+
  ggtitle("p-value: xxxx, Adjusted R-squared: xxxxx")+
  xlab("Mean log2FC Translation, ribosomal structure and biogenesis top 200")+
  ylab("Mean relative abundance in Roots (log-transformed)")+
  geom_text_repel(
    aes(label = Strain),
    color = "black",
    size = 9/.pt, # font size 9 pt
    point.padding = 0.8, 
    box.padding = 0,
    min.segment.length = 0,
    max.overlaps = 20,
    seed = 7654 # For reproducibility reasons
  )+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p3


