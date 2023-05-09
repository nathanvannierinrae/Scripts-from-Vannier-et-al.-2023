###############################Figure 5c boxplot#######################################################################
######################################################################################################
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(dunn.test)
library(PMCMR)
data=read.table("dataset_colonization.csv", sep=",",  h=T)
data=as.data.frame(data)
data=data[-which(data$genotype=="Negative"),]
data=data[-which(data$genotype=="EcpD"),]
data=data[-which(data$genotype=="FimC"),]
data=data[-which(data$genotype=="IMPH2"),]
data=data[-which(data$genotype=="PstB"),]
data=data[-which(data$genotype=="PstS"),]
data$genotype <- factor(data$genotype,
                        levels = c('WT',"PstABCS","ExbD1","TypA"),ordered = TRUE)
levels(data$genotype)[6]<-"PstABCS"
levels(data$genotype)[7]<-"ExbD1"
levels(data$genotype)[8]<-"TypA"

p<-data %>% ggplot( aes(x=data$genotype, y=data$cfu.mg, fill=data$genotype)) +geom_boxplot(,outlier.shape = NA, alpha=0.7) +
   scale_fill_manual("legend",values = c("WT" = "#440154FF", "PstABCS" = "#4ac16dFF","ExbD1" = "#a0da39FF","TypA" = "#fde725FF"), labels=c("WT","??PstABCS","??exbD1","??TypA")) +
   geom_jitter(color="black",size=0.95, alpha=0.8,shape = 19)+
   labs(y = "CFU per mg root freshweight", x="") 
p+ 
   geom_text(data=tibble(x=2, y=300),aes(x=x, y=y, label= "*"), inherit.aes = FALSE, size=8)+
   geom_text(data=tibble(x=3, y=300),aes(x=x, y=y, label= "*"), inherit.aes = FALSE, size=8)+
   geom_text(data=tibble(x=4, y=300),aes(x=x, y=y, label= "***"), inherit.aes = FALSE, size=8)+
   geom_hline(yintercept = 58.80252,  alpha=0.9, size=0.8, color="black",linetype='solid')+
   theme_bw() + 
   theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
         , axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
scale_x_discrete(labels=c("WT"="WT","PstABCS" ="??PstABCS", "ExbD1"="??ExbD1","TypA"="??TypA"))

###################################Figure 5e boxplot complementation ###################################################################
######################################################################################################

library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(dunn.test)
library(PMCMR)
cts=read.csv("Colonies_per_FW_complementation.txt", sep="\t", h=T)
cts=cts[-which(cts$R179=="cExbD2"),]
cts=cts[-which(cts$R179=="eVExbD"),]
cts=cts[-which(cts$R179=="eVTypA"),]

cts_ExbD=cts[-which(cts$R179=="cTypA"),]
cts_ExbD=cts_ExbD[-which(cts_ExbD$R179=="dTypA"),]
cts_ExbD=cts_ExbD[-which(cts_ExbD$R179=="eVTypA"),]

cts_typA=cts[-which(cts$R179=="cExbD"),]
cts_typA=cts_typA[-which(cts_typA$R179=="dExbD"),]
cts_typA=cts_typA[-which(cts_typA$R179=="eVExbD"),]

kruskal.test(cts$Colonies~cts$R179)
dunn.test(cts$Colonies,cts$R179,method="bh")
p.val.dunn<-c(0.1939,0.1479,0.2945,0.0035,0.0014,0.0697,0.0001)
p.adjust(p.val.dunn, method = "fdr", n = length(p.val.dunn))
?dunn.test


cts$R179 <- factor(cts$R179,
                   levels = c('WT','dExbD',"cExbD","dTypA","cTypA"),ordered = TRUE)

levels(cts$R179)[2]<-"ExbD1"
levels(cts$R179)[3]<-"ExbD1:ExbD1"
levels(cts$R179)[4]<-"TypA"
levels(cts$R179)[5]<-"TypA:TypA"


p<-cts %>%
   ggplot( aes(x=R179, y=Colonies, fill=R179)) +
   geom_boxplot(,outlier.shape = NA, alpha=0.7) +
   scale_fill_manual("legend",values = c("WT" = "#440154FF", "ExbD1" = "#a0da39FF","ExbD1:ExbD1" = "#a0da39FF", "TypA" = "#fde725FF","TypA:TypA" =  "#fde725FF"),
                     labels=c("WT","??ExbD1","??ExbD1:ExbD1","??TypA","??TypA:TypA"))+
   geom_jitter(color="black", size=0.95, alpha=0.8,shape = 19) +
   geom_hline(yintercept = 29.62253, alpha=0.9, size=0.8, color="black",linetype='solid')+  theme(
      legend.position="none",
   )+
   labs(y = "CFU per mg root freshweight", x="")
p+ 
   geom_text(data=tibble(x=2, y=195),aes(x=x, y=y, label= "*"), inherit.aes = FALSE, size=6)+
   geom_text(data=tibble(x=3, y=200),aes(x=x, y=y, label= "n.s"), inherit.aes = FALSE, size=5)+
   geom_text(data=tibble(x=4, y=195),aes(x=x, y=y, label= "**"), inherit.aes = FALSE, size=6)+
   geom_text(data=tibble(x=5, y=200),aes(x=x, y=y, label= "n.s"), inherit.aes = FALSE, size=5)+
   theme_bw() + 
   theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
         , axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
   scale_x_discrete(labels=c("WT"="WT","ExbD1"="??ExbD1","ExbD1:ExbD1"="??ExbD1:ExbD1","TypA"="??TypA","TypA:TypA"="??TypA:TypA"))

