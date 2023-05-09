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
data$genotype <- factor(data$genotype,
                        levels = c('WT','FimC',"IMPH2","PstB","PstS","PstABCS","ExbD1","TypA"),ordered = TRUE)

levels(data$genotype)[2]<-"??FimC"
levels(data$genotype)[3]<-"??IMPH2"
levels(data$genotype)[4]<-"??PstB"
levels(data$genotype)[5]<-"??PstS"
levels(data$genotype)[6]<-"??PstABCS"
levels(data$genotype)[7]<-"??ExbD1"
levels(data$genotype)[8]<-"??TypA"

kruskal.test(data$cfu.mg~data$genotype)
dunn.test(data$cfu.mg,data$genotype)
p.val.dunn<-c(0.2242,0.0086,0.1009,0.0844,0.0142,0.1714,0.2545,0.0000)
p.adjust(p.val.dunn, method = "fdr", n = length(p.val.dunn))
dunn.test.control (data$cfu.mg, data$genotype, p.adjust.method = p.adjust.methods, .)
p<-data %>% ggplot( aes(x=data$genotype, y=data$cfu.mg, fill=data$genotype)) +geom_boxplot(,outlier.shape = NA, alpha=0.7) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, labels=c("WT","??FimC","??IMPH2","??PstB","??PstS","??PstABCS","??exbD1","??TypA")) +
  geom_jitter(color="black",size=0.95, alpha=0.8,shape = 19) 

p+ 
  geom_text(data=tibble(x=2, y=615),aes(x=x, y=y, label= "n.s"), inherit.aes = FALSE, size=5)+
  geom_text(data=tibble(x=3, y=615),aes(x=x, y=y, label= "n.s"), inherit.aes = FALSE, size=5)+
  geom_text(data=tibble(x=4, y=615),aes(x=x, y=y, label= "n.s"), inherit.aes = FALSE, size=5)+
  geom_text(data=tibble(x=5, y=615),aes(x=x, y=y, label= "n.s"), inherit.aes = FALSE, size=5)+
  geom_text(data=tibble(x=6, y=600),aes(x=x, y=y, label= "*"), inherit.aes = FALSE, size=6)+
  geom_text(data=tibble(x=7, y=600),aes(x=x, y=y, label= "*"), inherit.aes = FALSE, size=6)+
  geom_text(data=tibble(x=8, y=600),aes(x=x, y=y, label= "***"), inherit.aes = FALSE, size=6)+
  geom_hline(yintercept = 58.80252,  alpha=0.9, size=0.8, color="black",linetype='solid')+
  theme_bw() + 
  theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
        , axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))



