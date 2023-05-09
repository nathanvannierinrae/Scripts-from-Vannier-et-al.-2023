############# Stack barplots abundances bacteria RNA #########
library(reshape2)
library(funrar)
library(ggplot2)
library(tidyverse)
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
abundances3 <- melt(relative_abundances, id.vars = "Strain")
abundances3<-inner_join(abundances3,taxonomy, by="Strain")

ggplot(abundances3, aes(x = variable, y = value, fill=Order  )) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=c("Rhizobiales"="#009933", "Xanthomonadales"="gold","Actinomycetales"="#995c00", "Caulobacterales"="darkslategray2", "Pseudomonadales"="#e68a00","Sphingomonadales"="#009999", "Burkholderiales"="#80ffaa", "Flavobacteriales"="#cccc00", "Bacillales"="#ffcc99"))


############# Stack barplots abundances Fungi RNA #########
library(reshape2)
library(funrar)
library(ggplot2)
library(tidyverse)
taxonomy=read.csv("taxo_fungi.csv", h=T, sep=";")
abundances=read.table("abundance_fungi.csv", h=T, row.names = 1)


colnames(abundances)<-c("Strain", "Matrix 1", "Matrix 2", "Roots 1","Matrix 3", "Roots 2", "Roots 3")
abundances=as.data.frame(abundances)
abundances=abundances[,c(1,4,6,7,2,3,5)]
relative_abundances=make_relative(as.matrix(t(abundances[,2:7])))
relative_abundances[is.na(relative_abundances)]<-0
relative_abundances=t(relative_abundances)
relative_abundances=as.data.frame(relative_abundances)
relative_abundances$Strain<-abundances$Strain
abundances3 <- melt(relative_abundances, id.vars = "Strain")
abundances3<-inner_join(abundances3,taxonomy, by="Strain")

ggplot(abundances3, aes(x = variable, y = value, fill =order  )) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=c("Mortierellales"="#2E2E2E","Glomerellales"="#B40486","Helotiales"="#FA5882","Hypocreales"="#B40404","Pleosporales"="#B18904","Sordiarales"="#FA8258","Xylariales"="#7c1fA4"))


############# Stack barplots abundances bacteria DNA #########
library(reshape2)
library(funrar)
library(ggplot2)
library(tidyverse)
taxonomy=read.csv("tax_table_DNA_bact.txt", h=T, sep="\t")
abundances=read.table("DNA_16S.txt", h=T, row.names = 1)

abundances[,7]<-rownames(abundances)
abundances=abundances[,c(7,1:6)]
colnames(abundances)<-c("Strain", "Roots 1","Roots 2","Roots 3", "Matrix 1", "Matrix 2", "Matrix 3")
abundances=as.data.frame(abundances)
relative_abundances=make_relative(as.matrix(t(abundances[,2:7])))
relative_abundances[is.na(relative_abundances)]<-0
relative_abundances=t(relative_abundances)
relative_abundances=as.data.frame(relative_abundances)
relative_abundances$Strain<-abundances$Strain
abundances3 <- melt(relative_abundances, id.vars = "Strain")
abundances3<-inner_join(abundances3,taxonomy, by="Strain")

ggplot(abundances3, aes(x = variable, y = value, fill=Order  )) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=c("Rhizobiales"="#009933", "Xanthomonadales"="gold","Actinomycetales"="#995c00", "Caulobacterales"="darkslategray2", "Pseudomonadales"="#e68a00","Sphingomonadales"="#009999", "Burkholderiales"="#80ffaa", "Flavobacteriales"="#cccc00", "Bacillales"="#ffcc99"))

############# Stack barplots abundances Fungi DNA #########
library(reshape2)
library(funrar)
library(ggplot2)
library(tidyverse)
taxonomy=read.csv("taxo_fungi.csv", h=T, sep=";")
abundances=read.table("DNA_ITS.txt", h=T, row.names = 1)
abundances[,10]<-rownames(abundances)

abundances=abundances[,c(10,1,2,4,7,5,6)]
colnames(abundances)<-c("old_names", "Roots 1", "Roots 2", "Roots 3","Matrix 1", "Matrix 2", "Matrix 3")
abundances=as.data.frame(abundances)
relative_abundances=make_relative(as.matrix(t(abundances[,2:7])))
relative_abundances[is.na(relative_abundances)]<-0
relative_abundances=t(relative_abundances)
relative_abundances=as.data.frame(relative_abundances)
relative_abundances$old_names<-abundances$old_names
abundances3 <- melt(relative_abundances, id.vars = "old_names")
abundances3$old_names<-as.factor(abundances3$old_names)
taxonomy$old_names<-as.factor(taxonomy$old_names)
abundances3<-inner_join(abundances3,taxonomy, by="old_names")

ggplot(abundances3, aes(x = variable, y = value, fill =order  )) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=c("Mortierellales"="#2E2E2E","Glomerellales"="#B40486","Helotiales"="#FA5882","Hypocreales"="#B40404","Pleosporales"="#B18904","Sordiarales"="#FA8258","Xylariales"="#7c1fA4"))





