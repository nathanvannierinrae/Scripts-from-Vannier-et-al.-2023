##############
## Bacteria ##
##############

### Deseq strain genes  ####
cts=read.csv("strain.genes.qf", sep="\t", h=T,row.names = 1)
library("DESeq2")
library("apeglm")
cts$Strain=as.factor(cts$Strain)

for (f in 1:length(levels(cts$Strain))) {
  tryCatch({
    counts=cts[which(cts$Strain==levels(cts$Strain)[f]),]
    counts=counts[,3:8]
    cts2=sapply(counts, as.integer)
    row.names(cts2)=row.names(counts)
    cond=factor(c("A","A","B","A","B","B"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)
    results.All <- lfcShrink(dds, coef=2, type="apeglm")
    results.de <- subset(results.All, padj < 0.05)
    results.oe <- subset(results.All, padj < 0.05 & log2FoldChange > 0)
    results.ue <- subset(results.All, padj < 0.05 & log2FoldChange < 0)
    norm_counts<-counts(dds)
    write.csv(norm_counts, paste("F:/data_14-03-20/nouvelle_analyse_07-10-2021/filesForDESeq_bacteria_koGrouped_nodup/norm_counts/", levels(cts$Strain)[f], "_deseq.csv", sep = ""))
    write.csv(results.All, paste("F:/data_14-03-20/nouvelle_analyse_07-10-2021/filesForDESeq_bacteria_koGrouped_nodup/all_genes/", levels(cts$Strain)[f], "_deseq.csv", sep = ""))
    write.csv(results.de, paste("F:/data_14-03-20/nouvelle_analyse_07-10-2021/filesForDESeq_bacteria_koGrouped_nodup/de_genes/",levels(cts$Strain)[f], "_deseq.csv", sep = ""))
    write.csv(results.oe, paste("F:/data_14-03-20/nouvelle_analyse_07-10-2021/filesForDESeq_bacteria_koGrouped_nodup/oe_genes/",levels(cts$Strain)[f], "_deseq.csv", sep = ""))
    write.csv(results.ue, paste("F:/data_14-03-20/nouvelle_analyse_07-10-2021/filesForDESeq_bacteria_koGrouped_nodup/ue_genes/",levels(cts$Strain)[f], "_deseq.csv", sep = ""))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} 



### Deseq strain OG  ####
cts=read.csv("strain.OG.qf", sep="\t", h=T)
cts$Strain=as.factor(cts$Strain)
for (f in 1:length(levels(cts$Strain))) {
  tryCatch({
    counts=cts[which(cts$Strain==levels(cts$Strain)[f]),]
    rownames(counts)=counts$OG
    counts=counts[,3:8]
    cts2=sapply(counts, as.integer)
    row.names(cts2)=row.names(counts)
    cond=factor(c("A","A","B","A","B","B"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(res, paste(levels(cts$Strain)[f], "OG_nofilter_deseq.csv", sep = ""))
    resLFC <- lfcShrink(dds, coef=2, type="apeglm")
    write.csv(resLFC, paste(levels(cts$Strain)[f], "OG_nofilter_deseq_LFC.csv", sep = ""))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} 

### Deseq community OG ####
cts=read.csv("community.OG.qf", sep="\t", h=T, row.names = 1)
library("DESeq2")
library("apeglm")

    cts2=sapply(cts, as.integer)
    row.names(cts2)=row.names(cts)
    cond=factor(c("A","A","B","A","B","B"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(res, "community_OG_unfiltered_deseq.csv")
    resLFC <- lfcShrink(dds, coef=2, type="apeglm")
    write.csv(resLFC, "community_OG_unfiltered_deseq_LFC.csv")
    
### Deseq community OG with duplicates####
   cts=read.csv("community.OG_duplicates.qf", sep="\t", h=T, row.names = 1)
    library("DESeq2")
    library("apeglm")
    
    cts2=sapply(cts, as.integer)
    row.names(cts2)=row.names(cts)
    cond=factor(c("A","A","B","A","B","B"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(res, "community_OG_unfiltered_deseq.csv")
    resLFC <- lfcShrink(dds, coef=2, type="apeglm")
    write.csv(resLFC, "community_OG_unfiltered_deseq_LFC.csv")
    

    
#              #
#    FUNGI     #
#              #
    
    ### Deseq strain genes  ####
    cts=read.csv("strain.genes.qf", sep="\t", h=T, row.names = 1)
    library("DESeq2")
    library("apeglm")
    
    
    for (f in 1:length(levels(cts$Strain))) {
      tryCatch({
        counts=cts[which(cts$Strain==levels(cts$Strain)[f]),]
        counts=counts[,3:8]
        cts2=sapply(counts, as.integer)
        row.names(cts2)=row.names(counts)
        cond=factor(c("A","A","B","A","B","B"))
        dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
        dds <- estimateSizeFactors(dds)
        dds <- DESeq(dds)
        res <- results(dds)
        write.csv(res, paste(levels(cts$Strain)[f], "unfiltered_deseq.csv", sep = ""))
        resLFC <- lfcShrink(dds, coef=2, type="apeglm")
        write.csv(resLFC, paste(levels(cts$Strain)[f], "unfiltered_deseq_LFC.csv", sep = ""))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    } 
    
    
    
    ### Deseq strain OG  ####
    cts=read.csv("strain.OG.qf", sep="\t", h=T)
    
    for (f in 1:length(levels(cts$Strain))) {
      tryCatch({
        counts=cts[which(cts$Strain==levels(cts$Strain)[f]),]
        rownames(counts)=counts$OG
        counts=counts[,3:8]
        cts2=sapply(counts, as.integer)
        row.names(cts2)=row.names(counts)
        cond=factor(c("A","A","B","A","B","B"))
        dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
        dds <- estimateSizeFactors(dds)
        dds <- DESeq(dds)
        res <- results(dds)
        write.csv(res, paste(levels(cts$Strain)[f], "OG_nofilter_deseq.csv", sep = ""))
        resLFC <- lfcShrink(dds, coef=2, type="apeglm")
        write.csv(resLFC, paste(levels(cts$Strain)[f], "OG_nofilter_deseq_LFC.csv", sep = ""))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    } 
    
    ### Deseq community OG ####
    library("DESeq2")
    library("apeglm")
    
    cts2=sapply(cts, as.integer)
    row.names(cts2)=row.names(cts)
    cond=factor(c("A","A","B","A","B","B"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(res, "community_OG_unfiltered_deseq.csv")
    resLFC <- lfcShrink(dds, coef=2, type="apeglm")
    write.csv(resLFC, "community_OG_unfiltered_deseq_LFC.csv")
    
    ### Deseq community OG with duplicates ####
    setwd("E:/data_14-03-20/filesForDESeq202004/filesForDESeq_bacteria")
    cts=read.csv("community.OG.qf", sep="\t", h=T, row.names = 1)
    library("DESeq2")
    library("apeglm")
    
    cts2=sapply(cts, as.integer)
    row.names(cts2)=row.names(cts)
    cond=factor(c("A","A","B","A","B","B"))
    dds <- DESeqDataSetFromMatrix(cts2, DataFrame(cond), ~cond)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)

    
    