
library(EnhancedVolcano)
library(tidyverse)
setwd("all_genes/")
strains=c("Root61","Root322","Root22","Root102","Root172","Root149","Root491","Root685","Root123D2","Root670","Root154","Root1277","Root569","Root68","Root1280","Root70","Root179","Root559","Root667","Root181")
strains=as.factor(strains)
for (i in 1:length(levels(strains))){
  tryCatch({
    dat=read.csv(paste(levels(strains)[i], "unfiltered_deseq_LFC.csv",sep=''))
    
minimum=min(dat$log2FoldChange, na.rm=T)
maximum=max(dat$log2FoldChange, na.rm=T)

#create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals.colour <- rep('black', nrow(dat))

# set the base name/label as 'Mid'
names(keyvals.colour) <- rep('Not enriched', nrow(dat))

# modify keyvals for variables with fold change > 2.5
keyvals.colour[which(dat$log2FoldChange > 1.5 & dat$padj<0.05)] <- 'forestgreen'
names(keyvals.colour)[which(dat$log2FoldChange > 1.5 & dat$padj<0.05)] <- 'Enriched in Roots'

# modify keyvals for variables with fold change < -2.5
keyvals.colour[which(dat$log2FoldChange < -1.5)] <- 'darkorange3'
names(keyvals.colour)[which(dat$log2FoldChange < -1.5 & dat$padj<0.05)] <- 'Enriched in Matrix'




pdf(file = paste("F:/figures_mutant_papers/figures_23-09-21/volcano plots/", levels(strains)[20],".pdf" ,sep=''),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

EnhancedVolcano(dat,lab=rownames(dat),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(minimum-1, maximum+1),
                ylim = c(0,26),
                title = paste(levels(strains)[i]) ,
                subtitle = "Matrix vs Roots",
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.5,
                labSize = 0,
                labCol = "white",
                cutoffLineType = "blank",
                vline=c(-1.5, 1.5),
                caption = "log2FC cutoff, 1.5; p-value cutoff, 0.05",
                legendPosition = "right",
                legendLabSize = 14,
                colAlpha = 0.5, drawConnectors = FALSE,
                colCustom = keyvals.colour,hline = c(0.05))
dev.off()

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
