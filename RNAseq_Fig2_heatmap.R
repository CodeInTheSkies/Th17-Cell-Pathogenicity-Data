# This code creates the heatmap seen in Fig. 2, 
# using cuffdiff expression output lists and FPKMs
# The two input files needed, "cuffdiff_TG-GTP_gene_exp.diff" and "genes.read_group_tracking",
# can be downloaded from the repo folder "RNAseq_inputs"

library(ggplot2)
library(cowplot)
library(pheatmap)
library(ggforce)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

rm(list=ls())
strng <- "T_GFP_PosVsNeg" ; lpath <- "./cuffdiff_TG-GTP/"
spath <- "./heatmaps/"
df <- read.delim(paste0(lpath,"cuffdiff_TG-GTP_gene_exp.diff"), quote="")
fpkms <- read.delim(paste0(lpath,"genes.read_group_tracking"), quote="")
numOfGenes <- 20

# Getting the names of the top X and bottom X, plus Chat and Bhlhe40
# subsetting, and cleaning up, ordering by descending LogFC
subdf <- df %>% 
  filter(significant == "yes", log2.fold_change.!="Inf", log2.fold_change.!="-Inf") %>% 
  arrange(desc(log2.fold_change.))

# selecting top and bottom 20, after ordering by decreasing logFC
slctdGenes <- subdf$gene[c(1:numOfGenes, (length(subdf$gene)-numOfGenes+1):length(subdf$gene))]
# Checking, and adding Chat and Bhlhe40, if not already there
splitNum=0
if(!("Chat" %in% slctdGenes)){
  slctdGenes <- c(slctdGenes, "Chat")
  splitNum <- splitNum + 1
}
if(!("Bhlhe40" %in% slctdGenes)){
  slctdGenes <- c(slctdGenes, "Bhlhe40")
  splitNum <- splitNum + 1
}
# Now, slctdGenes has the needed gene names
# Since the gene_exp.diff file did not have the FPKMs of the individual replicates, we have to load the other file and extract
# the individual FPKMs from there
# Odd format, the FPKMs there are in rows, rather than columns
# We do the required wrangling next
tt=fpkms[fpkms$tracking_id %in% slctdGenes,]
tt$sample_id <- paste0(tt$condition,"_",tt$replicate)
tt$gene <- tt$tracking_id
subtt <- tt %>% select(gene, FPKM, sample_id)
subtt <- reshape2::dcast(subtt,gene~sample_id, value.var = "FPKM")
row.names(subtt) <- subtt$gene
subtt <- subtt[ , !(names(subtt) %in% "gene")]
colnames(subtt) <- gsub("GTP","GFP",colnames(subtt))
# order in the same order as slctdGenes, since that was top and bottom 20 PLUS Chat and Bhlhe40 at the end
subtt <- subtt[slctdGenes,]

# Heatmap
subtt <-t(scale(t(subtt)))
splitvec <- c(rep(1,2*numOfGenes),rep(2,splitNum))
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
hh <- Heatmap(subtt, name = "Expression heatmap", col = col_fun, rect_gp = gpar(col= "white"), heatmap_legend_param = list(title="FPKM \n (z-scored)"), border=T,
              row_title = "", row_split = splitvec, row_gap = unit(3, "mm"),
              cluster_columns = F, cluster_rows = F, show_row_dend = F, column_names_rot = 30, column_names_side=c("bottom", "top"), column_names_gp = gpar(fontsize = 11) )
pdf(file = paste(spath,"ExprHtMap_",strng,"_TopBtm",numOfGenes,"_v1.pdf",sep=""),
    useDingbats = FALSE, width=4, height=10)
draw(hh, padding = unit(c(2, 5, 2, 2), "mm")) #bottom, left, top, right paddings
dev.off()

