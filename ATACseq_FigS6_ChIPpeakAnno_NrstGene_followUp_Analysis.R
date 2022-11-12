# This script should be run using R version 4.0.3 or above
#
# Uses ChIPpeakAnno package to annotate ATACseq peaks
#
# To get the most informative peaks out, we do the following. From comparing Pos
# vs Neg cases for the three pairs, we leave the big chunk of shared genes out
# in each case. Instead, we take the smaller sets of genes that stand out from
# all three Neg fractions separately and perform an overlap Venn diagram for
# those to answer the question, "which of the Neg-specific genes are shared
# across the three sample pairs?" Then perform the same for all three
# Pos-specific genes.
#
# This will give us a smaller number of genes that may be highly correlated with
# the Neg or Pos cell state.
#
# INPUTs: Six "overlap.optimal_peak.narrowPeak" files corresponding to Pos and
# Neg of each of the three cases, namely, CD3CD28, Soluble, and SodiumButyrate
# These files can be downloaded from the repo folder "ATACseq_inputs" in ".gz" format, which need to be 
# uncompressed before running this script
# 
# OUTPUTs: Five CSV files, (available in the repo) self-explanatory from their names:
# "annoTaTeD_NrstGeneTSS_All3_Pos_3CircleVenn_22481_UniqSoluble.csv"
# "annoTaTeD_NrstGeneTSS_All3_Pos_3CircleVenn_4309_UniqSodiumButyrate.csv"
# "annoTaTeD_NrstGeneTSS_All3_Pos_3CircleVenn_4222_UniqCD3CD28.csv"
# "annoTaTeD_NrstGeneTSS_All3_Pos_UniquePeaksFromPosvsNegOverlappingRun.csv"
# "annoTaTeD_NrstGeneTSS_All3_Neg_UniquePeaksFromPosvsNegOverlappingRun.csv"
# 
# OUTPUTs: Two Venn diagrams, self-explanatory from their names:
# "All3_Pos_UniquePeaksFromPosvsNegOverlappingRun_Venn"
# "All3_Neg_UniquePeaksFromPosvsNegOverlappingRun_Venn"
# 
# All the above OUTPUT files can be found in the repo folder "ATACseq_outputs".

library(ChIPpeakAnno)

# We take the peak output files from ENCODE pipeline
# Why does MACS output duplicate peak regions?  Due to different summits within the same peak
# In the below discussions, they say that you can just merge those duplicate peaks
# https://www.biostars.org/p/464618/
# https://www.biostars.org/p/361109/
# http://seqanswers.com/forums/showthread.php?t=50394
lpath <- "./croo_outputs/"
spath <- "./results/"
# The scripts expects the above two folders to be already present; modify these
# paths accordingly as needed according to local computer requirements

# Pos vs Neg in pairs
# CD3CD28
DatRead <- read.table(file=paste0(lpath,"CD3CD28_Pos_overlap.optimal_peak.narrowPeak"), sep="\t", skip=1)
tt=data.frame(Chr=DatRead[,1], Start=DatRead[,2], End=DatRead[,3])
CD3CD28Pos <- toGRanges(unique(tt))
cat("\nNo. of unique peaks in CD3CD28 Pos: ", dim(tt)[1],"\n")
# Neg
DatRead <- read.table(file=paste0(lpath,"CD3CD28_Neg_overlap.optimal_peak.narrowPeak"), sep="\t", skip=1)
tt=data.frame(Chr=DatRead[,1], Start=DatRead[,2], End=DatRead[,3])
CD3CD28Neg <- toGRanges(unique(tt))
cat("\nNo. of unique peaks in CD3CD28 Neg: ", dim(tt)[1],"\n")
ol <- findOverlapsOfPeaks(CD3CD28Pos, CD3CD28Neg, maxgap=1000, connectedPeaks="keepAll")
CD3CD28_ol_unqPks <- ol$uniquePeaks

# Soluble
DatRead <- read.table(file=paste0(lpath,"Soluble_Pos_overlap.optimal_peak.narrowPeak"), sep="\t", skip=1)
tt=data.frame(Chr=DatRead[,1], Start=DatRead[,2], End=DatRead[,3])
SolublePos <- toGRanges(unique(tt))
cat("\nNo. of unique peaks in Soluble Pos: ", dim(tt)[1],"\n")
# Neg
DatRead <- read.table(file=paste0(lpath,"Soluble_Neg_overlap.optimal_peak.narrowPeak"), sep="\t", skip=1)
tt=data.frame(Chr=DatRead[,1], Start=DatRead[,2], End=DatRead[,3])
SolubleNeg <- toGRanges(unique(tt))
cat("\nNo. of unique peaks in Soluble Neg: ", dim(tt)[1],"\n")
ol <- findOverlapsOfPeaks(SolublePos, SolubleNeg, maxgap=1000, connectedPeaks="keepAll")
Soluble_ol_unqPks <- ol$uniquePeaks

# SodiumButyrate
DatRead <- read.table(file=paste0(lpath,"SodiumButyrate_Pos_overlap.optimal_peak.narrowPeak"), sep="\t", skip=1)
tt=data.frame(Chr=DatRead[,1], Start=DatRead[,2], End=DatRead[,3])
SodiumButyratePos <- toGRanges(unique(tt))
cat("\nNo. of unique peaks in SodiumButyrate Pos: ", dim(tt)[1],"\n")
# Neg
DatRead <- read.table(file=paste0(lpath,"SodiumButyrate_Neg_overlap.optimal_peak.narrowPeak"), sep="\t", skip=1)
tt=data.frame(Chr=DatRead[,1], Start=DatRead[,2], End=DatRead[,3])
SodiumButyrateNeg <- toGRanges(unique(tt))
cat("\nNo. of unique peaks in SodiumButyrate Neg: ", dim(tt)[1],"\n")
ol <- findOverlapsOfPeaks(SodiumButyratePos, SodiumButyrateNeg, maxgap=1000, connectedPeaks="keepAll")
SodiumButyrate_ol_unqPks <- ol$uniquePeaks

# Now, finding overlap peaks between pairs of unique_peaks GRanges objects
# Just subsetting the GRanges objects
CD3CD28Pos_unqPks <- CD3CD28_ol_unqPks[grep(pattern = "^CD3CD28Pos", names(CD3CD28_ol_unqPks))]
CD3CD28Neg_unqPks <- CD3CD28_ol_unqPks[grep(pattern = "^CD3CD28Neg", names(CD3CD28_ol_unqPks))]

SolublePos_unqPks <- Soluble_ol_unqPks[grep(pattern = "^SolublePos", names(Soluble_ol_unqPks))]
SolubleNeg_unqPks <- Soluble_ol_unqPks[grep(pattern = "^SolubleNeg", names(Soluble_ol_unqPks))]

SodiumButyratePos_unqPks <- SodiumButyrate_ol_unqPks[grep(pattern = "^SodiumButyratePos", names(SodiumButyrate_ol_unqPks))]
SodiumButyrateNeg_unqPks <- SodiumButyrate_ol_unqPks[grep(pattern = "^SodiumButyrateNeg", names(SodiumButyrate_ol_unqPks))]

length(CD3CD28Pos_unqPks)
length(CD3CD28Neg_unqPks)

length(SolublePos_unqPks)
length(SolubleNeg_unqPks)

length(SodiumButyratePos_unqPks)
length(SodiumButyrateNeg_unqPks)
# The above lengths match
# > length(CD3CD28Pos_unqPks)
# 5827
# > length(CD3CD28Neg_unqPks)
# 14797
# >
# > length(SolublePos_unqPks)
# 24984
# > length(SolubleNeg_unqPks)
# 4999
# >
# > length(SodiumButyratePos_unqPks)
# 6361
# > length(SodiumButyrateNeg_unqPks)
# 11606

# Now, find the overlaps and try the Venn diagrams
ol_all3_Negs <- findOverlapsOfPeaks(CD3CD28Neg_unqPks, SolubleNeg_unqPks, SodiumButyrateNeg_unqPks, maxgap=1000, connectedPeaks="keepAll")
all3SharedNeg <- ol_all3_Negs$peaklist$`CD3CD28Neg_unqPks///SolubleNeg_unqPks///SodiumButyrateNeg_unqPks`
pdf(file = paste0(spath,"All3_Neg_UniquePeaksFromPosvsNegOverlappingRun_Venn.pdf"),
    useDingbats = FALSE, width=8, height=6)
makeVennDiagram(ol_all3_Negs, totalTest=1e6,
                fill=c("#CC79A7", "#56B4E9", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2", "#E69F00"), #circle border color
                cat.col=c("#D55E00", "#0072B2", "#E69F00"))
dev.off()


ol_all3_Pos <- findOverlapsOfPeaks(CD3CD28Pos_unqPks, SolublePos_unqPks, SodiumButyratePos_unqPks, maxgap=1000, connectedPeaks="keepAll")
all3SharedPos <- ol_all3_Pos$peaklist$`CD3CD28Pos_unqPks///SolublePos_unqPks///SodiumButyratePos_unqPks`
pdf(file = paste0(spath,"All3_Pos_UniquePeaksFromPosvsNegOverlappingRun_Venn.pdf"),
    useDingbats = FALSE, width=8, height=6)
makeVennDiagram(ol_all3_Pos, totalTest=1e6,
                fill=c("#CC79A7", "#56B4E9", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2", "#E69F00"), #circle border color
                cat.col=c("#D55E00", "#0072B2", "#E69F00"))
dev.off()

#
# Now, we try and add annotations to just those Ranges/genes that are shared between all 3 Unique sets
# All3 Neg
data(TSS.mouse.GRCm38)
wrkngDat.annotD <- annotatePeakInBatch(all3SharedNeg, AnnotationData=TSS.mouse.GRCm38,
                                       output="nearestLocation", select="arbitrary")
wrkngDat.annotD_DF <- as.data.frame(wrkngDat.annotD)

library('biomaRt')
mart <- useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl", host = "https://dec2021.archive.ensembl.org"))
genestt <- unique(as.vector(wrkngDat.annotD_DF$feature))
G_list <- biomaRt::getBM(attributes= c("ensembl_gene_id","mgi_symbol"), filters="ensembl_gene_id", values=genestt, mart=mart)
wrkngDat.annotD_DF_new <- merge(wrkngDat.annotD_DF,G_list,by.x="feature",by.y="ensembl_gene_id",all.x=T,sort=F)

# Renaming some columns for consistency before merging!
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "seqnames")] <- "Chr"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "start")] <- "Start"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "end")] <- "End"

colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "mgi_symbol")] <- "NearestGene_mgi_symbol"
# moving mgi_symbol to front using dplyr
library(dplyr)
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, NearestGene_mgi_symbol, everything())
wrkngDat.annotD_DF_new$peakNames <- gsub(",",";",wrkngDat.annotD_DF_new$peakNames, fixed=T)
# moving the peakNames column to the end
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, -peakNames, peakNames)
write.csv(wrkngDat.annotD_DF_new, file=paste0(spath, "annoTaTeD_NrstGeneTSS_All3_Neg_UniquePeaksFromPosvsNegOverlappingRun.csv"))

# All3 Pos
data(TSS.mouse.GRCm38)
wrkngDat.annotD <- annotatePeakInBatch(all3SharedPos, AnnotationData=TSS.mouse.GRCm38,
                                       output="nearestLocation", select="arbitrary")
wrkngDat.annotD_DF <- as.data.frame(wrkngDat.annotD)

mart <- useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl", host = "https://dec2021.archive.ensembl.org"))
genestt <- unique(as.vector(wrkngDat.annotD_DF$feature))
G_list <- biomaRt::getBM(attributes= c("ensembl_gene_id","mgi_symbol"), filters="ensembl_gene_id", values=genestt, mart=mart)
wrkngDat.annotD_DF_new <- merge(wrkngDat.annotD_DF,G_list,by.x="feature",by.y="ensembl_gene_id",all.x=T,sort=F)

# Renaming some columns for consistency before merging!
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "seqnames")] <- "Chr"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "start")] <- "Start"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "end")] <- "End"

colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "mgi_symbol")] <- "NearestGene_mgi_symbol"
# moving mgi_symbol to front using dplyr
library(dplyr)
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, NearestGene_mgi_symbol, everything())
wrkngDat.annotD_DF_new$peakNames <- gsub(",",";",wrkngDat.annotD_DF_new$peakNames, fixed=T)
# moving the peakNames column to the end
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, -peakNames, peakNames)
write.csv(wrkngDat.annotD_DF_new, file=paste0(spath, "annoTaTeD_NrstGeneTSS_All3_Pos_UniquePeaksFromPosvsNegOverlappingRun.csv"))

# -------------
# Now, we try and add annotations to those Ranges/genes that are Unique to each of the 3 Unique sets
# i.e., the numbers 4222, 22481, and 4309 in the Pos Venn 3-circle

data(TSS.mouse.GRCm38)
# All3 Pos - CD3CD28Pos_unqPks
wrkngDat.annotD <- annotatePeakInBatch(ol_all3_Pos$peaklist$CD3CD28Pos_unqPks, AnnotationData=TSS.mouse.GRCm38,
                                       output="nearestLocation", select="arbitrary")
wrkngDat.annotD_DF <- as.data.frame(wrkngDat.annotD)

mart <- useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl", host = "https://dec2021.archive.ensembl.org"))
genestt <- unique(as.vector(wrkngDat.annotD_DF$feature))
G_list <- biomaRt::getBM(attributes= c("ensembl_gene_id","mgi_symbol"), filters="ensembl_gene_id", values=genestt, mart=mart)
wrkngDat.annotD_DF_new <- merge(wrkngDat.annotD_DF,G_list,by.x="feature",by.y="ensembl_gene_id",all.x=T,sort=F)

# Renaming some columns for consistency before merging!
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "seqnames")] <- "Chr"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "start")] <- "Start"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "end")] <- "End"

colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "mgi_symbol")] <- "NearestGene_mgi_symbol"
# moving mgi_symbol to front using dplyr
library(dplyr)
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, NearestGene_mgi_symbol, everything())
wrkngDat.annotD_DF_new$peakNames <- gsub(",",";",wrkngDat.annotD_DF_new$peakNames, fixed=T)
# moving the peakNames column to the end
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, -peakNames, peakNames)
write.csv(wrkngDat.annotD_DF_new, file=paste0(spath, "annoTaTeD_NrstGeneTSS_All3_Pos_3CircleVenn_4222_UniqCD3CD28.csv"))

# All3 Pos - SodiumButyratePos_unqPks
wrkngDat.annotD <- annotatePeakInBatch(ol_all3_Pos$peaklist$SodiumButyratePos_unqPks, AnnotationData=TSS.mouse.GRCm38,
                                       output="nearestLocation", select="arbitrary")
wrkngDat.annotD_DF <- as.data.frame(wrkngDat.annotD)

mart <- useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl", host = "https://dec2021.archive.ensembl.org"))
genestt <- unique(as.vector(wrkngDat.annotD_DF$feature))
G_list <- biomaRt::getBM(attributes= c("ensembl_gene_id","mgi_symbol"), filters="ensembl_gene_id", values=genestt, mart=mart)
wrkngDat.annotD_DF_new <- merge(wrkngDat.annotD_DF,G_list,by.x="feature",by.y="ensembl_gene_id",all.x=T,sort=F)

# Renaming some columns for consistency before merging!
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "seqnames")] <- "Chr"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "start")] <- "Start"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "end")] <- "End"

colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "mgi_symbol")] <- "NearestGene_mgi_symbol"
# moving mgi_symbol to front using dplyr
library(dplyr)
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, NearestGene_mgi_symbol, everything())
wrkngDat.annotD_DF_new$peakNames <- gsub(",",";",wrkngDat.annotD_DF_new$peakNames, fixed=T)
# moving the peakNames column to the end
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, -peakNames, peakNames)
write.csv(wrkngDat.annotD_DF_new, file=paste0(spath, "annoTaTeD_NrstGeneTSS_All3_Pos_3CircleVenn_4309_UniqSodiumButyrate.csv"))


# All3 Pos - SolublePos_unqPks
wrkngDat.annotD <- annotatePeakInBatch(ol_all3_Pos$peaklist$SolublePos_unqPks, AnnotationData=TSS.mouse.GRCm38,
                                       output="nearestLocation", select="arbitrary")
wrkngDat.annotD_DF <- as.data.frame(wrkngDat.annotD)

mart <- useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl", host = "https://dec2021.archive.ensembl.org"))
genestt <- unique(as.vector(wrkngDat.annotD_DF$feature))
G_list <- biomaRt::getBM(attributes= c("ensembl_gene_id","mgi_symbol"), filters="ensembl_gene_id", values=genestt, mart=mart)
wrkngDat.annotD_DF_new <- merge(wrkngDat.annotD_DF,G_list,by.x="feature",by.y="ensembl_gene_id",all.x=T,sort=F)

# Renaming some columns for consistency before merging!
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "seqnames")] <- "Chr"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "start")] <- "Start"
colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "end")] <- "End"

colnames(wrkngDat.annotD_DF_new)[which(names(wrkngDat.annotD_DF_new) == "mgi_symbol")] <- "NearestGene_mgi_symbol"
# moving mgi_symbol to front using dplyr
library(dplyr)
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, NearestGene_mgi_symbol, everything())
wrkngDat.annotD_DF_new$peakNames <- gsub(",",";",wrkngDat.annotD_DF_new$peakNames, fixed=T)
# moving the peakNames column to the end
wrkngDat.annotD_DF_new <- dplyr::select(wrkngDat.annotD_DF_new, -peakNames, peakNames)
write.csv(wrkngDat.annotD_DF_new, file=paste0(spath, "annoTaTeD_NrstGeneTSS_All3_Pos_3CircleVenn_22481_UniqSoluble.csv"))

# In the above 3 csv files, the number of annotated rows is a bit less than the original starting numbers, 4309, etc.  I think this is because of 
# some overlapping going on, discussed here:
# https://support.bioconductor.org/p/63583/


# > sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.0.6          biomaRt_2.44.4       ChIPpeakAnno_3.22.4  GenomicRanges_1.40.0 GenomeInfoDb_1.24.2  Biostrings_2.56.0    XVector_0.28.0       IRanges_2.22.2      
# [9] S4Vectors_0.26.1     BiocGenerics_0.34.0 
# 
# loaded via a namespace (and not attached):
#   [1] ProtGenerics_1.20.0         bitops_1.0-6                matrixStats_0.57.0          bit64_4.0.5                 RColorBrewer_1.1-2          progress_1.2.2             
# [7] httr_1.4.2                  tools_4.0.3                 utf8_1.1.4                  R6_2.5.0                    lazyeval_0.2.2              DBI_1.1.0                  
# [13] colorspace_2.0-0            ade4_1.7-16                 GetoptLong_1.0.5            withr_2.3.0                 tidyselect_1.1.0            prettyunits_1.1.1          
# [19] VennDiagram_1.6.20          bit_4.0.4                   curl_4.3                    compiler_4.0.3              graph_1.66.0                Biobase_2.48.0             
# [25] formatR_1.7                 xml2_1.3.2                  DelayedArray_0.14.1         rtracklayer_1.48.0          RBGL_1.64.0                 askpass_1.1                
# [31] rappdirs_0.3.1              stringr_1.4.0               digest_0.6.27               Rsamtools_2.4.0             pkgconfig_2.0.3             ensembldb_2.12.1           
# [37] dbplyr_2.1.1                limma_3.44.3                BSgenome_1.56.0             regioneR_1.20.1             rlang_0.4.11                GlobalOptions_0.1.2        
# [43] rstudioapi_0.13             RSQLite_2.2.1               shape_1.4.5                 generics_0.1.0              BiocParallel_1.22.0         RCurl_1.98-1.2             
# [49] magrittr_2.0.1              GO.db_3.11.4                GenomeInfoDbData_1.2.3      futile.logger_1.4.3         Matrix_1.3-0                Rcpp_1.0.5                 
# [55] fansi_0.4.1                 lifecycle_1.0.0             stringi_1.5.3               MASS_7.3-53                 SummarizedExperiment_1.18.2 zlibbioc_1.34.0            
# [61] BiocFileCache_1.12.1        grid_4.0.3                  blob_1.2.1                  crayon_1.4.1                lattice_0.20-41             splines_4.0.3              
# [67] GenomicFeatures_1.40.1      multtest_2.44.0             circlize_0.4.11             hms_1.0.0                   ComplexHeatmap_2.4.3        pillar_1.6.0               
# [73] rjson_0.2.20                seqinr_4.2-5                futile.options_1.0.1        XML_3.99-0.5                glue_1.4.2                  lambda.r_1.2.4             
# [79] BiocManager_1.30.10         idr_1.2                     png_0.1-7                   vctrs_0.3.6                 openssl_1.4.3               purrr_0.3.4                
# [85] clue_0.3-58                 assertthat_0.2.1            xfun_0.22                   AnnotationFilter_1.12.0     survival_3.2-7              tibble_3.1.1               
# [91] GenomicAlignments_1.24.0    tinytex_0.31                AnnotationDbi_1.50.3        memoise_1.1.0               cluster_2.1.0               ellipsis_0.3.1 