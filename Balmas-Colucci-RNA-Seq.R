#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# RNA-Seq Analysis to accompany:
# Balmas et al (2017)
# Group 2 innate lymphoid cells prevent endotoxin-induced fetal demise        
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2017-Balmas-Colucci
#
# CTR Code: CTR_fc287_0002 
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics Facility
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

#
# initial install of packages, uncomment and mofify as required
#
#source("http://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
#install.packages("ggrepel")
#install.packages("gplots")

library('DESeq2')
library('ggplot2')
library('ggplot2')
library("pheatmap")
library("ggrepel")
library("reshape2")
library("biomaRt")
library("matrixStats")
library("plyr")


message("+-------------------------------------------------------------------------------")
message("+ Set up some constants e.g. base directories")
message("+-------------------------------------------------------------------------------")

Project  <- "2017-Balmas-Colucci"
Base.dir <- getwd() 
setwd(Base.dir)
HTSeq.dir <- paste(Base.dir,"/HTSeq_Counts", sep="")

elementTextSize <- 10

significance    <- 0.05



message("+-------------------------------------------------------------------------------")
message("+ Set up the sample table")
message("+-------------------------------------------------------------------------------")

#
# Other samples = Lung & Lymph node
#
sampleFiles      <- grep('SLX-91*',list.files(HTSeq.dir),value=TRUE)
sampleNames      <- gsub(".r_1_val_tophat.accepted_hits_merged.bam_htseq_counts.txt", "", sampleFiles)
sampleNames      <- gsub(".C5.*.s_.", "", sampleNames)
sampleNames      <- gsub("L2DRBC", "", sampleNames)
sampleBarcodes   <- sampleNames
sampleBatch      <- gsub("\\..*", "", sampleNames)
sampleTissue     <- c("Lymph_Node", "Lung", "Lymph_Node", "Lung", "Lymph_Node", "Lung")#, 
sampleComp_1     <- c("Other", "Other", "Other", "Other", "Other", "Other")
sampleComp_2     <- c("Other", "Other", "Other", "Other", "Other", "Other")
sampleCell       <- c("ILC2", "ILC2", "ILC2", "ILC2", "ILC2", "ILC2")

sampleTableOther <- data.frame(sampleName=sampleBarcodes, fileName=sampleFiles, condition=sampleTissue, cell=sampleCell, batch=sampleBatch, comp1=sampleComp_1, comp2=sampleComp_2)
print(sampleTableOther)

#
# Uterine Samples
#
sampleFiles      <- grep('SLX-10*',list.files(HTSeq.dir),value=TRUE)
sampleNames      <- gsub(".C9.*.htseq_counts.txt", "", sampleFiles)
sampleNames      <- gsub("L2DRBC", "", sampleNames)
sampleBatch      <- gsub("\\..*", "", sampleNames)
sampleBarcodes   <- sampleNames
sampleGestation  <- c("Virgin", "Virgin", "Virgin", "E9.5", "E9.5", "E9.5", "E18.5", "E18.5", "E18.5")
sampleComp_1     <- c("Virgin", "Virgin", "Virgin", "Other", "Other", "Other", "Other", "Other", "Other")
sampleComp_2     <- c("Uterine", "Uterine", "Uterine",   "Uterine", "Uterine", "Uterine",   "Uterine", "Uterine", "Uterine")
sampleTimePoint  <- c("ILC2", "ILC2", "ILC2", "ILC2", "ILC2", "ILC2", "ILC2", "ILC2", "ILC2")

sampleTableUterine <- data.frame(sampleName=sampleBarcodes, fileName=sampleFiles, condition=sampleGestation, cell=sampleTimePoint, batch=sampleBatch,  comp1=sampleComp_1, comp2=sampleComp_2)
print(sampleTableUterine)


sampleTable <- rbind(sampleTableUterine,sampleTableOther)
head(sampleTable, 15)


message("+-------------------------------------------------------------------------------")
message("+ Retrieve ensEMBL annotations")
message("+-------------------------------------------------------------------------------")
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)          
head(ensEMBL2id)

message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")

ddsHTSeq         <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design=~condition)
colnames(ddsHTSeq)

#
# Uterine only samples for the heatmap
#
ddsHTSeq.Uterine <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTableUterine, directory=HTSeq.dir, design=~condition)
colnames(ddsHTSeq.Uterine)

message("+-------------------------------------------------------------------------------")
message("+ Create dds object")
message("+-------------------------------------------------------------------------------")

dds         <- DESeq(ddsHTSeq)
dds.Uterine <- DESeq(ddsHTSeq.Uterine)

message("+-------------------------------------------------------------------------------")
message("+ Run transformations")
message("+-------------------------------------------------------------------------------")

rld         <- rlogTransformation(dds, blind=T)
rld.Uterine <- rlogTransformation(dds.Uterine, blind=T)
vsd         <- varianceStabilizingTransformation(dds, blind=T)

message("+-------------------------------------------------------------------------------")
message("+ Get the normalised read counts")
message("+-------------------------------------------------------------------------------")

normCounts                 <- as.data.frame(counts(dds, normalized=TRUE))
normCounts$ensembl_gene_id <- rownames(normCounts)
normCounts.ann             <- merge(normCounts, ensEMBL2id, by="ensembl_gene_id")
write.csv(normCounts.ann, file = paste(Project, "_DESeq2_NormalisedCounts.csv", sep=""))

message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

# All genes
pca = prcomp(t(assay(rld)))

pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(sampleTable$fileName, pca$x, sampleTable$condition, sampleTable$sampleName)

pdf(paste(Project, "_DESeq2_Annotated_PCA_Fig3A.pdf", sep=""),width=10,height=7)
par(bg=NA)
ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sampleTable$condition)) ) ) +
  #geom_point(size = 8, shape=1, color='black' ) +
  geom_jitter(size = 7, alpha=.5, width=4, height=4) + # to show overlapping points add a small random shift
  #geom_text_repel(aes(label=sampleTable$sampleName)) +
  xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste(Project, " PCA (All Genes)", sep="")) +
  scale_colour_manual(name="Tissue", values = c("purple", "#1B9E77", "#F27314", "red3", "grey38", "blue")) +
  theme(text = element_text(size=elementTextSize)) 
dev.off()

message("+-------------------------------------------------------------------------------")
message("+ Create results object")
message("+-------------------------------------------------------------------------------")


res_Virgin_E9.5                       <- results(dds, contrast=c("condition",  "Virgin", "E9.5"))
res_Virgin_E9.5                       <- res_Virgin_E9.5[order(res_Virgin_E9.5$padj),]
res_Virgin_E9.5.df                    <- as.data.frame(subset(res_Virgin_E9.5, padj <= significance  & padj > 0))
res_Virgin_E9.5.df                    <- data.frame(ensembl_gene_id=rownames(res_Virgin_E9.5.df),res_Virgin_E9.5.df)
res_Virgin_E9.5.ann                   <- merge(res_Virgin_E9.5.df, ensEMBL2id, by="ensembl_gene_id")
res_Virgin_E9.5.ann$description       <- gsub("..Source.*", "", res_Virgin_E9.5.ann$description)
res_Virgin_E9.5.ann.ranked            <- res_Virgin_E9.5.ann[order(abs(res_Virgin_E9.5.ann$log2FoldChange), decreasing = TRUE),] 
head(res_Virgin_E9.5.ann.ranked)
write.csv(res_Virgin_E9.5.ann[order(abs(res_Virgin_E9.5.ann$log2FoldChange), decreasing = TRUE),], file = paste(Project, "_DESeq2_DEGs_", 'Virgin_vs_E9.5.csv', sep=""))


res_Virgin_E18.5                      <- results(dds, contrast=c("condition",  "Virgin", "E18.5"))
res_Virgin_E18.5                      <- res_Virgin_E18.5[order(res_Virgin_E18.5$padj),]
res_Virgin_E18.5.df                   <- as.data.frame(subset(res_Virgin_E18.5, padj < significance  & padj > 0))
res_Virgin_E18.5.df                   <- data.frame(ensembl_gene_id=rownames(res_Virgin_E18.5.df),res_Virgin_E18.5.df)
res_Virgin_E18.5.ann                  <- merge(res_Virgin_E18.5.df, ensEMBL2id, by="ensembl_gene_id")
res_Virgin_E18.5.ann$description      <- gsub("..Source.*", "", res_Virgin_E18.5.ann$description)
res_Virgin_E18.5.ann.ranked           <- res_Virgin_E18.5.ann[order(abs(res_Virgin_E18.5.ann$log2FoldChange), decreasing = TRUE),] 
head(res_Virgin_E18.5.ann.ranked)
write.csv(res_Virgin_E18.5.ann[order(abs(res_Virgin_E18.5.ann$log2FoldChange), decreasing = TRUE),], file = paste(Project, "_DESeq2_DEGs_", 'Virgin_vs_E18.5.csv', sep=""))


ddsHTSeq.comp2 <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design=~comp2)
dds.comp2      <- DESeq(ddsHTSeq.comp2)
rld.comp2      <- rlogTransformation(dds.comp2, blind=T)
vsd.comp2      <- varianceStabilizingTransformation(dds.comp2, blind=T)

res_Uterine_Others  <- results(dds.comp2) #, contrast=c("comp1",  "Virgin", "Others"))
res_Uterine_Others <- res_Uterine_Others[order(res_Uterine_Others$padj),]
res_Uterine_Others.df <- as.data.frame(subset(res_Uterine_Others, padj <= significance  & padj > 0))
res_Uterine_Others.df  <- data.frame(ensembl_gene_id=rownames(res_Uterine_Others.df),res_Uterine_Others.df)
res_Uterine_Others.ann <- merge(res_Uterine_Others.df, ensEMBL2id, by="ensembl_gene_id")
res_Uterine_Others.ann$description <- gsub("..Source.*", "", res_Uterine_Others.ann$description)
res_Uterine_Others.ann.ranked <- res_Uterine_Others.ann[order(abs(res_Uterine_Others.ann$log2FoldChange), decreasing = TRUE),] 
head(res_Uterine_Others.ann.ranked)
write.csv(res_Uterine_Others.ann[order(abs(res_Uterine_Others.ann$log2FoldChange), decreasing = TRUE),], file = paste(Project, "_DESeq2_DEGs_", 'Uterine_vs_Others.csv', sep=""))


message("+-------------------------------------------------------------------------------")
message("+ Create MA Plots")
message("+-------------------------------------------------------------------------------")

functionCustomEnrichmentPlot <- function(results, Project, Title, significance, foldchange, favoritegenes) {
  # Get annotation infor for a favorite set of ensEMBL ids
  results$log2FoldChange[results$log2FoldChange > 12.5] <- 12.5
  results$log2FoldChange[results$log2FoldChange < -12.5] <- -12.5
  
  labeldata.rows            <- match(favoritegenes, row.names(results))
  labeldata                 <- results[labeldata.rows,]
  labeldata$ensembl_gene_id <- rownames(labeldata)
  labeldata.ann             <- merge(labeldata, ensEMBL2id, by="ensembl_gene_id")
  
  plt <- ggplot(data = results, aes(x=baseMean, y=log2FoldChange )) + 
    geom_abline(intercept = logFoldChanceCutOff, slope = 0, colour='red', alpha=0.25) + 
    geom_abline(intercept = -logFoldChanceCutOff, slope = 0, colour='red', alpha=0.25) +
    geom_point(size=0.5, alpha=0.5, col="black") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange >= 0)), size=1, alpha=0.5,  col="blue") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange < 0)),  size=1, alpha=0.5,  col="red") +
    
    geom_point(data=subset(labeldata.ann, log2FoldChange > 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkblue',  shape=21, stroke=0.5) +
    geom_point(data=subset(labeldata.ann, log2FoldChange < 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkred', shape=21, stroke=0.5) +
    geom_label_repel(data=subset(labeldata.ann, log2FoldChange > 0), 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='blue', colour='white', point.padding = unit(0.25, "lines"),  size=6, segment.color = 'darkblue',  nudge_y = 1) +
    geom_label_repel(data=subset(labeldata.ann, log2FoldChange < 0), 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='red', colour='white', point.padding = unit(0.25, "lines"), size=6, segment.color = 'darkred', nudge_y = -1) +
    scale_x_log10() + scale_y_reverse() +
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") + ggtitle(paste(Project, " DESeq2 MA ", Title, " [fc ", logFoldChanceCutOff, ", sig ", significance, "]", sep="")) 
  
  pdf(paste(Project, "_DESeq2_MA", "_fc", logFoldChanceCutOff, "_sig", significance, "_", Title, ".pdf", sep=""),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt })
  dev.off()
}


logFoldChanceCutOff <- 5

favoritegenesY <- ensEMBL2id[grep("Gzm", ensEMBL2id$external_gene_name), ]$ensembl_gene_id #All Granzymes
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000026011") #CTLA4
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000026012") #CD28
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000071552") #Tigit
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000019256") #Ahr
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000029378") #Areg
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000014599") #Csf1 
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000036117") #Il5 
favoritegenesY <- append(favoritegenesY, "ENSMUSG00000019987") #Arg1

functionCustomEnrichmentPlot(as.data.frame(res_Uterine_Others), Project, "res_Uterine_Others_Fig3B", significance, logFoldChanceCutOff, favoritegenesY)


favoritegenes4B <- ""
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000079186")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000015441")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000040284")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000022156")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000059256")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000037202")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000045827")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000062345")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000037411")
favoritegenes4B <- append(favoritegenes4B, "ENSMUSG00000026715")

functionCustomEnrichmentPlot(as.data.frame(res_Virgin_E9.5),Project,  "res_Virgin_E9.5_Fig4B",  significance, logFoldChanceCutOff, favoritegenes4B)

favoritegenes4D <- ""
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000022156")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000079186")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000026011")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000026012")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000058427")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000035042")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000014599")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000038067")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000018916")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000025068")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000015522")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000019256")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000039697")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000002603")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000039239")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000029378")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000026335")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000027750")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000034708")
favoritegenes4D <- append(favoritegenes4D, "ENSMUSG00000019768")

functionCustomEnrichmentPlot(as.data.frame(res_Virgin_E18.5),Project, "res_Virgin_E18.5_Fig4D", significance, logFoldChanceCutOff, favoritegenes4D)

message("+-------------------------------------------------------------------------------")
message("+ Create HeatMaps")
message("+-------------------------------------------------------------------------------")

# Uterine specific gene lists

res_Virgin_E9.5  <- results(dds.Uterine, contrast=c("condition",  "Virgin", "E9.5"))
res_Virgin_E9.5  <- res_Virgin_E9.5[order(res_Virgin_E9.5$padj),]

res_Virgin_E18.5 <- results(dds.Uterine, contrast=c("condition",  "Virgin", "E18.5"))
res_Virgin_E18.5 <- res_Virgin_E18.5[order(res_Virgin_E18.5$padj),]

res_E9.5_E18.5   <- results(dds.Uterine, contrast=c("condition",  "E9.5", "E18.5"))
res_E9.5_E18.5   <- res_E9.5_E18.5[order(res_E9.5_E18.5$padj),]


# Get some meta data for annotating the pheatmap
df <- as.data.frame(colData(dds.Uterine)[,c("condition")])
colnames(df) <- "Gestation"
rownames(df) <- rownames(colData(dds.Uterine))
rownames(df) <- gsub("SLX-", "SLX.", rownames(df) )
head(df, 9)

# Take all above log2FoldChange cut off
logFoldChanceCutOff <- 7.5

l1  <- res_Virgin_E18.5[ order( - abs(res_Virgin_E18.5$log2FoldChange) ), ] 
l2  <- res_Virgin_E9.5[  order( - abs(res_Virgin_E9.5$log2FoldChange)  ), ] 
l3  <- res_E9.5_E18.5[   order( - abs(res_E9.5_E18.5$log2FoldChange)   ), ] 

l1l <- rownames( l1[ !(is.na(l1$log2FoldChange)) & (abs(l1$log2FoldChange) >= logFoldChanceCutOff) & !(is.na(l1$padj)) & (l1$padj <= significance), ] )
l2l <- rownames( l2[ !(is.na(l2$log2FoldChange)) & (abs(l2$log2FoldChange) >= logFoldChanceCutOff) & !(is.na(l2$padj)) & (l2$padj <= significance), ] )
l3l <- rownames( l3[ !(is.na(l3$log2FoldChange)) & (abs(l3$log2FoldChange) >= logFoldChanceCutOff) & !(is.na(l3$padj)) & (l3$padj <= significance), ] )

genes2plot <- unique( c(l1l, l2l, l3l) )

# Add external gene names to add to the plot and tidy up the table for plotting
rows               <- match(genes2plot, row.names(rld.Uterine))
mat                <- assay(rld.Uterine)[rows,c(1,2,3,4,5,6,7,8,9)]
mat2               <- mat
mat2.df            <- data.frame(ensembl_gene_id=rownames(mat2),mat2)
mat2.ann           <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
rownames(mat2.ann) <- mat2.ann$external_gene_name
mat2.new           <- mat2.ann[, -c(1,11,12) ]

pdf(paste(Project, "_DESeq2_CountMatrixHeatmap_topDEGs_lf", logFoldChanceCutOff, "sig", significance, "_Fig4A.pdf", sep=""), onefile=FALSE, width=10, height=7) 
par(bg=NA)
pheatmap(mat2.new, fontsize=5, annotation_col=df, fontsize_row=3, main=paste(Project," ::: DEGs log2FC >= ", logFoldChanceCutOff, " sign=", significance, sep=""))
dev.off()

message("+-------------------------------------------------------------------------------")
message("+ DESeq2 Analysis Complete")
message("+-------------------------------------------------------------------------------")
#------------------------------------------------------------------------------
# FIN
#------------------------------------------------------------------------------
