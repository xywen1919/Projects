# load packages

library(GEOquery)
library(GOstats)
library(GO.db)
library(Category)

library(AnnotationDbi) # retrieve info from dataset
library(annotate) 
library(org.Hs.eg.db)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ggridges)
# BiocManager::install("ggupset")
library(ggupset)

library(dendextend)
library(vsn) # for visualization
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(ggbeeswarm)
library(apeglm)

library(PoiClaClu)
library(glmpca)
library(M3C)
library(Rtsne)
library(cluster)
# install_github("vqv/ggbiplot")
library(ggbiplot)
# install.packages("gridExtra")
library(gridExtra)
library(devtools)
library(rgl)

# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
# install.packages("ggpubr")
library(ggpubr)

# =====================================================
# Data prepare
# =====================================================
# reference package vignettes
# browseVignettes(package = "DESeq2")

# load the count data
fileURL <- paste("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE153921&format=file&file=GSE153921%5FAndrew%5FXPO7%5Fmerged%5Fgene%5Fcounts%2Ecsv%2Egz")

download.file(fileURL,"GSE153921_Andrew_XPO7_merged_gene_counts.csv.gz")

gene_counts <- read.csv("GSE153921_Andrew_XPO7_merged_gene_counts.csv.gz", row.names=1)

dim(gene_counts)

# filtering out zero counts
keep <- apply(gene_counts, 1, sum) > 10
gene_counts <- gene_counts[keep,]

dim(gene_counts)


# obtain corresponding gene symbol from ENSEMBL IDs (rownames)
keytypes(org.Hs.eg.db) 
 
# create a data frame for gene annotation, with ensembl id as row name
gene_anno <- data.frame(symbol = mapIds(org.Hs.eg.db, rownames(gene_counts),
                                        keytype = "ENSEMBL",
                                        column = "SYMBOL"
)
)                          

rownames(gene_anno) <- rownames(gene_counts)

# add genename column
gene_anno$geneName <- mapIds(org.Hs.eg.db, rownames(gene_counts),
                              keytype = "ENSEMBL", 
                              column = "GENENAME")
# add entrez id column
gene_anno$ENTREZID <- mapIds(org.Hs.eg.db, rownames(gene_counts),
                             keytype = "ENSEMBL", 
                             column = "ENTREZID") 
head(gene_anno) 

keytypes(org.Hs.eg.db)

# update gene_counts and gene_anno with real symbol
keep_anno <- which(gene_anno$symbol != "NA")

gene_counts <- gene_counts[keep_anno,]
head(gene_counts)

gene_anno <- gene_anno[keep_anno,]
head(gene_anno)


# study design
gse <- getGEO("GSE153921")
gsedta <- gse$GSE153921_series_matrix.txt.gz


coldata <- data.frame(Expgroup = gsedta$source_name_ch1)
rownames(coldata) <- names(gene_counts)
coldata$GEO <- gsedta$geo_accession
coldata$se <- relevel(factor(rep(c("non_induced", "induced"), 
                                 c(3, 9))), 
                      "non_induced")
coldata$ko <- relevel(factor(rep(c("control", "no_ko", "ko_shXPO7_3","ko_shXPO7_5"), 
                                 c(3,3,3,3))),
                      "control")



#creating the DESeq2 object from the matrix of counts 
dds <- DESeqDataSetFromMatrix(countData = gene_counts, 
                              colData = coldata,
                              design = ~ ko)

colData(dds)

# make a copy for dds
cnts <- dds

# ========================================================================
# Normalization of dataset
# before normalization, perform size factor estimation
cnts <- estimateSizeFactors(cnts)
cnts <- estimateDispersions(cnts)

plotDispEsts(cnts)


# three normalization methods: vst, rlog, and log2
# ntd
ntd_cnts <- normTransform(cnts)
# head(assay(ntd_cnts))

# vst
vst_cnts <- vst(cnts, blind = FALSE)
# head(assay(vst_cnts))

# rlog
rld_cnts <- rlog(cnts, blind = FALSE)
# head(assay(rld_cnts))


# difference between normalized and non-normalized read counts
par(mfrow=c(2,2))
boxplot(log2(counts(cnts)+1), notch=TRUE,
        main="log2(Non-normalized read counts+1)",
        cex=.6, xaxt="n")
boxplot(assay(ntd_cnts), notch=TRUE,
        main="Normal transformed read counts", 
        cex=.6, xaxt="n")
boxplot(assay(vst_cnts), notch=TRUE,
        main="vst transformed read counts", 
        cex=.6, xaxt="n")
boxplot(assay(rld_cnts), notch=TRUE,
        main="rlog transformed read counts", 
        cex=.6, xaxt="n")


# plot the expression measurements for all 3 methods
par(mfrow=c(1,4))
plot(counts(cnts), main="non_normalized")
plot(assay(ntd_cnts), main="log2")
plot(assay(vst_cnts), main="VST")
plot(assay(rld_cnts), main="rlog")


# plot meanSD 
cntp <- meanSdPlot(counts(cnts))$gg + ggtitle("non_normalized counts")+ scale_y_continuous(limits = c(0, 3))
ntdp <- meanSdPlot(assay(ntd_cnts))$gg + ggtitle("log2")+ scale_y_continuous(limits = c(0, 3))
vstp <- meanSdPlot(assay(vst_cnts))$gg + ggtitle("VST")+ scale_y_continuous(limits = c(0, 3)) 
rldp <- meanSdPlot(assay(rld_cnts))$gg + ggtitle("rlog")+ scale_y_continuous(limits = c(0, 3))  
ggarrange(cntp, ntdp, vstp, rldp, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
# =========================================================================
# differential expressed analysis (one factor within ko, compared to control)
dds <- DESeq(dds)
colData(dds)
resultsNames(dds)

# To count-in the influence of senescence induce, we perform a group-interaction analysis

# modify design
paste0(dds$se, dds$ko)
dds$group <- factor(paste0(dds$se, dds$ko))
design(dds) <- ~group

# differential expression analysis
dds <- DESeq(dds)
colData(dds)
resultsNames(dds)

# Results
res_ko5ko3 <- results(dds, 
                      name = "group_inducedko_shXPO7_5_vs_inducedko_shXPO7_3")
res_nkoko3 <- results(dds, 
                      name = "group_inducedno_ko_vs_inducedko_shXPO7_3")
res_ctrko3 <- results(dds, 
                      name = "group_non_inducedcontrol_vs_inducedko_shXPO7_3")

res_ctrnko <- results(dds, 
                      contrast = c("group", "inducedno_ko", "non_inducedcontrol"))
res_ctrko5 <- results(dds, 
                      contrast = c("group", "inducedko_shXPO7_5", "non_inducedcontrol"))
res_nkoko5 <- results(dds, 
                      contrast = c("group", "inducedko_shXPO7_5", "inducedno_ko"))

summary(res_ko5ko3, na.rm=T); head(res_ko5ko3)
summary(res_nkoko3, na.rm=T); head(res_nkoko3)
summary(res_ctrko3, na.rm=T); head(res_ctrko3)
summary(res_ctrnko, na.rm=T); head(res_ctrnko)
summary(res_ctrko5, na.rm=T); head(res_ctrko5)
summary(res_nkoko5, na.rm=T); head(res_nkoko5)


# select significant genes based on criteria
resSig_ko5ko3_up <- subset(res_ko5ko3, padj<.01 & log2FoldChange > 1)
resSig_ko5ko3_down <- subset(res_ko5ko3, padj<.01 & log2FoldChange < -1)

resSig_nkoko3_up <- subset(res_nkoko3, padj<.01 & log2FoldChange > 1)
resSig_nkoko3_down <- subset(res_nkoko3, padj<.01 & log2FoldChange < -1)

resSig_ctrko3_up <- subset(res_ctrko3, padj<.01 & log2FoldChange > 1)
resSig_ctrko3_down <- subset(res_ctrko3, padj<.01 & log2FoldChange < -1)

resSig_ctrnko_up <- subset(res_ctrnko, padj<.01 & log2FoldChange > 1)
resSig_ctrnko_down <- subset(res_ctrnko, padj<.01 & log2FoldChange < -1)

resSig_ctrko5_up <- subset(res_ctrko5, padj<.01 & log2FoldChange > 1)
resSig_ctrko5_down <- subset(res_ctrko5, padj<.01 & log2FoldChange < -1)

resSig_nkoko5_up <- subset(res_nkoko5, padj<.01 & log2FoldChange > 1)
resSig_nkoko5_down <- subset(res_nkoko5, padj<.01 & log2FoldChange < -1)


# significant gene IDs
geneSig <- c( rownames(resSig_ko5ko3_up) ,
              rownames(resSig_ko5ko3_down) ,
              rownames(resSig_nkoko3_up) ,
              rownames(resSig_nkoko3_down) ,
              rownames(resSig_ctrko3_up) ,
              rownames(resSig_ctrko3_down) ,
              rownames(resSig_ctrnko_up) ,
              rownames(resSig_ctrnko_down) ,
              rownames(resSig_ctrko5_up) ,
              rownames(resSig_ctrko5_down) ,
              rownames(resSig_nkoko5_up) ,
              rownames(resSig_nkoko5_down) )
length(geneSig)

# Significant genes statistic data
 
log2foldchangeSig <- c(
  resSig_ko5ko3_up$log2FoldChange ,
  resSig_ko5ko3_down$log2FoldChange ,
  resSig_nkoko3_up$log2FoldChange ,
  resSig_nkoko3_down$log2FoldChange ,
  resSig_ctrko3_up$log2FoldChange ,
  resSig_ctrko3_down$log2FoldChange ,
  resSig_ctrnko_up$log2FoldChange ,
  resSig_ctrnko_down$log2FoldChange ,
  resSig_ctrko5_up$log2FoldChange ,
  resSig_ctrko5_down$log2FoldChange ,
  resSig_nkoko5_up$log2FoldChange ,
  resSig_nkoko5_down$log2FoldChange
)


padjSig <- c(
  resSig_ko5ko3_up$padj ,
  resSig_ko5ko3_down$padj ,
  resSig_nkoko3_up$padj ,
  resSig_nkoko3_down$padj ,
  resSig_ctrko3_up$padj ,
  resSig_ctrko3_down$padj ,
  resSig_ctrnko_up$padj ,
  resSig_ctrnko_down$padj ,
  resSig_ctrko5_up$padj ,
  resSig_ctrko5_down$padj ,
  resSig_nkoko5_up$padj ,
  resSig_nkoko5_down$padj
)

geneSig_sta <- as.data.frame(cbind(geneSig,log2foldchangeSig, padjSig))
geneSig_sta$cgrp <- c(
  rep("ko5ko3", length(rownames(resSig_ko5ko3_up))) ,
  rep("ko5ko3", length(rownames(resSig_ko5ko3_down))) ,
  rep("nkoko3", length(rownames(resSig_nkoko3_up))) ,
  rep("nkoko3", length(rownames(resSig_nkoko3_down))) ,
  rep("ctrko3", length(rownames(resSig_ctrko3_up))) ,
  rep("ctrko3", length(rownames(resSig_ctrko3_down))) ,
  rep("ctrnko", length(rownames(resSig_ctrnko_up))) ,
  rep("ctrnko", length(rownames(resSig_ctrnko_down))) ,
  rep("ctrko5", length(rownames(resSig_ctrko5_up))) ,
  rep("ctrko5", length(rownames(resSig_ctrko5_down))) ,
  rep("nkoko5", length(rownames(resSig_nkoko5_up))) ,
  rep("nkoko5", length(rownames(resSig_nkoko5_down)))
)

geneSig_sta$updown <- c(
  rep("up", length(rownames(resSig_ko5ko3_up))) ,
  rep("down", length(rownames(resSig_ko5ko3_down))) ,
  rep("up", length(rownames(resSig_nkoko3_up))) ,
  rep("down", length(rownames(resSig_nkoko3_down))) ,
  rep("up", length(rownames(resSig_ctrko3_up))) ,
  rep("down", length(rownames(resSig_ctrko3_down))) ,
  rep("up", length(rownames(resSig_ctrnko_up))) ,
  rep("down", length(rownames(resSig_ctrnko_down))) ,
  rep("up", length(rownames(resSig_ctrko5_up))) ,
  rep("down", length(rownames(resSig_ctrko5_down))) ,
  rep("up", length(rownames(resSig_nkoko5_up))) ,
  rep("down", length(rownames(resSig_nkoko5_down)))
)

head(geneSig_sta)


# check DE results by plotting counts
#
#ctrnko
top_ctrnko <- rownames(res_ctrnko)[which.min(res_ctrnko$padj)]
d1 <- plotCounts(dds, gene = top_ctrnko, 
                 intgroup = c("se", "ko"),
                 returnData = TRUE)

#ctrko3
top_ctrko3 <- rownames(res_ctrko3)[which.min(res_ctrko3$padj)]
d2 <- plotCounts(dds, gene = top_ctrko3, 
                 intgroup = c("se", "ko"),
                 returnData = TRUE)

#ctrko5
top_ctrko5 <- rownames(res_ctrko5)[which.min(res_ctrko5$padj)]
d3 <- plotCounts(dds, gene = top_ctrko5, 
                 intgroup = c("se", "ko"),
                 returnData = TRUE)


# plotting
ggplot(d1, aes(x=ko, y=count, color=se), x.text.angle = 90 )+
  scale_y_log10()+
  geom_beeswarm(cex=3)+
  labs(title = "Top DE Gene - control vs. no_ko") -> tp1
ggplot(d2, aes(x=ko, y=count, color=se), x.text.angle = 90)+
  scale_y_log10()+
  geom_beeswarm(cex=3)+
  labs(title = "Top DE Gene - control vs. ko3") -> tp2
ggplot(d3, aes(x=ko, y=count, color=se), x.text.angle = 90)+
  scale_y_log10()+
  geom_beeswarm(cex=3)+
  labs(title = "Top DE Gene - control vs. ko5") -> tp3

ggarrange(tp1 + rremove("x.text"),
          tp2 + rremove("x.text"),
          tp3 + rremove("x.text"),
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)




# plotMA
resultsNames(dds)

res_ko5ko3_lfcs <- lfcShrink(dds, coef = "group_inducedko_shXPO7_5_vs_inducedko_shXPO7_3", type = "apeglm")

res_nkoko3_lfcs <- lfcShrink(dds, coef = "group_inducedno_ko_vs_inducedko_shXPO7_3", type = "apeglm")

res_ctrko3_lfcs <- lfcShrink(dds, coef = "group_non_inducedcontrol_vs_inducedko_shXPO7_3", type = "apeglm")

par(mfrow=c(1,3))
plotMA(res_ko5ko3_lfcs, ylim=c(-5, 5))+title("ko5 vs.ko3")
plotMA(res_nkoko3_lfcs, ylim=c(-5, 5))+title("no_ko vs.ko3")
plotMA(res_ctrko3_lfcs, ylim=c(-5, 5))+title("control vs.ko3")


# check DE results by volcano plot

# volcano plot
# obtain gene symbols
 symb_ko5ko3 <- gene_anno[which(rownames(gene_anno)==rownames(res_ko5ko3_lfcs)),1]

symb_nkoko3 <- gene_anno[which(rownames(gene_anno)==rownames(res_nkoko3_lfcs)),1]

symb_ctrko3 <- gene_anno[which(rownames(gene_anno)==rownames(res_ctrko3_lfcs)),1]

# volcano plotting

EnhancedVolcano(res_ko5ko3_lfcs,
                # lab = rownames(res_ko5ko3_lfcs),
                lab = symb_ko5ko3,
                x = 'log2FoldChange',
                y = 'pvalue', title = "ko5 vs. ko3: Volcano Plot") -> vp1

EnhancedVolcano(res_nkoko3_lfcs,
                # lab = rownames(res_nkoko3_lfcs),
                lab = symb_nkoko3,
                x = 'log2FoldChange',
                y = 'pvalue', title = "no_ko vs. ko3: Volcano Plot") -> vp2

EnhancedVolcano(res_ctrko3_lfcs,
                # lab = rownames(res_nkoko3_lfcs),
                lab = symb_ctrko3,
                x = 'log2FoldChange',
                y = 'pvalue', title = "control vs. ko3: Volcano Plot") -> vp3

grid.arrange(vp1, vp2, vp3, nrow = 1)


# ======================================================================
# clustering of significant genes 
# using rlog normalized counts

rld_cnts_Sig <- subset(assay(rld_cnts),
                       rownames(assay(rld_cnts)) %in% geneSig)
dim(rld_cnts_Sig)


# define calculate distance function
calc_dist <- function(df){
  tempcor <- cor(df, method = "pearson")
  tempdit <- as.dist(1-tempcor/2)
  return(tempdit)
}

# distance for samples
cnts_dist <- calc_dist(rld_cnts_Sig)

# clustering for samples
cnts_clust <- hclust(cnts_dist, method = "ave")

par(mfrow=c(1,2))
cnts_clust$labels <- coldata$se
plot(cnts_clust, cex = .9)

cnts_clust$labels <- coldata$ko
plot(cnts_clust, cex = .9)


# distance for genes
cnts_g_dist <- calc_dist(t(rld_cnts_Sig))

# perform clustering
cnts_g_clust <- hclust(cnts_g_dist, method = "ave")


# define function to compute the number of clustering using silhouetteWide values
calic_k4clust <- function(in_clust, in_dist){
  mean_sil <- NULL
  mean_sil[1] <- 0
  for (i in 2:20){
    t_clusters <- cutree(in_clust, k=i)
    t_clusters_sil <- silhouette(t_clusters, dist=in_dist)
    mean_sil[i] <- mean(t_clusters_sil[,"sil_width"])
  }
  plot(mean_sil)
}

# use silhouette width values to guide the number of clustering
calic_k4clust(cnts_g_clust, cnts_g_dist)
abline(v=2, col="red")


# validation of the model for k=2
k <- 2 
cnts_g_cutg <- cutree(as.dendrogram(cnts_g_clust), k = k)
sil <- silhouette(cnts_g_cutg, dist = cnts_g_dist)
rownames(sil) <- row.names(rld_cnts_Sig)

cnts_g_cutgr <- data.frame(cnts_g_cutg, row.names(rld_cnts_Sig)) 
names(cnts_g_cutgr)  <- c("cluster", "ENSEMBL")
rownames(cnts_g_cutgr) <- cnts_g_cutgr$ENSEMBL
head(cnts_g_cutgr)

pdf()
plot(sil, main = "Silhouette plot", 
     cex.names = 0.8, col = 2:(k + 1), nmax = 100)
dev.off()

# dendrogram with k-mean clustering
plot(as.dendrogram(cnts_g_clust), leaflab = "none",
     main ="Dendrogram for significant differential expressed genes")
rect.hclust(cnts_g_clust, k = 2, border = "red")

# gene number for each cluster
table(cnts_g_cutg)


# visualize clusters using heatmap

clustgrp <- data.frame(cluster = ifelse(test = cnts_g_cutg == 1, yes = "cluster1", no = "cluster2"))

pheatmap(rld_cnts_Sig, 
         clustering_distance_cols = cnts_dist,
         clustering_distance_rows = cnts_g_dist, 
         annotation_col = coldata[ ,c(3:4)],
         annotation_row = clustgrp,
         cutree_rows = 2,
         show_rownames = F)

# heatmap for cluster
# cluster1
# distance for cluster1
cnts_clust1_dist <- calc_dist(rld_cnts_Sig[cnts_g_cutg==1,])
cnts_clust1_g_dist <- calc_dist(t(rld_cnts_Sig[cnts_g_cutg==1,]))

pheatmap(rld_cnts_Sig[cnts_g_cutg==1,], 
         clustering_distance_cols = cnts_clust1_dist,
         clustering_distance_rows = cnts_clust1_g_dist, 
         annotation_col = coldata[ ,c(3:4)],
         show_rownames = F)

# cluster2
# distance for cluster2
cnts_clust2_dist <- calc_dist(rld_cnts_Sig[cnts_g_cutg==2,])
cnts_clust2_g_dist <- calc_dist(t(rld_cnts_Sig[cnts_g_cutg==2,]))

pheatmap(rld_cnts_Sig[cnts_g_cutg==2,], 
         clustering_distance_cols = cnts_clust2_dist,
         clustering_distance_rows = cnts_clust2_g_dist, 
         annotation_col = coldata[ ,c(3:4)],
         show_rownames = F)

# perform PCA

# for samples
dataProcomp <- prcomp(t(rld_cnts_Sig), center = T, scale = T)
summary(dataProcomp)

pcadata <- cbind(as.data.frame(dataProcomp$x), coldata)
ggplot(pcadata, aes(PC1, PC2, col=se, shape=ko))+
  geom_point()+
  ggtitle("PCA plot for samples") -> pca_sp

# for genes
dataProcomp_gene <- prcomp(rld_cnts_Sig, center = T, scale = T)
summary(dataProcomp_gene)

# anno_Sig <- gene_anno[rownames(rld_cnts_Sig),1]
pcadata_gene <- as.data.frame(dataProcomp_gene$x)

ggplot(pcadata_gene, aes(PC1, PC2))+
  geom_point()+
  ggtitle("PCA plot for genes") -> pca_gp

grid.arrange(pca_sp,pca_gp, nrow = 1)


# cluster1
dataProcomp_clust1_gene <- prcomp(rld_cnts_Sig[cnts_g_cutg==1,], center = T, scale = T)

pcadata_clust1_gene <- as.data.frame(dataProcomp_clust1_gene$x)

ggplot(pcadata_clust1_gene, aes(PC1, PC2))+
  geom_point()+
  ggtitle("PCA plot for cluster1 genes") 

# cluster2
dataProcomp_clust2_gene <- prcomp(rld_cnts_Sig[cnts_g_cutg==2,], center = T, scale = T)

pcadata_clust2_gene <- as.data.frame(dataProcomp_clust2_gene$x)

ggplot(pcadata_clust2_gene, aes(PC1, PC2))+
  geom_point()+
  ggtitle("PCA plot for cluster2 genes") 

##############################################################################
# Gene Enrichment Analysis

head(gene_anno) 
head(cnts_g_cutgr)
head(geneSig_sta)
head(rld_cnts_Sig)


# cluster gene names
geneSig_clust1 <- cnts_g_cutgr[cnts_g_cutgr$cluster==1, 2]
geneSig_clust2 <- cnts_g_cutgr[cnts_g_cutgr$cluster==2, 2]
length(unique(geneSig))==length(geneSig_clust1)+length(geneSig_clust2)

# up and down gene names
geneSig_up <- unique(geneSig_sta[geneSig_sta$updown=="up", 1])
geneSig_down <- unique(geneSig_sta[geneSig_sta$updown=="down", 1])


# ===
# significant genes enrichment analysis
# all significant genes
geneSig_param <- new("GOHyperGParams",
                       geneIds = gene_anno[geneSig,3],
                       universeGeneIds=gene_anno$ENTREZID,
                       annotation="org.Hs.eg.db",
                       ontology="BP",
                       pvalueCutoff=0.001,
                       testDirection="over")

geneSig_overRepresented <- hyperGTest(geneSig_param)
head(summary(geneSig_overRepresented),10)

# clust1
geneSig_clust1_param <- new("GOHyperGParams",
                     geneIds = gene_anno[geneSig_clust1,3],
                     universeGeneIds=gene_anno$ENTREZID,
                     annotation="org.Hs.eg.db",
                     ontology="BP",
                     pvalueCutoff=0.001,
                     testDirection="over")

geneSig_clust1_overRepresented <- hyperGTest(geneSig_clust1_param)
head(summary(geneSig_clust1_overRepresented),10)

# clust2
geneSig_clust2_param <- new("GOHyperGParams",
                            geneIds = gene_anno[geneSig_clust2,3],
                            universeGeneIds=gene_anno$ENTREZID,
                            annotation="org.Hs.eg.db",
                            ontology="BP",
                            pvalueCutoff=0.001,
                            testDirection="over")

geneSig_clust2_overRepresented <- hyperGTest(geneSig_clust2_param)
head(summary(geneSig_clust2_overRepresented),10)

# ===
# Geneset enrichment analysis

# get a gene list(entrez id) of decreasing sorted foldchanges
# for significant genes
geneSig_foldchag_sort <- geneSig_sta[order(as.numeric(geneSig_sta$log2foldchangeSig), decreasing = T), c(1:2)]
geneSig_foldchag_sort$entrezid <- gene_anno[geneSig_foldchag_sort$geneSig,3]


geneSig_foldchang_list <- as.numeric(geneSig_foldchag_sort$log2foldchangeSig)
names(geneSig_foldchang_list) <- geneSig_foldchag_sort$entrezid

geneSig_gse <- gseGO(geneList = geneSig_foldchang_list,
                    ont = "ALL",
                    keyType = "ENTREZID", 
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = .05,
                    verbose = TRUE,
                    OrgDb = org.Hs.eg.db,
                    pAdjustMethod = "fdr")

# output
require(DOSE)
dotplot(geneSig_gse, showCategory=10, split=".sign") + 
  facet_grid(.~.sign)
# dotplot(geneSig_gse)


# pathway enrichment distribution with p-value
ridgeplot(geneSig_gse, label_format = 20, showCategory = 10)+
  labs(x="Enrichment distribution")


# GSEA plot
gseaplot(geneSig_gse, by="all", 
         title=geneSig_gse$Description[1], geneSetID = 1)


# === cluster 1
# get a gene list(entrez id) of decreasing sorted foldchanges
# for cluster1 genes
geneSig_clust1_sta <- geneSig_sta[geneSig_sta$geneSig %in% geneSig_clust1,]

geneSig_clust1_foldchag_sort <- 
  geneSig_clust1_sta[order(as.numeric(geneSig_clust1_sta$log2foldchangeSig), decreasing = T), c(1:2)]

geneSig_clust1_foldchag_sort$entrezid <- gene_anno[geneSig_clust1_foldchag_sort$geneSig,3]

geneSig_clust1_foldchang_list <- as.numeric(geneSig_clust1_foldchag_sort$log2foldchangeSig)

names(geneSig_clust1_foldchang_list) <- geneSig_clust1_foldchag_sort$entrezid  

geneSig_clust1_gse <- gseGO(geneList = geneSig_clust1_foldchang_list,
                     ont = "ALL",
                     keyType = "ENTREZID", 
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = .05,
                     verbose = TRUE,
                     OrgDb = org.Hs.eg.db,
                     pAdjustMethod = "fdr")

# output
require(DOSE)
dotplot(geneSig_clust1_gse, showCategory=10, split=".sign") + 
  facet_grid(.~.sign)


# pathway enrichment distribution with p-value
ridgeplot(geneSig_clust1_gse, label_format = 20, showCategory = 10)+
  labs(x="Enrichment distribution")


# GSEA plot
gseaplot(geneSig_clust1_gse, by="all", 
         title=geneSig_clust1_gse$Description[1], geneSetID = 1)
 
# 
#  === cluster 2
# 
geneSig_clust2_sta <- geneSig_sta[geneSig_sta$geneSig %in% geneSig_clust2,]

geneSig_clust2_foldchag_sort <- 
  geneSig_clust2_sta[order(as.numeric(geneSig_clust2_sta$log2foldchangeSig), decreasing = T), c(1:2)]

geneSig_clust2_foldchag_sort$entrezid <- gene_anno[geneSig_clust2_foldchag_sort$geneSig,3]

geneSig_clust2_foldchang_list <- as.numeric(geneSig_clust2_foldchag_sort$log2foldchangeSig)

names(geneSig_clust2_foldchang_list) <- geneSig_clust2_foldchag_sort$entrezid  

geneSig_clust2_gse <- gseGO(geneList = geneSig_clust2_foldchang_list,
                            ont = "ALL",
                            keyType = "ENTREZID", 
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = .05,
                            verbose = TRUE,
                            OrgDb = org.Hs.eg.db,
                            pAdjustMethod = "fdr")

# output
# require(DOSE)
dotplot(geneSig_clust2_gse, showCategory=10, split=".sign") + 
  facet_grid(.~.sign)


# pathway enrichment distribution with p-value
ridgeplot(geneSig_clust2_gse, label_format = 20, showCategory = 10)+
  labs(x="Enrichment distribution")


# GSEA plot
gseaplot(geneSig_clust2_gse, by="all", 
         title=geneSig_clust2_gse$Description[1], geneSetID = 1)



# ============= Over-Representation Analysis ===================================
# for significant genes

# obtain ENTREZID list for significant genes
# geneSig_entz <- gene_anno[unique(geneSig),3]

geneSig_enrichgo <- enrichGO(gene = gene_anno[unique(geneSig),3],
                              universe = gene_anno$ENTREZID,
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

## use simplify to remove redundant terms
geneSig_enrichgo <- simplify(geneSig_enrichgo, 
                             cutoff=0.3, by="p.adjust",
                             select_fun=min)

# === cluster1

geneSig_clust1_enrichgo <- enrichGO(gene = gene_anno[geneSig_clust1,3],
                             universe = gene_anno$ENTREZID,
                             OrgDb = org.Hs.eg.db, 
                             keyType = "ENTREZID",
                             readable = T,
                             ont = "BP",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.10)

## use simplify to remove redundant terms
geneSig_clust1_enrichgo <- simplify(geneSig_clust1_enrichgo, 
                                    cutoff=0.3, by="p.adjust",
                                    select_fun=min)
# === cluster2

geneSig_clust2_enrichgo <- enrichGO(gene = gene_anno[geneSig_clust2,3],
                                    universe = gene_anno$ENTREZID,
                                    OrgDb = org.Hs.eg.db, 
                                    keyType = "ENTREZID",
                                    readable = T,
                                    ont = "BP",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

## use simplify to remove redundant terms
geneSig_clust2_enrichgo <- simplify(geneSig_clust2_enrichgo, 
                                    cutoff=0.3, by="p.adjust",
                                    select_fun=min)


## plotting
upsetplot(geneSig_enrichgo)
upsetplot(geneSig_clust1_enrichgo)
upsetplot(geneSig_clust2_enrichgo)


barplot(geneSig_enrichgo, drop=T, showCategory = 10, 
        title = "GoBiologicalPathwaysfor significant genes", 
        font.size = 8)

barplot(geneSig_clust1_enrichgo, drop=T, showCategory = 10, 
        title = "GoBiologicalPathwaysfor cluster 1", 
        font.size = 8)
barplot(geneSig_clust2_enrichgo, drop=T, showCategory = 10, 
        title = "GoBiologicalPathwaysfor cluster 2", 
        font.size = 8)


dotplot(geneSig_enrichgo)
dotplot(geneSig_clust1_enrichgo)
dotplot(geneSig_clust2_enrichgo)


goplot(geneSig_enrichgo,showCategory = 10)
goplot(geneSig_clust1_enrichgo,showCategory = 10)
goplot(geneSig_clust2_enrichgo, showCategory = 10)

# top 10 genes for geneSig, cluster1, and cluster2
print("Top fold change genes in significant gene list")
head((gene_anno[unique(geneSig_foldchag_sort$geneSig), c(1:2)]),10)
print("================")
print("Top fold change genes in cluster 1 gene list")
head((gene_anno[unique(geneSig_clust1_foldchag_sort$geneSig), c(1:2)]),10)
print("================")
print("Top fold change genes in cluster 2 gene list")
head((gene_anno[unique(geneSig_clust2_foldchag_sort$geneSig), c(1:2)]),10)
session_info()
