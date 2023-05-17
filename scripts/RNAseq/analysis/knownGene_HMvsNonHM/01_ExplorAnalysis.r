#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
#function
plotPCA.alt <- function (object, intgroup = "condition", ntop = Inf, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + 
    geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}
#folder
base.dir<- "/omics/groups/OE0219/internal/jmmlc_rnaseq/220805_rnaseq_knownDEG"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#DESeq2 Analysis
dds <-readRDS(file.path(results.dir, "dds.rds"))

#Extracting transformed values
vst <- vst(dds)
anno <- colData(vst)
saveRDS(vst,file = file.path(results.dir,"vst.rds"))
#vst <- readRDS(file.path(results.dir,"vst.rds"))

#anno
pbat_col = list(
  Tissue = c(tumor = "#99a637", prenatal = "#252525", cordblood = "#737373", adult_bonemarrow = "#ababab"), 
  Celltype = c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"), 
  Genotype = c(neg = "#529a51", KRAS = "#99a637", PTPN11 = "#007458", wildtype ="#ababab"),
  Patient = c(D117 = "#0058b4", D129 = "#2188c9", 
            D217 = "#fbbb25", I217 = "#fca349", 
            D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
            normal ="#ababab"), 
  Epigenotype = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4", wildtype ="#ababab")
  )
col <- pbat_col$Epigenotype
anno_colors <- list(consensusCluster3=col)
anno <- colData(vst)
annovst <- anno[,c("consensusCluster3" ), drop=FALSE]

#plot pca
pcaDatavst<- plotPCA(vst, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarvst <- round(100* attr(pcaDatavst, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst, x="PC1", y="PC2",
            size=5,
            color = "consensusCluster3", 
            label= "EWOG_ID",
            repel=TRUE,
            legend = "right",palette=col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst[2], "% variance")) )
dev.off()

for(i in colnames(anno)[-c(23,25)]){
    pdf(file.path(PreDE.dir, paste0("PCA12_", name,".pdf")), height = 5, width = 6)
    print(ggscatter(pcaDatavst, x="PC1", y="PC2",
                size=5,
                color = i, 
                label= "EWOG_ID",
                repel=TRUE,
                legend = "right",
                ellipse = F ,mean.point = FALSE,
                star.plot = F,  xlab=(paste0("PC1: ", percentvarvst[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst[2], "% variance")) ))
    dev.off()
}


#pca of 5000 mv genes
pcaDatavst<- plotPCA(vst, intgroup =  colnames(anno), returnData =TRUE, ntop=5000)
percentvarvst <- round(100* attr(pcaDatavst, "percentVar"))
pdf(file.path(PreDE.dir, "PCA12_5000mvGenes.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst, x="PC1", y="PC2",
            size=5,
            color = "consensusCluster3", 
            label= "EWOG_ID",
            repel=TRUE,
            legend = "right",palette=col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst[2], "% variance")) )
dev.off()

#function to plot PC2 and 3
pcaDatavst.alt<- plotPCA.alt(vst, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
percentvarvst.alt <- round(100* attr(pcaDatavst.alt, "percentVar"))
pdf(file.path(PreDE.dir, "PCA23.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst.alt, x="PC2", y="PC3",            
          size=5,
            color = "consensusCluster3", 
            label= "EWOG_ID",
            repel=TRUE,
            legend = "right",palette=col,
            ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarvst.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarvst.alt[2], "% variance")) )
dev.off()

#pca 2 and 3 of 5000 mv genes
pcaDatavst.alt<- plotPCA.alt(vst, intgroup = colnames(anno), returnData =TRUE, ntop=5000)
percentvarvst.alt <- round(100* attr(pcaDatavst.alt, "percentVar"))
pdf(file.path(PreDE.dir, "PCA23_5000mvGenes.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst.alt, x="PC2", y="PC3",            
          size=5,
            color = "consensusCluster3", 
            label= "EWOG_ID",
            repel=TRUE,
            legend = "right",palette=col,
            ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarvst.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarvst.alt[2], "% variance")) )
dev.off()

#Sample Clustering correlation
mat <- assay(vst)
dissimilarity <- 1 - cor(mat, use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
col1 <- hrld$labels
names(col1)<- as.character(anno[col1,]$consensusCluster3)
col1 <- col1[order.dendrogram(dend)]
col1 <- col[names(col1)]
labels(dend)<- as.character(anno[labels(dend),]$consensusCluster3)
dend <- dend %>% 
set("branches_lwd", 2) %>%
set("labels_colors",col1) %>% 
set("labels_cex", .6 )%>%
set("leaves_pch", 19)%>% 
set("leaves_cex", 1.5)%>% 
set("leaves_col", col1)
pdf(file.path(PreDE.dir, "Clustering_correlation.pdf"), height = 5, width = 5)
dend %>% plot
dev.off()

#Gene clustering: Heatmaps
#500 most variable repeats
topVarrepeatsvst<- head(order(rowVars(assay(vst)), decreasing=TRUE),  1000)
matvst <- assay(vst)[topVarrepeatsvst,]
pdf(file.path(PreDE.dir,"Heatmap1000vst_Scale.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors)
dev.off()

topVarrepeatsvst<- head(order(rowVars(assay(vst)), decreasing=TRUE),  100)
matvst <- assay(vst)[topVarrepeatsvst,]
pdf(file.path(PreDE.dir,"Heatmap100vst_Scale.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

#look if metadata is sign. associated with pcs
ref_meth = assay(vst)
sampleAnnotation = as.data.frame(sample_anno)
sampleAnnotation <- sampleAnnotation[, -c(1, 2,3,25,28,29)]
scale = T; center = T
pr = prcomp(ref_meth, scale. = scale, center = center)
vars = 100 * pr$sdev^2/sum(pr$sdev^2)
res = list()
for (i in 4:ncol(sampleAnnotation)) {
  print(i)
  if (is.factor(sampleAnnotation[, i]) | is.character(sampleAnnotation[, i]) | is.logical(sampleAnnotation[, i])) {
    annot = as.factor(sampleAnnotation[, i])
    res[[names(sampleAnnotation)[i]]] = apply(
      pr$rotation, 2, function(x) {
        if (length(unique(annot)) > 2) {
          kruskal.test(split(x, annot))$p.value
        }else {
          wilcox.test(split(x, annot)[[1]], split(x, annot)[[2]])$p.value
        }
  })
  }
  else {
    annot = sampleAnnotation[, i]
    res[[names(sampleAnnotation)[i]]] = apply(
      pr$rotation, 2, function(x) {
        cor.test(x, annot)$p.value
  })
  }
}
pc_res = list(pcs = pr$rotation, vars = vars, association = do.call("rbind", res))
pc_res$association[,1:3]< 0.05
sapply(pc_res[,c(1:3)], function(x)x <0.05)