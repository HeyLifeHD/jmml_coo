#Libraries
library(DSS)
library(bsseq)
library(doParallel)
library(ChIPseeker)
library(foreach)
library(SummarizedExperiment)
library(doMC)
library(rtracklayer)
library(HDF5Array)
library(org.Mm.eg.db)
library(ggpubr)
library(randomcoloR)
library(RColorBrewer)
library(VennDiagram)
library(ChIPpeakAnno)
library(dendextend)
library(pheatmap)
library("ape")
library(dplyr)

#load my own data
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200803_DMR_array_HSC_comb"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno_reduced.rds"))
#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
dmrs_final$all <- dmrs_red
mcols(dmrs_final$all)$direction <- "hypo"



#load anands data
comb_norm_dat = readRDS(file = "/home/heyj/c010-projects/Anand/HSC_tree_datasets_for_Jsochka/combined_arrays_GRset.Rds")
source("/home/heyj/c010-projects/Anand/HSC_tree_datasets_for_Jsochka/accessory_scripts.R")

#prepare color and data
col_pal = c(ggsci::pal_d3(palette = "category10", alpha = 1)(10), "gray70", "darkred")
comb_norm_dat$pch = ifelse(test = comb_norm_dat$GSE == "GSE110554", yes = 21, no = 22)
names(col_pal) = unique(comb_norm_dat$Cell_type)

pdf(file.path(analysis.dir, "PCA_healthyHierachy_array.pdf"), width=10, height=10)
temp <- pcaPlot(dat = comb_norm_dat, sampGroups = comb_norm_dat$Cell_type, pal = col_pal, legendPos = "bottomright", legendNCol = 2, numPositions = 25000, pch = comb_norm_dat$pch)
dev.off()

#cell-type Annotations
samp_anno = as.data.frame(pData(object = comb_norm_dat))
samp_anno$Cell_Type_colors = col_pal[samp_anno$Cell_type]
#trajectory analysis
set.seed(seed = 1024)
phyl_dist = stats::dist(t(getBeta(object = comb_norm_dat)[1:25000,]), method = "man")
dist_nj = ape::nj(X = phyl_dist)
pdf(file = file.path(analysis.dir, "Trajectory_healthyHierachy_array.pdf"), width = 3, height = 3, bg = "white", pointsize = 9)
par(mar = c(0, 0, 0, 0))
plot_tree(tree = dist_nj, col_df = samp_anno[,c("Cell_Type_colors"), drop = FALSE], col_id = 1, type = 'unrooted', rotate = 270, cex = 0.01)
legend(x = "topright", legend = names(col_pal), col = col_pal, ncol = 3, cex = 1.2, pch = 19)
dev.off()

#Same trees with DMPs (each cell type v/s HSc)
#Results from comparing each cell type against HSC
comb_norm_dat_dmp = readRDS(file = "/home/heyj/c010-projects/Anand/HSC_tree_datasets_for_Jsochka/DMP_list_everything_vs_HSC.Rds")
col_pal = c(ggsci::pal_d3(palette = "category10", alpha = 1)(10), "gray70", "darkred")
names(col_pal) = unique(comb_norm_dat_dmp$Cell_type)
pdf(file = file.path(analysis.dir, "PCA_healthyHierachyDMPs_array.pdf"), width = 4, height = 4, bg = "white", pointsize = 9)
temp = pcaPlot(dat = comb_norm_dat_dmp, sampGroups = comb_norm_dat_dmp$Cell_type, pal = col_pal, legendPos = "bottomright", legendNCol = 2, numPositions = nrow(comb_norm_dat_dmp), pt_cex = 1.5)
dev.off()

set.seed(seed = 1024)
phyl_dist = stats::dist(t(getBeta(object = comb_norm_dat_dmp)), method = "man")
dist_nj = ape::nj(X = phyl_dist)

pdf(file = file.path(analysis.dir, "Trajectory_healthyHierachyDMP_array.pdf"), width = 3, height = 3, bg = "white", pointsize = 9)
par(mar = c(0, 0, 0, 0))
plot_tree(tree = dist_nj, col_df = samp_anno[,c("Cell_Type_colors"), drop = FALSE], col_id = 1, type = 'unrooted', rotate = 270, cex = 0.01)
legend(x = "bottomright", legend = names(col_pal), col = col_pal, ncol = 2, cex = 1, pch = 19)
dev.off()



#extract cpg annotation for comb_norm_dat_dmp
cpg_oi <- rownames(getBeta(object = comb_norm_dat_dmp))
beta_oi <- getBeta(object = comb_norm_dat_dmp)
anno <- granges(comb_norm_dat_dmp)
anno_oi <- anno[cpg_oi,]
meth_oi<- bsseq::getMeth(bsseq_all, regions= anno_oi, type = "raw", what=c("perRegion"))
rownames(meth_oi)<- cpg_oi
#subset tumor
pheno <- pData(bsseq_all)
pheno_tumor <- pheno[pheno$Sample_Type =="tumor",]
meth_oi_tumor <- meth_oi[,rownames(pheno_tumor)]
#get complete cases
idx <- complete.cases(meth_oi_tumor)
table(idx)
meth_oi_tumor_sub <- meth_oi_tumor[idx,]
beta_oi_sub <- beta_oi[idx,]

#get color anntation normal
col_pal_normal = c(ggsci::pal_d3(palette = "category10", alpha = 1)(10), "gray70", "darkred")
names(col_pal_normal) = unique(comb_norm_dat_dmp$Cell_type)
comb_norm_dat_dmp$Cell_type
colors_final_normal <- samp_anno[,c("Cell_Type_colors"), drop = FALSE]
colnames(colors_final_normal)<- "color"
#get color annotation tumor
pheno_tumor$group <- as.character(pheno_tumor$Epigenotype)
colors <- data.frame(group=unique(pheno_tumor$group), color=c("#0041a5", "#b5241c"))
pheno_color <- left_join(as.data.frame(pheno_tumor), colors)
rownames(pheno_color)<- pheno_color$Name
colors_final <- pheno_color[,"color", drop=FALSE]
col_pal_tumor <- as.character(colors$color)
names(col_pal_tumor) <- colors$group

#combine data
comb_data <- cbind(meth_oi_tumor_sub, beta_oi_sub)
groups <- c( pheno_tumor$group , comb_norm_dat_dmp$Cell_type)
col_pal_comb <- c(col_pal_tumor,col_pal_normal)
colors_final_comb <- rbind(colors_final,colors_final_normal )

#plot pca
pdf(file = file.path(analysis.dir, "PCA_healthyHierachyDMPs_array_jmmlIncl_Epigenotype.pdf"), width = 4, height = 4, bg = "white", pointsize = 9)
temp = pcaPlot(dat = comb_data, sampGroups = groups, pal = col_pal_comb, legendPos = "bottomright", legendNCol = 2, numPositions = nrow(comb_norm_dat_dmp), pt_cex = 1.5)
dev.off()

#tree
phyl_dist = stats::dist(t(comb_data), method = "man")
dist_nj = ape::nj(X = phyl_dist)


pdf(file = file.path(analysis.dir, "PhylTree_healthyHierachyDMPs_array_jmmlIncl_Epigenotype.pdf"), width = 7, height = 7, bg = "white", pointsize = 9)
plot_tree(tree = dist_nj, col_df = colors_final_comb, col_id = 1, type = 'unrooted', rotate = 80, cex = 0.01)
legend(x = "bottomleft", legend = names(col_pal_comb), col = col_pal_comb, ncol = 2, cex = 1.2, pch = 19)
dev.off()