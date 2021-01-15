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

#load my own data: Epigenotype vs HSC_cb
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
input_dmr.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"
analysis.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200826_expressionIntegration_changes"
dir.create(analysis.dir)

#load data
#bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(input_dmr.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
dmrs_red<- readRDS(file.path(input_dmr.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))


#load expression data
files <- list.files(path="/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/expression/200825_DEG/HM_vs_LM/", pattern=".csv", full.names=TRUE)
expr <- as.data.frame(data.table::fread(files))
expr$SYMBOL <- expr$gene

merged <-  dplyr::inner_join(as.data.frame(dmrs_final$EpigenotypeHM), expr, by="SYMBOL")
#get annotation of sign
merged$Col <- NA
merged$Col <- ifelse(merged$p_val_adj < 0.05 & merged$avg_logFC >0.25 & merged$pct.1 >0.1 & merged$pct.2 >0.1, "sign. upregulated",ifelse(merged$p_val_adj < 0.05 & merged$avg_logFC<(-0.25)& merged$pct.1 >0.1 & merged$pct.2 >0.1, "sign. downregulated", "no sign. change"))
merged$Col <- ifelse(is.na(merged$Col), "no sign. change", merged$Col)
merged$Col<- as.factor(merged$Col)
merged$Plot<-NA
merged[merged$Col!="no sign. change",]$Plot <- merged[merged$Col!="no sign. change",]$SYMBOL
dir.create(file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM"), recursive=TRUE)
write.table( merged, file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM" , paste0("expression_methylation",".txt")), quote=FALSE, sep="\t", row.names=FALSE)


#plot
#all DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
dir.create(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", "EpigenotypeHM"), recursive=TRUE)

name <-"HM_vs_LM"
g <- ggscatter(merged, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM" , paste0("expression_methylation.pdf")), width=4,height=4)
    print(g)
dev.off()
g <- ggscatter(merged, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    font.label = c(3, "plain","black"),
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM"  , paste0("expression_methylation_label.pdf")), width=7,height=7)
    print(g)
dev.off()


#plot
#PromoterDMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
merged_sub <- merged[grep("Promoter", merged$annotation),]
name <-"HM_vs_LM"
g <- ggscatter(merged_sub, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", "EpigenotypeHM" , paste0("expression_methylation.pdf")), width=4,height=4)
    print(g)
dev.off()
g <- ggscatter(merged_sub, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    font.label = c(3, "plain","black"),
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", "EpigenotypeHM"  , paste0("expression_methylation_label.pdf")), width=7,height=7)
    print(g)
dev.off()




#all data

#load expression data
files <- list.files(path="/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/expression/200825_DEG/HM_vs_LM/", pattern=".rds", full.names=TRUE)
expr <- as.data.frame(readRDS(files))
expr$SYMBOL <- expr$gene

merged <-  dplyr::inner_join(as.data.frame(dmrs_final$EpigenotypeHM), expr, by="SYMBOL")
print(paste(i, j, dim(merged)))
#get annotation of sign
merged$Col <- NA
merged$Col <- ifelse(merged$p_val_adj < 0.05 & merged$avg_logFC >0.25 & merged$pct.1 >0.1 & merged$pct.2 >0.1, "sign. upregulated",ifelse(merged$p_val_adj < 0.05 & merged$avg_logFC<(-0.25)& merged$pct.1 >0.1 & merged$pct.2 >0.1, "sign. downregulated", "no sign. change"))
merged$Col <- ifelse(is.na(merged$Col), "no sign. change", merged$Col)
merged$Col<- as.factor(merged$Col)
merged$Plot<-NA
merged[merged$Col!="no sign. change",]$Plot <- merged[merged$Col!="no sign. change",]$SYMBOL
dir.create(file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM"), recursive=TRUE)
write.table( merged, file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM" , paste0("expression_methylation",".txt")), quote=FALSE, sep="\t", row.names=FALSE)


#plot
#all DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
dir.create(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", "EpigenotypeHM"), recursive=TRUE)

name <-"HM_vs_LM"
g <- ggscatter(merged, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM" , paste0("expression_methylation_allGenes.pdf")), width=4,height=4)
    print(g)
dev.off()
g <- ggscatter(merged, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    font.label = c(3, "plain","black"),
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "allDMRs", "EpigenotypeHM"  , paste0("expression_methylation_label_allGenes.pdf")), width=7,height=7)
    print(g)
dev.off()


#plot
#PromoterDMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
merged_sub <- merged[grep("Promoter", merged$annotation),]
name <-"HM_vs_LM"
g <- ggscatter(merged_sub, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", "EpigenotypeHM" , paste0("expression_methylation.pdf")), width=4,height=4)
    print(g)
dev.off()
g <- ggscatter(merged_sub, x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Methyladion difference","\n",name),
    color = "Col" ,
    shape = 16, size = 1.5,
    palette=col, # Points color, shape and size 61
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    font.label = c(3, "plain","black"),
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
    )+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")
pdf(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", "EpigenotypeHM"  , paste0("expression_methylation_label.pdf")), width=7,height=7)
    print(g)
dev.off()

paste0(unique(merged[merged$avg_logFC>0 & merged$diff.Methy < 0,]$SYMBOL), collapse=", ")
paste0(unique(merged[merged$avg_logFC<0 & merged$diff.Methy > 0,]$SYMBOL), collapse=", ")