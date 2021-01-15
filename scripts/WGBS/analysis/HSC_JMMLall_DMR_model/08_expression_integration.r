
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
library(LOLA)

#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

#create further subset
dmrs_final_sub <- lapply(dmrs_final, function(x){
    x <- x[x %outside% dmrs_HSC_red,]
    x
})
dmrs_final_sub2 <- lapply(dmrs_final, function(x){
    x <- x[x %over% dmrs_HSC_red,]
    x
})
names(dmrs_final_sub)<- paste0(names(dmrs_final_sub), "_hierachy_DMRs_removed")
names(dmrs_final_sub2)<- paste0(names(dmrs_final), "_only_hierachy_DMRs")

lapply(dmrs_final, length)
lapply(dmrs_final_sub, length)
lapply(dmrs_final_sub2, length)
dmrs_final <- c(dmrs_final, dmrs_final_sub,dmrs_final_sub2)

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
#dmrs_final$all <- dmrs_red
#mcols(dmrs_final$all)$direction <- "hypo"

#load expression data
expr <- readRDS("/icgc/dkfzlsdf/analysis/C010/JMMLC/scRNA_Data/New/DEG/JMML-HCA_CB_Sampled-panJMML/JMML-HCA_CB_Sampled-CD34-HSC.DEG.rds")
expr$SYMBOL <- expr$gene

#merge data
library(dplyr)
merged <- list()
for(i in names(dmrs_final)){
        merged[[i]]<-  inner_join(as.data.frame(dmrs_final[[i]]), expr, by="SYMBOL")
        print(paste(dim(merged[[i]])))
        #get annotation of sign
        merged[[i]]$Col <- NA
        merged[[i]]$Col <- ifelse(merged[[i]]$p_val_adj < 0.05 & merged[[i]]$avg_logFC >1 , "sign. upregulated",ifelse(merged[[i]]$p_val_adj < 0.05 & merged[[i]]$avg_logFC<(-1), "sign. downregulated", "no sign. change"))
        merged[[i]]$Col <- ifelse(is.na(merged[[i]]$Col), "no sign. change", merged[[i]]$Col)
        merged[[i]]$Col<- as.factor(merged[[i]]$Col)
        merged[[i]]$Plot<-NA
        merged[[i]][merged[[i]]$Col!="no sign. change",]$Plot <- merged[[i]][merged[[i]]$Col!="no sign. change",]$SYMBOL
        dir.create( file.path(analysis.dir, i, "expr_integration", "allDMRs"), recursive=TRUE)
        write.table( merged[[i]], file.path(analysis.dir, i, "expr_integration", "allDMRs",  paste0("expression_methylation_pan_JMML-HCA_CB_Sampled-CD34_",i,".txt")), quote=FALSE, sep="\t", row.names=FALSE)
}


#plot
#all DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
names(col)<- c("no sign. change",  "sign. downregulated", "sign. upregulated" )
g <- list()
for(i in names(merged)){
    dir.create(file.path(analysis.dir, "scatter_plots", "allDMRs", i ), recursive=TRUE)
        name <- gsub("_", " ", i)
        g <- ggscatter(merged[[i]], x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf( file.path(analysis.dir, i, "expr_integration", "allDMRs",  paste0("scatter_expression_methylation_pan_JMML-HCA_CB_Sampled-CD34_",i,".pdf")), width=4,height=4)
            print(g)
        dev.off()
        g <- ggscatter(merged[[i]], x = "avg_logFC", y = "diff.Methy",xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf(file.path(analysis.dir, i, "expr_integration", "allDMRs",  paste0("scatter_expression_methylation_pan_JMML-HCA_CB_Sampled-CD34_",i,"_labeled.pdf")), width=4,height=4)
            print(g)
        dev.off()
    }


#Promoter DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
names(col)<- c("no sign. change",  "sign. downregulated", "sign. upregulated" )
g <- list()
for(i in names(merged)){
    dir.create(file.path(analysis.dir, i, "expr_integration", "PromoterDMRs" ), recursive=TRUE)
        name <- gsub("_", " ", i)
        merged_sub <- merged[[i]]
        merged_sub <- merged_sub[grep("Promoter", merged_sub$annotation),]
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
        pdf( file.path(analysis.dir, i, "expr_integration", "PromoterDMRs",  paste0("scatter_expression_methylation_pan_JMML-HCA_CB_Sampled-CD34_",i,".pdf")), width=4,height=4)
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
        pdf( file.path(analysis.dir, i, "expr_integration", "PromoterDMRs",  paste0("scatter_expression_methylation_pan_JMML-HCA_CB_Sampled-CD34_",i,"_label.pdf")), width=4,height=4)
            print(g)
        dev.off()
        print(i)
    }


