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
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
input_dmr.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200805_DMR_Model_CelltypeGroup_sub_cbHSC"
analysis.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200826_expressionIntegration_changes"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(input_dmr.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
dmrs_red<- readRDS(file.path(input_dmr.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))
dmrs_final$all <- dmrs_red

#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

#create further subset
dmrs_final_sub <- lapply(dmrs_final, function(x){
    x <- x[x %outside% dmrs_HSC_red,]
    x
})
lapply(dmrs_final, length)
lapply(dmrs_final_sub, length)
names(dmrs_final_sub)<- paste0(names(dmrs_final_sub), "_hierachy_DMRs_removed")
dmrs_final <- c(dmrs_final, dmrs_final_sub)

#annotate dmrs
dmrs_anno<- list()
dmrs_anno_df<- list()
#plot files and annotate
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
for(i in names(dmrs_final)){
#annotate dmrs
dmrs_anno[[i]] <- annotatePeak(peak= dmrs_final[[i]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

#extract annotated granges
dmrs_anno_df[[i]] <- as.data.frame(dmrs_anno[[i]])
mcols(dmrs_final[[i]]) <- dmrs_anno_df[[i]]
}

#adjust bsseq object
pheno <- pData(bsseq_all)
pheno$group <-as.character(pheno$Patient)
pheno[pheno$Sample_Type=="normal", ]$group <-"cordblood_HSC"
pheno$group <- as.factor(pheno$group )
pheno$group <- droplevels(pheno$group)
pheno$group <- relevel(pheno$group , "cordblood_HSC")
pheno$Celltype <- droplevels(pheno$Celltype)
pheno$Celltype <- relevel(pheno$Celltype , "HSC")
pData(bsseq_all)<- pheno

#define contrasts to calculate methylation difference for every patient - cord blood hsc
contrasts <- t(expand.grid(as.character(unique(pheno[pheno$Sample_Type=="tumor", ]$Patient)), "cordblood_HSC"))
colnames(contrasts)<- paste0(contrasts[1,], "_vs_", contrasts[2,])

#calculate methylation differences for hm as well as lm dmrs
meth_diff <- list()
for(i in names(dmrs_final)){
    meth_diff[[i]]<- list()
    for(j in colnames(contrasts)){
         meth_diff[[i]][[j]]<- rowMeans( as.matrix(mcols(dmrs_final[[i]][, rownames(pheno[pheno$group==contrasts[1,j],])])),na.rm=TRUE) - rowMeans( as.matrix(mcols(dmrs_final[[i]][, rownames(pheno[pheno$group==contrasts[2,j],])])),na.rm=TRUE)
    }
   meth_diff[[i]] <-  do.call("cbind",meth_diff[[i]])
   mcols(dmrs_final[[i]]) <- cbind(mcols(dmrs_final[[i]]),  meth_diff[[i]])
}

#load expression data
files <- list.files(path="/omics/groups/OE0219/internal/jmmlc_pbat/data/expression/200825_DEG/DEG_HCA_CB_HSC_HC/", pattern=".rds", full.names=TRUE)
temp <- sapply(strsplit(files, "-", ,fixed=TRUE), "[", 6)
temp <- sapply(strsplit(temp, ".", ,fixed=TRUE), "[", 1)
names(files)<- paste0(temp, "_vs_","cordblood_HSC")
expr <- lapply(files, function(x){
    x <-readRDS(x)
    x$SYMBOL <-x$gene 
    x
    }
)

#merge data
library(dplyr)
merged <- list()
for(i in names(dmrs_final)){
    merged[[i]]<- list()
    for(j in colnames(contrasts)){
        merged[[i]][[j]] <-  inner_join(as.data.frame(dmrs_final[[i]]), expr[[j]], by="SYMBOL")
        print(paste(i, j, dim(merged[[i]][[j]])))
        #get annotation of sign
        merged[[i]][[j]]$Col <- NA
        merged[[i]][[j]]$Col <- ifelse(merged[[i]][[j]]$p_val_adj < 0.05 & merged[[i]][[j]]$avg_logFC >1 & merged[[i]][[j]]$pct.1 >0.1 & merged[[i]][[j]]$pct.2 >0.1, "sign. upregulated",ifelse(merged[[i]][[j]]$p_val_adj < 0.05 & merged[[i]][[j]]$avg_logFC<(-1)& merged[[i]][[j]]$pct.1 >0.1 & merged[[i]][[j]]$pct.2 >0.1, "sign. downregulated", "no sign. change"))
        merged[[i]][[j]]$Col <- ifelse(is.na(merged[[i]][[j]]$Col), "no sign. change", merged[[i]][[j]]$Col)
        merged[[i]][[j]]$Col<- as.factor(merged[[i]][[j]]$Col)
        merged[[i]][[j]]$Plot<-NA
        merged[[i]][[j]][merged[[i]][[j]]$Col!="no sign. change",]$Plot <- merged[[i]][[j]][merged[[i]][[j]]$Col!="no sign. change",]$SYMBOL
        write.table( merged[[i]][[j]], file.path(analysis.dir, "scatter_plots", "allDMRs", i , paste0("expression_methylation_",j,".txt")), quote=FALSE, sep="\t", row.names=FALSE)
    }
}

#plot
#all DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged)){
    dir.create(file.path(analysis.dir, "scatter_plots", "allDMRs", i ), recursive=TRUE)
    for(j in names(merged[[i]])){
        name <- gsub("_", " ", j)
        g <- ggscatter(merged[[i]][[j]], x = "avg_logFC", y = j,xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf(file.path(analysis.dir, "scatter_plots", "allDMRs", i , paste0("expression_methylation_",j,".pdf")), width=4,height=4)
            print(g)
        dev.off()
        g <- ggscatter(merged[[i]][[j]], x = "avg_logFC", y = j,xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf(file.path(analysis.dir, "scatter_plots", "allDMRs", i , paste0("expression_methylation_",j,"_label.pdf")), width=14,height=14)
            print(g)
        dev.off()
        print(j)
    }
    print(i)
}

#Promoter DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged)){
    dir.create(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", i ), recursive=TRUE)
    for(j in names(merged[[i]])){
        name <- gsub("_", " ", j)
        merged_sub <- merged[[i]][[j]]
        merged_sub <- merged_sub[grep("Promoter", merged_sub$annotation),]
        g <- ggscatter(merged_sub, x = "avg_logFC", y = j,xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", i , paste0("expression_methylation_",j,".pdf")), width=4,height=4)
            print(g)
        dev.off()
        g <- ggscatter(merged_sub, x = "avg_logFC", y = j,xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf(file.path(analysis.dir, "scatter_plots", "PromoterDMRs", i , paste0("expression_methylation_",j,"_label.pdf")), width=14,height=14)
            print(g)
        dev.off()
        print(j)
    }
    print(i)
}


