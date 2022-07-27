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
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

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
files <- list.files(path="/omics/groups/OE0219/internal/JMMLC/scRNA_Data/New/DEG/JMML-HCA_CB_Sampled-Patient", pattern=".rds",recursive=TRUE, full.names=TRUE)
temp <- strsplit(files, "-", ,fixed=TRUE)
temp  <- sapply(temp, function(x)x[length(x)])
temp <- sapply(strsplit(temp, ".", ,fixed=TRUE), "[", 1)
cell <- sapply(strsplit(files, "-", ,fixed=TRUE), "[", 5)

names(files)<- paste0(temp, "_vs_","cordblood_HSC", ".",cell)
expr <- lapply(files, function(x){
    x <-readRDS(x)
    x$SYMBOL <-x$gene 
    x
    }
)
expr_sep <- list()
for(i in unique(cell)){
    expr_sep[[i]]<- list()
    expr_sep[[i]]<- expr[grep(i, names(files))]
    names(expr_sep[[i]]) <- sapply(strsplit(names(expr_sep[[i]]), ".", fixed=TRUE),"[",1)
}


#merge data
library(dplyr)
merged <- list()
for(i in names(dmrs_final)){
    merged[[i]]<- list()
    for(cell in names(expr_sep)){
        merged[[i]][[cell]]<- list()
        dir.create(file.path(analysis.dir, i, "expr_integration", "allDMRs", cell), recursive=TRUE)
            for(j in names(expr_sep[[cell]])){
                merged[[i]][[cell]][[j]] <-  inner_join(as.data.frame(dmrs_final[[i]]), expr_sep[[cell]][[j]], by="SYMBOL")
                print(paste(i, j, dim(merged[[i]][[cell]][[j]])))
                #get annotation of sign
                merged[[i]][[cell]][[j]]$Col <- NA
                merged[[i]][[cell]][[j]]$Col <- ifelse(merged[[i]][[cell]][[j]]$p_val_adj < 0.05 & merged[[i]][[cell]][[j]]$avg_logFC >1 , "sign. upregulated",ifelse(merged[[i]][[cell]][[j]]$p_val_adj < 0.05 & merged[[i]][[cell]][[j]]$avg_logFC<(-1), "sign. downregulated", "no sign. change"))
                merged[[i]][[cell]][[j]]$Col <- ifelse(is.na(merged[[i]][[cell]][[j]]$Col), "no sign. change", merged[[i]][[cell]][[j]]$Col)
                merged[[i]][[cell]][[j]]$Col<- as.factor(merged[[i]][[cell]][[j]]$Col)
                merged[[i]][[cell]][[j]]$Plot<-NA
                merged[[i]][[cell]][[j]][merged[[i]][[cell]][[j]]$Col!="no sign. change",]$Plot <- merged[[i]][[cell]][[j]][merged[[i]][[cell]][[j]]$Col!="no sign. change",]$SYMBOL
                write.table( merged[[i]][[cell]][[j]], file.path(analysis.dir, i, "expr_integration", "allDMRs", cell, paste0("expression_methylation_JMML-HCA_CB_Sampled-Patient_",j,".txt")), quote=FALSE, sep="\t", row.names=FALSE)
            }
    }
}

#plot
#all DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
names(col)<- c("no sign. change",  "sign. downregulated", "sign. upregulated" )

g <- list()
for(i in names(merged)){
    for(cell in names(merged[[i]])){
            dir.create(file.path(analysis.dir, "scatter_plots", "allDMRs", i, cell ), recursive=TRUE)
    for(j in names(merged[[i]][[cell]])){
        name <- gsub("_", " ", j)
        g <- ggscatter(merged[[i]][[cell]][[j]], x = "avg_logFC", y = j,xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf( file.path(analysis.dir, i, "expr_integration", "allDMRs", cell, paste0("scatter_expression_methylation_JMML-HCA_CB_Sampled-Patient_",j,".pdf")), width=4,height=4)
            print(g)
        dev.off()
        g <- ggscatter(merged[[i]][[cell]][[j]], x = "avg_logFC", y = j,xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
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
        pdf( file.path(analysis.dir, i, "expr_integration", "allDMRs", cell, paste0("scatter_expression_methylation_JMML-HCA_CB_Sampled-Patient_",j,"_label.pdf")), width=4,height=4)
            print(g)
        dev.off()
        print(j)
    }
    print(i)
    }
}



#Promoter DMRs
#plot
#all DMRs
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged)){
    for(cell in names(merged[[i]])){
            dir.create( file.path(analysis.dir, i, "expr_integration", "PromoterDMRs", cell), recursive=TRUE)
    for(j in names(merged[[i]][[cell]])){
        name <- gsub("_", " ", j)
        merged_sub <- merged[[i]][[cell]][[j]]
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
        pdf( file.path(analysis.dir, i, "expr_integration", "PromoterDMRs", cell, paste0("scatter_expression_methylation_JMML-HCA_CB_Sampled-Patient_",j,".pdf")), width=4,height=4)
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
        pdf( file.path(analysis.dir, i, "expr_integration", "PromoterDMRs", cell, paste0("scatter_expression_methylation_JMML-HCA_CB_Sampled-Patient_",j,"_label.pdf")), width=4,height=4)
            print(g)
        dev.off()
        print(j)
    }
    print(i)
    }
}


