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
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load methylation data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#extract promoter methylation
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ENS<-EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
genes <- genes(ENS)
prom <- promoters(genes, upstream=500, downstream=500)
prom <- prom[seqnames(prom) %in%  c(paste0("chr", 1:21),"chrY", "chrX")]
meth_prom <- bsseq::getMeth(bsseq_all, regions= prom, type = "raw", what=c("perRegion"))
mcols(prom)<- cbind(mcols(prom), meth_prom)
#saveRDS(prom,file.path(input.dir ,"bsseq", "PromoterMeth_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
#prom <- readRDS(file.path(input.dir ,"bsseq", "PromoterMeth_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))

#load expression data
expr <- readRDS("/icgc/dkfzlsdf/analysis/C010/JMMLC/scRNA_Data/New/DEG/Statistics/HSC.PseudoBulk.rds")
#select expressed promoters
expr_sub <- expr[rownames(expr) %in%  prom$gene_id,]
prom_sub <- as.data.frame(mcols(prom))
prom_sub<- prom_sub[rownames(expr_sub),]
#get average promoter methylation per patient and cordblood hsc samples
pheno <-pData(bsseq_all)
group <- unique(pheno$Donor)
group <- lapply(group, function(x){
x <-rownames(pheno[pheno$Donor %in% x,])
x
})
names(group)<-  unique(pheno$Donor)
prom_sub_mean <- lapply(group, function(x){
    if(length(x) >1){
    x <-  rowMeans(prom_sub[,x], na.rm=TRUE)
    }else{
    x <- prom_sub[,x, drop=TRUE]
    }
   x
})
prom_sub_mean <- do.call("cbind", prom_sub_mean)
#get sample annotation the same
expr_sub <- expr_sub[,c("D129", "D213", "I217", "D217", "D360", "D123", "D124", "D117", "HCA_CB_Sampled")]
colnames(expr_sub)<- c("D129", "D213", "I217", "D217", "D360", "D123", "D124", "D117", "cordblood")
expr_sub <- expr_sub[, colnames(prom_sub_mean)]
dim(expr_sub)
dim(prom_sub_mean)
table(colnames(prom_sub_mean)==colnames(expr_sub))
table(rownames(prom_sub_mean)==rownames(expr_sub))

#subset only complete observations
idx <- complete.cases(prom_sub_mean)
prom_sub_mean <- prom_sub_mean[idx,]
expr_sub <- expr_sub[idx,]

#remove cordblood
expr_sub <- expr_sub[, -9]
prom_sub_mean <- prom_sub_mean[, -9]

#perform correlation test
corr <- list()
    for (j in rownames(expr_sub)){
        corr[[j]]<- cor.test(as.numeric(expr_sub[j,]), as.numeric(prom_sub_mean[j,]), method="pearson", alternative = "two.sided", exact=F)
    }
#extract results
corr_res <- data.frame(p.value=unlist(lapply(corr, function(x)x$p.value)), estimate=unlist(lapply(corr, function(x)x$estimate)))
corr_res$symbol <- prom[rownames(corr_res),]$symbol
#subset meaningful results and adjust for multiple testing
corr_res_sub <- corr_res[complete.cases(corr_res),]
corr_res_sub$p.adjust <- p.adjust(corr_res_sub$p.value, "hochberg")
corr_res_sub<- corr_res_sub[order(corr_res_sub$p.adjust, decreasing=FALSE),]
head(corr_res_sub, 20)

#plot results as histogramm
corr_res_sub$sign <- ifelse(corr_res_sub$p.adjust<0.1, "adj. p-value < 0.1", "adj. p-value > 0.1")
dir.create(file.path(analysis.dir,"correlation_analysis", "promoter_meth_500bp"), recursive=TRUE)
pdf(file=file.path(analysis.dir,"correlation_analysis", "promoter_meth_500bp", "correlation_res_histogram.pdf"), height=3.5, width=5)
ggpubr::gghistogram(corr_res_sub, x="estimate", bins=30, add="mean",
    xlab="Pearson correlation estimate\nAverage DMR methylation vs closest gene expression",ylab="Number of observations",
    fill="sign", palette=c("red", "black"))
dev.off()


#plot correlation of significant genes
corr_res_sub_sig <- corr_res_sub[which(corr_res_sub$p.adjust<0.1), ]

anno_col <- c(D117 = "#0058b4", D129 = "#2188c9", 
            D217 = "#fbbb25", I217 = "#fca349", 
            D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
            cordblood ="#ababab") 

dir.create(file.path(analysis.dir,"correlation_analysis", "promoter_meth_500bp","scatter_plots"))
for(i in rownames(corr_res_sub_sig)){
    temp <- data.frame(Expression=as.numeric(expr_sub[i,]), Methylation=as.numeric(prom_sub_mean[i,]), Donor=colnames(prom_sub_mean))
    pdf(file=file.path(analysis.dir,"correlation_analysis", "promoter_meth_500bp","scatter_plots", paste0("scatter_plot_", i, ".pdf")), height=7, width=7)
        print(ggpubr::ggscatter(temp, x="Expression",y="Methylation",
            xlab="Normalized gene epression",ylab="Average promoter methylation",color="Donor", palette=anno_col,
            shape = 16, size = 2, title=paste0(i, "\n", "estimate = ", round(corr_res_sub_sig[i,]$estimate , 3),"; p.adjust = ", round(corr_res_sub_sig[i,]$p.adjust , 5) ),
            add = "reg.line",  # Add regressin line
            conf.int = T, # Add confidence interval
            add.params = list(color = "grey", fill = "lightgray"),
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")))
    dev.off()
}
