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

#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb"

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
#dmrs_final$all <- dmrs_red
#mcols(dmrs_final$all)$direction <- "hypo"

#prepare matrix for overlap comparison
comp <- dmrs_red  
mcols(comp)<- NULL
ol <- list()
for(i in 1:length(dmrs_final)){
    ol[[i]] <- comp %over% dmrs_final[[i]]

}
ol <- do.call("cbind", ol)
colnames(ol)<- names(dmrs_final)
ol <- as.data.frame(ol)

ol[ol==TRUE]<- 1
dir.create(file.path(analysis.dir, "comparison"))
pdf(file.path(analysis.dir, "comparison", "comparison_DMR_regions_allSamples_allRegions.pdf"), width=14, height=10)
upset(ol,nsets= 20,nintersects=100,order.by = "freq")
dev.off()

pdf(file.path(analysis.dir, "comparison", "comparison_DMR_regions_allSamples_top10Int.pdf"), width=14, height=10)
upset(ol,nsets= 20,nintersects=10,order.by = "freq")
dev.off()

