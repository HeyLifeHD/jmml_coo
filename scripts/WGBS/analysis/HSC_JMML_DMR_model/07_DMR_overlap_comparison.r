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
library(UpSetR)
#load my own data: Epigenotype vs HSC_cb
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

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
pdf(file.path(analysis.dir, "comparison", "comparison_DMR_regions_mainComp_allRegions.pdf"), width=7, height=7)
upset(ol,nsets= 20,nintersects=100,order.by = "freq")
dev.off()

pdf(file.path(analysis.dir, "comparison", "comparison_DMR_regions_mainComp_top10Int.pdf"), width=7, height=7)
upset(ol,nsets= 20,nintersects=10,order.by = "freq")
dev.off()


#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

##create further subset
#dmrs_final_sub <- lapply(dmrs_final, function(x){
#    x <- x[x %outside% dmrs_HSC_red,]
#    x
#})
#lapply(dmrs_final, length)
#lapply(dmrs_final_sub, length)
#names(dmrs_final_sub)<- paste0(names(dmrs_final_sub), "_hierachy_DMRs_removed")
dmrs_final$hierachyDMRs <- dmrs_HSC_red
dmrs_final<- lapply(dmrs_final, function(x){
    mcols(x)<- NULL
    x
})
#prepare matrix for overlap comparison
comp <- reduce(unlist(GRangesList(dmrs_final)))
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
pdf(file.path(analysis.dir, "comparison", "comparison_DMR_regions_subComp_allRegions.pdf"), width=7, height=7)
upset(ol,nsets= 20,nintersects=100,order.by = "freq")
dev.off()

pdf(file.path(analysis.dir, "comparison", "comparison_DMR_regions_subComp_top10Int.pdf"), width=7, height=7)
upset(ol,nsets= 20,nintersects=10,order.by = "freq")
dev.off()

