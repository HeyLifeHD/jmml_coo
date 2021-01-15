#compare group comparison with model results
dmrs_model <- readRDS("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200608_DMR_model/dmrs_gr_sub_MethDiff_anno.rds")
dmrs_model_reduced <- readRDS("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200608_DMR_model/dmrs_gr_sub_MethDiff_anno_reduced.rds")
dmrs_group <- readRDS("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM/sig_dmrs_5in4_sub_anno.rds")

#check for overlap ov everything
#create new reduced bacckground with everything
BG <- reduce(unlist(GRangesList(dmrs_model_reduced, dmrs_group)))
FG<- dmrs_model
FG$group_HM_vs_LM <- dmrs_group
length(BG)
for(i in names(FG)){
    mcols(BG)[,i] <- 0
    mcols(BG[queryHits(findOverlaps(BG, FG[[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(BG))
#do plotting
library(UpSetR)
dir.create("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200608_DMR_model/comparison_group/")
pdf("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200608_DMR_model/comparison_group/upset_all_comp.pdf", width=10, height=7)
upset(regions, order.by = "freq", nsets=10, nintersects=20)
dev.off()

#check for overlap group vs model
#create new reduced bacckground
BG <- reduce(unlist(GRangesList(dmrs_model$EpigenotypeHM, dmrs_group)))
FG<- list(dmrs_model$EpigenotypeHM,dmrs_group )
names(FG)<- c("model_HM_vs_LM", "group_HM_vs_LM")
length(BG)
for(i in names(FG)){
    mcols(BG)[,i] <- 0
    mcols(BG[queryHits(findOverlaps(BG, FG[[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(BG))
#do plotting
library(UpSetR)
pdf("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200608_DMR_model/comparison_group/upset_Epigenotype_comp.pdf", width=10, height=7)
upset(regions, order.by = "freq", nsets=10, nintersects=20)
dev.off()