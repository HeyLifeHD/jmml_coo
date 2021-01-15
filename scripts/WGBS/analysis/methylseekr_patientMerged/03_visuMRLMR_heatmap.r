#Libraries
library(MethylSeekR)
library(parallel)
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPseeker)
library(ggpubr)
library(pheatmap)
#load my own data
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
methylseekr.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr"
dir.create(analysis.dir)

#load data
bsseq_merged <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_PatientMerged_CelltypeMerged.rds"))
UMRLMRsegments.gr <- readRDS(file.path(methylseekr.dir,  paste0("UMRsLMRs", ".rds")))

#subset references of interest for heatmap
samples <- names(UMRLMRsegments.gr)
samples_sub <- samples[1:9]
UMRLMRsegments.gr_sub <- UMRLMRsegments.gr[samples_sub]
#look at umr
UMRsegments.gr <- lapply(UMRLMRsegments.gr_sub , function(x){
    x <- x[x$type=="UMR",]
    x
})
#reduce ranges
red_UMR <- reduce(unlist(GRangesList(UMRsegments.gr )))
length(red_UMR )
#get methylation values
meth_UMR <- as.data.frame(bsseq::getMeth(bsseq_merged, regions= red_UMR, type = "raw", what=c("perRegion")))
meth_UMR_sub <- meth_UMR[complete.cases(meth_UMR),1:9]
dim(meth_UMR_sub)
#plot heatmap
pbat_col = list(
    Patient = c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000", normal ="#e9e9e9"), 
    Genotype = c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458", wildtype ="#e9e9e9"),
    Epigenotype = c(HM = "#b5241c", LM = "#0041a5", control_CB ="#e9e9e9", control_adult ="#c6c6c6"), 
    Tumor = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d",
                HSC_adult ="#252525", MPP_adult = "#737373", CMP_adult = "#ffb86f", GMP_adult = "#e27e37", MEP_adult="#d2624a",
                HSC_CB ="#252525"), 
    Celltype = c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#c6a27f", MEP = "#d2624a", CMP = "#ffb86f", GMP = "#e27e37")
    )
annovst <- as.data.frame(colData(bsseq_merged))[, c("Patient", "Genotype", "Epigenotype")] 
    
pheatmap(meth_UMR_sub,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
    clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
    scale="row",  annotation_color=pbat_col,
    filename=file.path(methylseekr.dir,"Heatmap_UMRs_rowScale_WarD2_Manhatten_sub.pdf"))
pheatmap(meth_UMR_sub,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
    clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
    scale="none",fontsize_row=5,  annotation_color=pbat_col,
    filename=file.path(methylseekr.dir,"Heatmap_UMRs_noScale_WarD2_Manhatten_sub.pdf"))



#look at lmr
LMRsegments.gr <- lapply(UMRLMRsegments.gr_sub , function(x){
    x <- x[x$type=="LMR",]
    x
})
#reduce ranges
red_LMR <- reduce(unlist(GRangesList(LMRsegments.gr )))
length(red_LMR )
#get methylation values
meth_LMR <- as.data.frame(bsseq::getMeth(bsseq_merged, regions= red_LMR, type = "raw", what=c("perRegion")))
meth_LMR_sub <- meth_UMR[complete.cases(meth_LMR),1:9]
dim(meth_LMR_sub)
#plot heatmap
pbat_col = list(
    Patient = c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000", normal ="#e9e9e9"), 
    Genotype = c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458", wildtype ="#e9e9e9"),
    Epigenotype = c(HM = "#b5241c", LM = "#0041a5", control_CB ="#e9e9e9", control_adult ="#c6c6c6"), 
    Tumor = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d",
                HSC_adult ="#252525", MPP_adult = "#737373", CMP_adult = "#ffb86f", GMP_adult = "#e27e37", MEP_adult="#d2624a",
                HSC_CB ="#252525"), 
    Celltype = c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#c6a27f", MEP = "#d2624a", CMP = "#ffb86f", GMP = "#e27e37")
    )
annovst <- as.data.frame(colData(bsseq_merged))[, c("Patient", "Genotype", "Epigenotype")] 
    
pheatmap(meth_LMR_sub,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
    clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
    scale="row",  annotation_color=pbat_col,
    filename=file.path(methylseekr.dir,"Heatmap_LMRs_rowScale_WarD2_Manhatten_sub.pdf"))
pheatmap(meth_LMR_sub,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
    clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
    scale="none",fontsize_row=5,  annotation_color=pbat_col,
    filename=file.path(methylseekr.dir,"Heatmap_LMRs_noScale_WarD2_Manhatten_sub.pdf"))
