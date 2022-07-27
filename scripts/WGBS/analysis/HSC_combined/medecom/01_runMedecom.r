#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb"


#Librarieslibrary(data.table)
library(bsseq)
library(MeDeCom)
library(pheatmap)
library(RColorBrewer)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

#select dmrs
dmrs <- dmrs_red

#get dmr methylation values
meth_all<-bsseq::getMeth(bsseq_all, regions=dmrs, what="perRegion", type="raw")
dim(meth_all)
meth_all <- meth_all[complete.cases(meth_all),]
dim(meth_all)
dir.create(file.path(analysis.dir  , "medecom_allSamples"))
saveRDS(meth_all,file.path(analysis.dir  , "medecom_allSamples","meth_medecom_input.rds"))

#run medecom
meth_all <- readRDS(file.path(analysis.dir  , "medecom_allSamples","meth_medecom_input.rds"))
medecom.result<-runMeDeCom(as.matrix(meth_all), 5:20, c(0,10^(-5:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=8)
saveRDS(medecom.result, file.path(analysis.dir  , "medecom_allSamples", "allDMRred_Medecom_results.rds"))

#plot CV for each K
pdf(file.path(analysis.dir, "medecom_allSamples", "allDMRs_Prameters.pdf"))
plotParameters(medecom.result)
dev.off()
#plot RSME for each K
pdf(file.path(analysis.dir, "medecom_allSamples", "allDMRs_Prameters_RMSE.pdf"))
plotParameters(medecom.result, statistic="RMSE")
dev.off()
#Plot CV, RSME and objective for each K
for (i in 5:20){
    dir.create(file.path(analysis.dir,"medecom_allSamples",paste0("K_",as.character(i))))
    pdf(file.path(analysis.dir, "medecom_allSamples",paste0("K_",i),  "all_DMRred_allParameters.pdf"))
    plotParameters(medecom.result,  K=i, lambdaScale="log")
    dev.off()
}
#get pheno 
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
    annovst <- as.data.frame(colData(bsseq_all))[, c("Patient", "Genotype", "Epigenotype", "Tumor", "Celltype")] 

col <- RColorBrewer::brewer.pal(12, "Reds")
#Plot Proportions for each K
for (i in 5:20){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.001)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom_allSamples",paste0("K_",i),  "all_DMRred_Proportions_l0001.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

for (i in 5:20){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.01)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom_allSamples",paste0("K_",i),  "all_DMRred_Proportions_l001.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

for (i in 5:20){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.1)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom_allSamples",paste0("K_",i),  "all_DMRred_Proportions_l01.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

#save proportions of interest
#medecom.result <- readRDS(file.path(analysis.dir, "medecom_allSamples", "allDMRred_Medecom_results.rds"))
#proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.01)
#saveRDS(proportions,file.path(analysis.dir, "medecom_allSamples","proportions_k5_001.rds"))


#same for 100 000 random regions
#get dmr methylation values
meth_all<-bsseq::getMeth(bsseq_all,  type="raw")
dim(meth_all)
meth_all <- meth_all[complete.cases(meth_all),]
dim(meth_all)
#select 100 000 random regions
meth_all <- meth_all[sample( nrow(meth_all),100000),]
dim(meth_all)

dir.create(file.path(analysis.dir  , "medecom_allSamples_randomCpGs"))
saveRDS(meth_all,file.path(analysis.dir  , "medecom_allSamples_randomCpGs","meth_medecom_input.rds"))

#run medecom
meth_all <- readRDS(file.path(analysis.dir  , "medecom_allSamples_randomCpGs","meth_medecom_input.rds"))
medecom.result<-runMeDeCom(as.matrix(meth_all), 5:20, c(0,10^(-5:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=8)
saveRDS(medecom.result, file.path(analysis.dir  , "medecom_allSamples_randomCpGs", "allDMRred_Medecom_results.rds"))

#plot CV for each K
pdf(file.path(analysis.dir, "medecom_allSamples_randomCpGs", "allDMRs_Prameters.pdf"))
plotParameters(medecom.result)
dev.off()
#plot RSME for each K
pdf(file.path(analysis.dir, "medecom_allSamples_randomCpGs", "allDMRs_Prameters_RMSE.pdf"))
plotParameters(medecom.result, statistic="RMSE")
dev.off()
#Plot CV, RSME and objective for each K
for (i in 5:20){
    dir.create(file.path(analysis.dir,"medecom_allSamples_randomCpGs",paste0("K_",as.character(i))))
    pdf(file.path(analysis.dir, "medecom_allSamples_randomCpGs",paste0("K_",i),  "all_DMRred_allParameters.pdf"))
    plotParameters(medecom.result,  K=i, lambdaScale="log")
    dev.off()
}
#get pheno 
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
    annovst <- as.data.frame(colData(bsseq_all))[, c("Patient", "Genotype", "Epigenotype", "Tumor", "Celltype")] 

col <- RColorBrewer::brewer.pal(12, "Reds")
#Plot Proportions for each K
for (i in 5:20){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.001)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom_allSamples_randomCpGs",paste0("K_",i),  "all_DMRred_Proportions_l0001.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

for (i in 5:20){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.01)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom_allSamples_randomCpGs",paste0("K_",i),  "all_DMRred_Proportions_l001.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

for (i in 5:20){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.1)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom_allSamples_randomCpGs",paste0("K_",i),  "all_DMRred_Proportions_l01.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}
