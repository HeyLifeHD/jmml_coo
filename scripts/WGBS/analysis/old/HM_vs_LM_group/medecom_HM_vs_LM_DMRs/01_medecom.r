
#Directories
odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"



#Librarieslibrary(data.table)
library(bsseq)
library(MeDeCom)
library(pheatmap)
library(RColorBrewer)

#load data
#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub_cov.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

#get dmr methylation values
meth_all<-bsseq::getMeth(bsseq_all, regions=dmrs, what="perRegion", type="raw")
meth_all <- meth_all[complete.cases(meth_all),]
dim(meth_all)
dir.create(file.path(analysis.dir  , "medecom"))
saveRDS(meth_all,file.path(analysis.dir  , "medecom","meth_medecom_input.rds"))

medecom.result<-runMeDeCom(as.matrix(meth_all), 2:15, c(0,10^(-5:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=10)
saveRDS(medecom.result, file.path(analysis.dir  , "medecom", "allDMRred_Medecom_results.rds"))

#plot CV for each K
pdf(file.path(analysis.dir, "medecom", "allDMRs_Prameters.pdf"))
plotParameters(medecom.result)
dev.off()
#plot RSME for each K
pdf(file.path(analysis.dir, "medecom", "allDMRs_Prameters_RMSE.pdf"))
plotParameters(medecom.result, statistic="RMSE")
dev.off()
#Plot CV, RSME and objective for each K
for (i in 2:15){
    dir.create(file.path(analysis.dir,"medecom",paste0("K_",as.character(i))))
    pdf(file.path(analysis.dir, "medecom",paste0("K_",i),  "all_DMRred_allParameters.pdf"))
    plotParameters(medecom.result,  K=i, lambdaScale="log")
    dev.off()
}
#get pheno 
pbat_col = list(
    Patient = c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"), 
    Genotype = c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458"),
    Epigenotype = c(HM = "#b5241c", LM = "#0041a5"), 
    tumor = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d"))
annovst <- as.data.frame(colData(bsseq_all))[, c("Epigenotype", "Patient", "tumor", "Genotype")] 

col <- RColorBrewer::brewer.pal(12, "Reds")
#Plot Proportions for each K
for (i in 2:15){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.001)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom",paste0("K_",i),  "all_DMRred_Proportions_l0001.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

for (i in 2:15){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.01)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom",paste0("K_",i),  "all_DMRred_Proportions_l001.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

for (i in 2:15){
proportions <- MeDeCom::getProportions(medecom.result, K=i, lambda=0.1)
colnames(proportions)<- colnames(meth_all)
pdf(file.path(analysis.dir, "medecom",paste0("K_",i),  "all_DMRred_Proportions_l01.pdf"), height=5)
print(pheatmap(proportions, scale="none", show_colnames=F, color=col,
annotation_col=annovst, annotation_color=pbat_col))
dev.off()
}

#save proportions of interest
medecom.result <- readRDS(file.path(analysis.dir, "medecom", "allDMRred_Medecom_results.rds"))
proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.01)
saveRDS(proportions)
