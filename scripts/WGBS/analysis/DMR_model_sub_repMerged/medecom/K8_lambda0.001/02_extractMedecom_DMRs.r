##### Joschka Hey 
##### 
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
library(pheatmap)
library(randomcoloR)
library(MeDeCom)

#directories
output.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"

#load data
bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))
dmrs_red_anal <- dmrs_final$EpigenotypeHM
meth_all <- readRDS(file.path(analysis.dir  , "medecom","meth_medecom_input.rds"))

#get dmr methylation values
meth_all<-bsseq::getMeth(bsseq_all, regions=dmrs_red_anal, what="perRegion", type="raw")
dim(meth_all)
idx <- complete.cases(meth_all)
meth_all <- meth_all[idx,]
dim(meth_all)
dmrs_red_anal <- dmrs_red_anal[idx,]
length(dmrs_red_anal)
saveRDS(dmrs_red_anal,file.path(analysis.dir, "medecom","dmrs_red_anal.rds" ))
#load medecom results 
medecom.result <- readRDS(file.path(analysis.dir, "medecom", "allDMRred_Medecom_results.rds"))

#extract T matrix for k05
T <-medecom.result@outputs$`1`$T
T <- T[,"lambda_0.001"]$K_8
#extract top 0.5 difference for all K 
dmr_LMC_05 <- list()
for(i in 1:8){
    T_diff <- rowMeans(T[,-i])-T[,i]
    T_hypo <- T_diff>0.5
    T_hyper <- T_diff<(-0.5)
    T_changed <- list(dmrs_red_anal[which(T_hypo),],dmrs_red_anal[which(T_hyper),] )
    names(T_changed)<- c("hypo", "hyper")
    dmr_LMC_05[[i]]<- T_changed
}
names(dmr_LMC_05)<- c("LMC1","LMC2", "LMC3", "LMC4", "LMC5", "LMC6", "LMC7","LMC8")
lapply(dmr_LMC_05, function(x)lapply(x, function(y)length(y)))
#unique(as.factor(dmr_LMC_05$LMC1$hypo$SYMBOL))
dir.create(file.path(analysis.dir, "medecom","downstream_K8_lambda0.001"))
saveRDS(dmr_LMC_05, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_05.rds"))


#extract top 0.4 difference for all K 
dmr_LMC_04 <- list()
for(i in 1:8){
    T_diff <- rowMeans(T[,-i])-T[,i]
    T_hypo <- T_diff>0.4
    T_hyper <- T_diff<(-0.4)
    T_changed <- list(dmrs_red_anal[which(T_hypo),],dmrs_red_anal[which(T_hyper),] )
    names(T_changed)<- c("hypo", "hyper")
    dmr_LMC_04[[i]]<- T_changed
}
names(dmr_LMC_04)<- c("LMC1","LMC2", "LMC3", "LMC4", "LMC5", "LMC6", "LMC7","LMC8")
lapply(dmr_LMC_04, function(x)lapply(x, function(y)length(y)))
#unique(as.factor(dmr_LMC_08$LMC1$hypo$SYMBOL))
saveRDS(dmr_LMC_04, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_04.rds"))

#extract top 0.3 difference for all K 
dmr_LMC_03 <- list()
for(i in 1:8){
    T_diff <- rowMeans(T[,-i])-T[,i]
    T_hypo <- T_diff>0.3
    T_hyper <- T_diff<(-0.3)
    T_changed <- list(dmrs_red_anal[which(T_hypo),],dmrs_red_anal[which(T_hyper),] )
    names(T_changed)<- c("hypo", "hyper")
    dmr_LMC_03[[i]]<- T_changed
}
names(dmr_LMC_03)<- c("LMC1","LMC2", "LMC3", "LMC4", "LMC5", "LMC6", "LMC7","LMC8")
lapply(dmr_LMC_03, function(x)lapply(x, function(y)length(y)))
#unique(as.factor(dmr_LMC_08$LMC1$hypo$SYMBOL))
saveRDS(dmr_LMC_03, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_03.rds"))

#extract top 0.8 difference for all K 
dmr_LMC_08 <- list()
for(i in 1:8){
    T_diff <- rowMeans(T[,-i])-T[,i]
    T_hypo <- T_diff>0.8
    T_hyper <- T_diff<(-0.8)
    T_changed <- list(dmrs_red_anal[which(T_hypo),],dmrs_red_anal[which(T_hyper),] )
    names(T_changed)<- c("hypo", "hyper")
    dmr_LMC_08[[i]]<- T_changed
}
names(dmr_LMC_08)<- c("LMC1","LMC2", "LMC3", "LMC4", "LMC5", "LMC6", "LMC7","LMC8")
lapply(dmr_LMC_08, function(x)lapply(x, function(y)length(y)))
#unique(as.factor(dmr_LMC_08$LMC1$hypo$SYMBOL))
saveRDS(dmr_LMC_08, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_08.rds"))


#extract top 0.7 difference for all K 
dmr_LMC_07 <- list()
for(i in 1:8){
    T_diff <- rowMeans(T[,-i])-T[,i]
    T_hypo <- T_diff>0.7
    T_hyper <- T_diff<(-0.7)
    T_changed <- list(dmrs_red_anal[which(T_hypo),],dmrs_red_anal[which(T_hyper),] )
    names(T_changed)<- c("hypo", "hyper")
    dmr_LMC_07[[i]]<- T_changed
}
names(dmr_LMC_07)<- c("LMC1","LMC2", "LMC3", "LMC4", "LMC5", "LMC6", "LMC7","LMC8")
lapply(dmr_LMC_07, function(x)lapply(x, function(y)length(y)))
#unique(as.factor(dmr_LMC_07$LMC1$hypo$SYMBOL))
saveRDS(dmr_LMC_07, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_07.rds"))

#combine as a list
dmr_LMC_list <- list(dmr_LMC_03, dmr_LMC_04, dmr_LMC_05, dmr_LMC_07, dmr_LMC_08)
names(dmr_LMC_list)<- c("dmr_LMC_03", "dmr_LMC_04", "dmr_LMC_05", "dmr_LMC_07", "dmr_LMC_08")

#subset only promoters
dmr_LMC_list_promoters <- list()
for(cut in names(dmr_LMC_list)){
    dmr_LMC_list_promoters[[cut]]<-list()
    for(lmc in names(dmr_LMC_list[[cut]])){
         dmr_LMC_list_promoters[[cut]][[lmc]]<-list()
         for(dir in names(dmr_LMC_list[[cut]][[lmc]])){
                dmr_LMC_list_promoters[[cut]][[lmc]][[dir]] <- dmr_LMC_list[[cut]][[lmc]][[dir]][grep("Promoter", dmr_LMC_list[[cut]][[lmc]][[dir]]$annotation),]
                print(paste(cut, lmc, dir, length(dmr_LMC_list_promoters[[cut]][[lmc]][[dir]]), sep=": "))
        }
    }
}

#all dmrs for only 0.7 difference
#extract top 0.7 difference for all K 
dmr_LMC_03_all <- list()
for(i in 1:8){
    T_diff <- rowMeans(T[,-i])-T[,i]
    T_all<- abs(T_diff)>0.3
    dmr_LMC_03_all[[i]]<- dmrs_red_anal[which(T_all),]
}
names(dmr_LMC_03_all)<- c("LMC1","LMC2", "LMC3", "LMC4", "LMC5", "LMC6", "LMC7","LMC8")
lapply(dmr_LMC_03_all, function(x)length(x))
#unique(as.factor(dmr_LMC_07$LMC1$hypo$SYMBOL))
saveRDS(dmr_LMC_03_all, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_03_all.rds"))


#unlist sets
dmr_LMC_list_unlist <- unlist(dmr_LMC_list)
dmr_LMC_list_promoters_unlist <- unlist(dmr_LMC_list_promoters)

#save files
dir.create(file.path(analysis.dir, "medecom","downstream_K8_lambda0.001"))
saveRDS(dmr_LMC_list_unlist, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_list_unlist.rds"))
saveRDS(dmr_LMC_list_promoters_unlist, file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_list_promoters_unlist.rds"))


