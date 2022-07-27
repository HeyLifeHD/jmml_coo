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
dir.create(analysis.dir)

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

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
#dmrs_final$all <- dmrs_red
#mcols(dmrs_final$all)$direction <- "hypo"

#split in hypo and hyper
dmrs_final_split <- lapply(dmrs_final, function(x){
    x <- split(x, x$direction)
})

#extract methylation levels for each 
meth_dmrs<- list()
for(comp in names(dmrs_final_split)){
    meth_dmrs[[comp]] <- list()
    meth_dmrs[[comp]]$hypo <- bsseq::getMeth(bsseq_all, regions= dmrs_final_split[[comp]]$hypo, type = "raw", what=c("perRegion"))
    meth_dmrs[[comp]]$hyper <- bsseq::getMeth(bsseq_all, regions= dmrs_final_split[[comp]]$hyper, type = "raw", what=c("perRegion"))
    print(comp)
}
#get methylation median
pheno <- as.data.frame(pData(bsseq_all))
meth_dmrs_median<- list()
for(comp in names(dmrs_final_split)){
    meth_dmrs_median[[comp]] <- list()
    temp1 <- pheno
    temp1$direction <-"hypo" 
    temp1$Methylation_level <- colMedians(meth_dmrs[[comp]]$hypo, na.rm=TRUE)
    temp2 <- pheno
    temp2$direction <-"hyper" 
    temp2$Methylation_level <- colMedians(meth_dmrs[[comp]]$hyper, na.rm=TRUE)
    meth_dmrs_median[[comp]] <- rbind(temp1, temp2)
    #meth_dmrs_median[[comp]]$methylation_level_hypo <- colMedians(meth_dmrs[[comp]]$hypo, na.rm=TRUE)
    #meth_dmrs_median[[comp]]$methylation_level_hyper <- colMedians(meth_dmrs[[comp]]$hyper, na.rm=TRUE)
}

#Plot average Methylation 
dir.create(file.path(analysis.dir ,"dmr_methylation_levels"))
for(i in names(meth_dmrs_median)){
#plotting
#Donor
compare_means(Methylation_level ~ Donor ,  data = meth_dmrs_median[[i]], method = "t.test")
#med <- aggregate(meth_dmrs_median[[i]][,"Methylation_level",drop=FALSE], list(meth_dmrs_median[[i]]$Donor), median)
#ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
ord <- c("cordblood", "D117", "D129", "D217", "I217", "D213", "D123", "D360", "D124")

pdf(file.path(analysis.dir ,"dmr_methylation_levels",paste0(i, "_AvMeth_Donor.pdf")), height=4, width=5)
print(ggpubr::ggboxplot(meth_dmrs_median[[i]], x="Donor", y="Methylation_level", color="Donor", facet.by="direction",ylim=c(0,1), order=ord,
    palette =c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#Epigenotype
compare_means(Methylation_level ~ Epigenotype ,  data = meth_dmrs_median[[i]], method = "t.test")
#med <- aggregate(meth_dmrs_median[[i]][,"Methylation_level",drop=FALSE], list(meth_dmrs_median[[i]]$Epigenotype), median)
#ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
ord <- c("wildtype", "LM", "IM", "HM")

pdf(file.path(analysis.dir ,"dmr_methylation_levels",paste0(i, "_AvMeth_Epigenotype.pdf")), height=4, width=5)
print(ggpubr::ggboxplot(meth_dmrs_median[[i]], x="Epigenotype", y="Methylation_level",  facet.by="direction", ylim=c(0,1), order=ord,
    color="Epigenotype", palette =c(wildtype ="#ababab", LM = "#0058b4", IM = "#fbbb25", HM = "#c33126"),title=i, ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#genotype
compare_means(Methylation_level ~ Genotype,  data = meth_dmrs_median[[i]], method = "t.test")
med <- aggregate(meth_dmrs_median[[i]][,"Methylation_level",drop=FALSE], list(meth_dmrs_median[[i]]$Genotype), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir ,"dmr_methylation_levels",paste0(i, "_AvMeth_Genotype.pdf")), height=4, width=5)
print(ggpubr::ggboxplot(meth_dmrs_median[[i]], x="Genotype", y="Methylation_level", color="Genotype",facet.by="direction", ylim=c(0,1), order=ord,
    palette =c(wildtype ="#ababab", neg = "#529a51", KRAS = "#99a637", PTPN11 = "#007458"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#tumor
compare_means(Methylation_level ~ Tumor,  data = meth_dmrs_median[[i]], method = "t.test")
med <- aggregate(meth_dmrs_median[[i]][,"Methylation_level",drop=FALSE], list(meth_dmrs_median[[i]]$Tumor), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir ,"dmr_methylation_levels",paste0(i, "_AvMeth_tumor.pdf")), height=4, width=5)
print(ggpubr::ggboxplot(meth_dmrs_median[[i]], x="Tumor", y="Methylation_level", color="Tumor",facet.by="direction",ylim=c(0,1), order=ord,
    palette =c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#c6a27f", 
                       HSC_adult ="#252525", MPP_adult = "#737373", CMP_adult = "#ffb86f", GMP_adult = "#e27e37", MEP_adult="#d2624a",
                       HSC_CB ="#252525"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#Celltype
compare_means(Methylation_level ~ Celltype,  data = meth_dmrs_median[[i]], method = "t.test")
med <- aggregate(meth_dmrs_median[[i]][,"Methylation_level",drop=FALSE], list(meth_dmrs_median[[i]]$Celltype), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir ,"dmr_methylation_levels",paste0(i, "_AvMeth_Celltype.pdf")), height=4, width=5)
print(ggpubr::ggboxplot(meth_dmrs_median[[i]], x="Celltype", y="Methylation_level", color="Celltype",facet.by="direction",ylim=c(0,1), order=ord,
    palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#Sample_Type
compare_means(Methylation_level ~ Sample_Type,  data = meth_dmrs_median[[i]], method = "t.test")
med <- aggregate(meth_dmrs_median[[i]][,"Methylation_level",drop=FALSE], list(meth_dmrs_median[[i]]$Sample_Type), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir ,"dmr_methylation_levels",paste0(i, "_AvMeth_Sample_Type.pdf")), height=4, width=5)
print(ggpubr::ggboxplot(meth_dmrs_median[[i]], x="Sample_Type", y="Methylation_level", color="Sample_Type",facet.by="direction",ylim=c(0,1), order=ord,
    palette =c(normal = "#ababab", tumor = "#99a637"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter") +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
    
}


