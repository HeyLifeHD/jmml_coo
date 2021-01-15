##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##run homer on stratified DMRs via ESC chromHMM

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
library(data.table)
library(parallel)
#Directories
#Directories
output.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_ExplorAnal"

#load data
bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
#dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
#dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#change annotation
pheno <- pData(bsseq_all)
pheno$Tissue<- c( rep("adult_bonemarrow", 15))
pheno$Donor <- c(as.character(pheno$Patient[1:15])) 
pheno$Epigenotype <- as.character(pheno$Epigenotype)
pheno[pheno$Patient %in% c("I217", "D217"),]$Epigenotype <- "IM"
pheno$Epigenotype <- as.factor(pheno$Epigenotype)
pheno$Celltype <- c("MPP","MPP","MPP","MPP" ,"MPP","HSC", "HSC", "HSC","LMPP","LMPP","LMPP","CD45RACD90","CD45RACD90","CD45RACD90","CD45RACD90")
pheno$Sample_Type <- "tumor"
pheno$Tumor<- pheno$tumor
pData(bsseq_all) <- pheno

#new output directory
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_model_sub_repMerged"
dir.create(analysis.dir)

#prepare data for gene regulatory strat
bed_files <- paste0("/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/roadmap_tracks_jmml/",list.files("/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/roadmap_tracks_jmml"))
bed_files<-lapply(bed_files,function(x){
    x <- fread(x)
    x <- x[,1:4]
    x <- as.data.frame(x)
    colnames(x)<- c("chromosome", "start","end", "state")
    x <-makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
    x
})
names(bed_files) <- c("ESC.H1", "STRM.MRW.MSC", "BLD.CD14.PC", "BLD.CD15.PC", "BLD.CD19.PPC","BLD.CD3.PPC","BLD.CD34.PC", "BLD.CD56.PC","BLD.PER.MONUC.PC", "SPLN", "BLD.GM12878", "BLD.K562.CNCR")
#split them by state
bed_files <-lapply(bed_files, function(x){
    x <- split(x, x$state)
    x
})
#modify naming
bed_files <- lapply(bed_files, function(x){
    names(x)<- gsub("/","-", names(x))
    x
})

#summarize average methylation for each region
meth <- list()
for(i in names(bed_files)){
    meth[[i]] <- mclapply(bed_files[[i]], function(x){
        x<- bsseq::getMeth(bsseq_all, regions=x, what = "perRegion", type="raw")
    }, mc.cores=3)
}
dir.create(file.path(analysis.dir,"AverageMeth_ChromHMMStates"))
saveRDS(meth, file.path(analysis.dir,"AverageMeth_ChromHMMStates", "meth.rds"))
#summarize and combine with pheno data
av_meth <- list()
pheno <- pData(bsseq_all)
for(i in names(meth)){
    av_meth[[i]] <- mclapply(meth[[i]], function(x){
        x<- as.data.frame(colMeans(x, na.rm = TRUE))
        colnames(x)<- "AverageMethylation"
        x <- as.data.frame(cbind(x, pheno))
        x
    }, mc.cores=3)

}
saveRDS(av_meth, file.path(analysis.dir,"AverageMeth_ChromHMMStates", "av_meth.rds"))

#plot results
for(i in names(av_meth)){
    dir.create(file.path(analysis.dir,"AverageMeth_ChromHMMStates", i), recursive=TRUE)
    for(j in names(av_meth[[i]])){
        name<- gsub("/","-", j)

        #per epigenotype
        compare_means(AverageMethylation ~ Epigenotype,  data = av_meth[[i]][[j]], method = "t.test")
        my_comparisons <- list( c("HM", "LM"), c("HM", "IM"), c("LM", "IM") )
        med <- aggregate(av_meth[[i]][[j]][,"AverageMethylation",drop=FALSE], list(av_meth[[i]][[j]]$Epigenotype), median)
        ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Epigenotype.pdf")), height=4, width=4)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,1), palette = c(wildtype ="#ababab", LM = "#0058b4", IM = "#fbbb25", HM = "#c33126")
        ,add = "jitter", order=ord
        )  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(0.8,0.9,1))+stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
        #per patient
        med <- aggregate(av_meth[[i]][[j]][,"AverageMethylation",drop=FALSE], list(av_meth[[i]][[j]]$Patient), median)
        ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Patient.pdf")), height=4, width=7)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="patient", y="AverageMethylation", color="patient",ylim=c(0,1), palette =c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220")
        ,add = "jitter", order=ord
        )  +stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
        #per quadrant
        compare_means(AverageMethylation ~ tumor,  data = av_meth[[i]][[j]], method = "t.test")
        med <- aggregate(av_meth[[i]][[j]][,"AverageMethylation",drop=FALSE], list(av_meth[[i]][[j]]$tumor), median)
        ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Tumor.pdf")), height=4, width=7)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="tumor", y="AverageMethylation", color="tumor",ylim=c(0,1), palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12", 
                       tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#99a637")
        ,add = "jitter", order=ord
        )  +stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
        #per celltype
        compare_means(AverageMethylation ~ Celltype,  data = av_meth[[i]][[j]], method = "t.test")
        med <- aggregate(av_meth[[i]][[j]][,"AverageMethylation",drop=FALSE], list(av_meth[[i]][[j]]$Celltype), median)
        ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Celltype.pdf")), height=4, width=7)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="Celltype", y="AverageMethylation", color="tumor",ylim=c(0,1), palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12", 
                       tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#99a637")
        ,add = "jitter", order=ord
        )  +stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
    }
}
