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
odcf.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_ExplorAnal"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub.rds"))


#prepare data for gene regulatory strat
bed_files <- paste0("/omics/groups/OE0219/internal/jmmlc_pbat/data/roadmap_tracks_jmml/",list.files("/omics/groups/OE0219/internal/jmmlc_pbat/data/roadmap_tracks_jmml"))
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
    names(x)<- gsub("/","-", names(ESC))
    x
})

#summarize average methylation for each region
meth <- list()
for(i in names(bed_files)){
    meth[[i]] <- mclapply(bed_files[[i]], function(x){
        x<- getMeth(bsseq_all, regions=x, what = "perRegion", type="raw")
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
    dir.create(file.path(analysis.dir,"AverageMeth_ChromHMMStates", i))
    for(j in names(av_meth[[i]])){
        name<- gsub("/","-", j)

        #per epigenotype
        my_comparisons <- list( c("HM", "LM") )
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Epigenotype.pdf")), height=4, width=4)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,1), palette = c(HM = "#b5241c", LM = "#0041a5")
        ,add = "jitter"
        )  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
        #per patient
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Patient.pdf")), height=4, width=7)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="patient", y="AverageMethylation", color="patient",ylim=c(0,1), palette =c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                            D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000")
        ,add = "jitter"
        )  +stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
        #per quadrant
        pdf(file.path(analysis.dir,"AverageMeth_ChromHMMStates",i, paste0("AverageMeth_",name,"_Tumor.pdf")), height=4, width=7)
        print(ggpubr::ggboxplot(av_meth[[i]][[j]], x="tumor", y="AverageMethylation", color="tumor",ylim=c(0,1), palette =c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d")
        ,add = "jitter"
        )  +stat_compare_means(label.y = .5) +rremove("legend"))
        dev.off()
    }
}
