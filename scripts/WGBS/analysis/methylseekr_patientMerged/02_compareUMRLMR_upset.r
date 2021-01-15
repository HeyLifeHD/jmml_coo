#Libraries
library(MethylSeekR)
library(parallel)
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPseeker)
library(ggpubr)
#load my own data
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
methylseekr.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr"
dir.create(analysis.dir)

UMRLMRsegments.gr <- readRDS(file.path(methylseekr.dir,  paste0("UMRsLMRs", ".rds")))
head(UMRLMRsegments.gr[[i]])
lapply(UMRLMRsegments.gr, length)

#plot number of umrs and lmrs
LMR<- sapply(UMRLMRsegments.gr, function(x){
    x<- length(x[x$type=="LMR",])
    x
})
UMR<- sapply(UMRLMRsegments.gr, function(x){
    x<- length(x[x$type=="UMR",])
    x
})
numb <- as.data.frame(rbind(LMR, UMR))
numb$type <- as.character(rownames(numb))
numb <- reshape2::melt(numb)
pdf(file.path(methylseekr.dir, "numberOf_UMRLMR_all.pdf"))
ggbarplot(numb, x="variable", y="value", fill="type", 
    palette=c("black", "lightgray"),
    ylab="Number of segments")+rotate_x_text(angle = 45)+rremove("xlab")
dev.off()
#plot number of umrs
numb <- as.data.frame(t(UMR))
numb$type <-"UMR"
numb <- reshape2::melt(numb)
pdf(file.path(methylseekr.dir, "numberOf_UMR_all.pdf"))
ggbarplot(numb, x="variable", y="value", fill="type", 
    palette=c("black", "lightgray"),
    ylab="Number of segments")+rotate_x_text(angle = 45)+rremove("xlab")
dev.off()

#plot number of lmrs
numb <- as.data.frame(t(LMR))
numb$type <-"LMR"
numb <- reshape2::melt(numb)
pdf(file.path(methylseekr.dir, "numberOf_LMR_all.pdf"))
ggbarplot(numb, x="variable", y="value", fill="type", 
    palette=c("black", "lightgray"),
    ylab="Number of segments")+rotate_x_text(angle = 45)+rremove("xlab")
dev.off()

#compare overlapping umr/lmrs
RLMRsegments.gr[[i]])

#look at overlap of lmr and umr
red_UMRLMR <- reduce(unlist(GRangesList(UMRLMRsegments.gr)))
length(red_UMRLMR)
for(i in names(UMRLMRsegments.gr)){
    mcols(red_UMRLMR)[,i] <- 0
    mcols(red_UMRLMR[queryHits(findOverlaps(red_UMRLMR, UMRLMRsegments.gr[[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(red_UMRLMR))
library(UpSetR)
pdf(file.path(methylseekr.dir, "UPSET_comparison_all_UMR_LMR.pdf"), width=14, height=7)
upset(regions,nsets = ncol(regions), nintersects = 50, order.by = "freq")
dev.off()

#same for only hsc-cb as healthy reference
samples <- names(UMRLMRsegments.gr)
samples_sub <- samples[1:9]
UMRLMRsegments.gr_sub <- UMRLMRsegments.gr[samples_sub]
#look at overlap of lmr
LMRsegments.gr <- lapply(UMRLMRsegments.gr_sub , function(x){
    x <- x[x$type=="LMR",]
    x
})

red_LMR <- reduce(unlist(GRangesList(LMRsegments.gr )))
length(red_LMR)
for(i in names(LMRsegments.gr)){
    mcols(red_LMR)[,i] <- 0
    mcols(red_LMR[queryHits(findOverlaps(red_LMR, LMRsegments.gr [[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(red_LMR))
library(UpSetR)
pdf(file.path(methylseekr.dir, "UPSET_comparison_sub_LMR.pdf"), width=14, height=7)
upset(regions,nsets = ncol(regions), nintersects = 50, order.by = "freq")
dev.off()


#look at overlap of umr
UMRsegments.gr <- lapply(UMRLMRsegments.gr_sub , function(x){
    x <- x[x$type=="UMR",]
    x
})

red_UMR <- reduce(unlist(GRangesList(UMRsegments.gr )))
length(red_UMR)
for(i in names(UMRsegments.gr)){
    mcols(red_UMR)[,i] <- 0
    mcols(red_UMR[queryHits(findOverlaps(red_UMR, UMRsegments.gr [[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(red_UMR))
library(UpSetR)
pdf(file.path(methylseekr.dir, "UPSET_comparison_sub_UMR.pdf"), width=14, height=7)
upset(regions,nsets = ncol(regions), nintersects = 50, order.by = "freq")
dev.off()



#look at overlap of lmr and umr
red_UMRLMR <- reduce(unlist(GRangesList(UMRLMRsegments.gr_sub )))
length(red_UMRLMR)
for(i in names(UMRLMRsegments.gr_sub )){
    mcols(red_UMRLMR)[,i] <- 0
    mcols(red_UMRLMR[queryHits(findOverlaps(red_UMRLMR, UMRLMRsegments.gr_sub [[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(red_UMRLMR))
library(UpSetR)
pdf(file.path(methylseekr.dir, "UPSET_comparison_sub_UMR_LMR.pdf"), width=14, height=7)
upset(regions,nsets = ncol(regions), nintersects = 50, order.by = "freq")
dev.off()


#look at overlap of lmr
LMRsegments.gr <- lapply(UMRLMRsegments.gr_sub , function(x){
    x <- x[x$type=="LMR",]
    x
})

red_LMR <- reduce(unlist(GRangesList(LMRsegments.gr )))
length(red_LMR)
for(i in names(LMRsegments.gr)){
    mcols(red_LMR)[,i] <- 0
    mcols(red_LMR[queryHits(findOverlaps(red_LMR, LMRsegments.gr [[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(red_LMR))
library(UpSetR)
pdf(file.path(methylseekr.dir, "UPSET_comparison_sub_LMR.pdf"), width=14, height=7)
upset(regions,nsets = ncol(regions), nintersects = 50, order.by = "freq")
dev.off()


#look at overlap of umr
UMRsegments.gr <- lapply(UMRLMRsegments.gr_sub , function(x){
    x <- x[x$type=="UMR",]
    x
})

red_UMR <- reduce(unlist(GRangesList(UMRsegments.gr )))
length(red_UMR)
for(i in names(UMRsegments.gr)){
    mcols(red_UMR)[,i] <- 0
    mcols(red_UMR[queryHits(findOverlaps(red_UMR, UMRsegments.gr [[i]])),])[,i]  <- 1
}
regions <- as.data.frame(mcols(red_UMR))
library(UpSetR)
pdf(file.path(methylseekr.dir, "UPSET_comparison_sub_UMR.pdf"), width=14, height=7)
upset(regions,nsets = ncol(regions), nintersects = 50, order.by = "freq")
dev.off()
