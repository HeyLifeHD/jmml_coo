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
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
methylseekr.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr"
dir.create(analysis.dir)

#load data
bsseq_merged <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_PatientMerged_CelltypeMerged.rds"))
UMRLMRsegments.gr <- readRDS(file.path(methylseekr.dir,  paste0("UMRsLMRs", ".rds")))

#reduce ranges
red_UMRLMR<- reduce(unlist(GRangesList(UMRLMRsegments.gr)))
length(red_UMRLMR)

#same only for umr
UMRsegments <- lapply(UMRLMRsegments.gr, function(x){
    x <- x[which(x$type=="UMR"),]
    x
})
red_UMR<- reduce(unlist(GRangesList(UMRsegments)))
length(red_UMR)

#same for lmr
LMRsegments <- lapply(UMRLMRsegments.gr, function(x){
    x <- x[which(x$type=="LMR"),]
    x
})
red_LMR<- reduce(unlist(GRangesList(LMRsegments)))
length(red_LMR)

#export lists UMRLMR
dmrs_final_df<-list()
for(i in names(UMRLMRsegments.gr)){
    dmrs_final_df[[i]] <- UMRLMRsegments.gr[[i]]
    dmrs_final_df[[i]]$peak_id <- 1:length(dmrs_final_df[[i]])
    dir.create(file.path(methylseekr.dir, i, "homer"), recursive=TRUE)
    #all
    write.table(data.frame(chr=as.character(seqnames(dmrs_final_df[[i]])), start=as.character(start(dmrs_final_df[[i]])), end=as.character(end(dmrs_final_df[[i]])), id=dmrs_final_df[[i]]$peak_id, 
        notUsed=NA, strand="+"),
        file.path(methylseekr.dir, i,  "homer", paste0("UMRLMR.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    }
#and bg
dmrs_final_df <- red_UMRLMR
dmrs_final_df$peak_id <- 1:length(dmrs_final_df)
write.table(data.frame(chr=as.character(seqnames(dmrs_final_df)), start=as.character(start(dmrs_final_df)), end=as.character(end(dmrs_final_df)), id=dmrs_final_df$peak_id, 
    notUsed=NA, strand="+"),
    file.path(methylseekr.dir,  paste0("UMRLMR_red.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)


#export lists UMR
dmrs_final_df<-list()
for(i in names(UMRsegments)){
    dmrs_final_df[[i]] <- UMRsegments[[i]]
    dmrs_final_df[[i]]$peak_id <- 1:length(dmrs_final_df[[i]])
    dir.create(file.path(methylseekr.dir, i, "homer"), recursive=TRUE)
    #all
    write.table(data.frame(chr=as.character(seqnames(dmrs_final_df[[i]])), start=as.character(start(dmrs_final_df[[i]])), end=as.character(end(dmrs_final_df[[i]])), id=dmrs_final_df[[i]]$peak_id, 
        notUsed=NA, strand="+"),
        file.path(methylseekr.dir, i,  "homer", paste0("UMR.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    }
#and bg
dmrs_final_df <- red_UMR
dmrs_final_df$peak_id <- 1:length(dmrs_final_df)
write.table(data.frame(chr=as.character(seqnames(dmrs_final_df)), start=as.character(start(dmrs_final_df)), end=as.character(end(dmrs_final_df)), id=dmrs_final_df$peak_id, 
    notUsed=NA, strand="+"),
    file.path(methylseekr.dir,  paste0("UMR_red.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

#export lists LMR
dmrs_final_df<-list()
for(i in names(LMRsegments)){
    dmrs_final_df[[i]] <- LMRsegments[[i]]
    dmrs_final_df[[i]]$peak_id <- 1:length(dmrs_final_df[[i]])
    #all
    write.table(data.frame(chr=as.character(seqnames(dmrs_final_df[[i]])), start=as.character(start(dmrs_final_df[[i]])), end=as.character(end(dmrs_final_df[[i]])), id=dmrs_final_df[[i]]$peak_id, 
        notUsed=NA, strand="+"),
        file.path(methylseekr.dir, i,  "homer", paste0("LMR.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    }
#and bg
dmrs_final_df <- red_LMR
dmrs_final_df$peak_id <- 1:length(dmrs_final_df)
write.table(data.frame(chr=as.character(seqnames(dmrs_final_df)), start=as.character(start(dmrs_final_df)), end=as.character(end(dmrs_final_df)), id=dmrs_final_df$peak_id, 
    notUsed=NA, strand="+"),
    file.path(methylseekr.dir,  paste0("LMR_red.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)


#run in command line
cd /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr/
chmod 777  -R ./
conda activate homer2

for file in `ls */homer/UMRLMR.bed`
do
    echo ${file}
    path=`dirname ${file}`
    mkdir ${path}/UMR_vs_redUMRLMR/
    findMotifsGenome.pl ${file} hg19 ${path}/UMR_vs_redUMRLMR/ -size given -preparsedDir ${path}/ -p 6 \
    -bg UMR_red.bed
done

cd /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr/
conda activate homer2

for file in `ls */homer/UMR.bed`
do
    echo ${file}
    path=`dirname ${file}`
    mkdir ${path}/UMR_vs_redUMR/
    findMotifsGenome.pl ${file} hg19 ${path}/UMR_vs_redUMR/ -size given -preparsedDir ${path}/ -p 6 \
    -bg UMR_red.bed
done

cd /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr/
conda activate homer2

for file in `ls */homer/LMR.bed`
do
    echo ${file}
    path=`dirname ${file}`
    mkdir ${path}/LMR_vs_redLMR/
    findMotifsGenome.pl ${file} hg19 ${path}/LMR_vs_redLMR/ -size given -preparsedDir ${path}/ -p 6 \
    -bg LMR_red.bed
done

