
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
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC"
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
dmrs_final_sub2 <- lapply(dmrs_final, function(x){
    x <- x[x %over% dmrs_HSC_red,]
    x
})
names(dmrs_final_sub)<- paste0(names(dmrs_final_sub), "_hierachy_DMRs_removed")
names(dmrs_final_sub2)<- paste0(names(dmrs_final), "_only_hierachy_DMRs")

lapply(dmrs_final, length)
lapply(dmrs_final_sub, length)
lapply(dmrs_final_sub2, length)
dmrs_final <- c(dmrs_final, dmrs_final_sub,dmrs_final_sub2)

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})


#export lists
dmrs_final_df<-list()
for(i in names(dmrs_final)){
    dmrs_final_df[[i]] <- as.data.frame(dmrs_final[[i]])
    dmrs_final_df[[i]]$peak_id <- 1:nrow(dmrs_final_df[[i]])

    dir.create(file.path(analysis.dir,i, "homer", "all"), recursive=TRUE)
    dir.create(file.path(analysis.dir,i, "homer", "hypo"), recursive=TRUE)
    dir.create(file.path(analysis.dir,i, "homer", "hyper"), recursive=TRUE)
    #all
    write.table(data.frame(chr=dmrs_final_df[[i]]$seqnames, start=dmrs_final_df[[i]]$start, end=dmrs_final_df[[i]]$end, id=dmrs_final_df[[i]]$peak_id, 
        notUsed=NA, strand="+"),
        file.path(analysis.dir, i, "homer", "all", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #hypo
    write.table(data.frame(chr=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]$seqnames, 
        start=dmrs_final_df[[i]][dmrs_final[[i]]$direction=="hypo",]$start, end=dmrs_final_df[[i]][dmrs_final[[i]]$direction=="hypo",]$end, 
        id=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,i, "homer", "hypo", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #hyper
    write.table(data.frame(chr=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$seqnames, 
        start=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$start, end=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$end, 
        id=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, i,"homer", "hyper", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}


#run in command line
cd /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC
chmod 777 -R ./
conda activate homer2

for file in `ls */homer/*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparse -preparsedDir ${path}/ -p 4
    echo ${path}
done

