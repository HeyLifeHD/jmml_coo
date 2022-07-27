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
#Directories
input.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/DMR_sub"
input.dir <- "icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200512_DMR_HM_vs_LM"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

dmrs_final<- list(HM_vs_LM=dmrs)
#export lists
dmrs_final_df<-list()
for(i in names(dmrs_final)){
    dmrs_final_df[[i]] <- as.data.frame(dmrs_final[[i]])
    dmrs_final_df[[i]]$peak_id <- 1:nrow(dmrs_final_df[[i]])

    dir.create(file.path(analysis.dir, "homer", "all"), recursive=TRUE)
    dir.create(file.path(analysis.dir, "homer", "hypo"), recursive=TRUE)
    dir.create(file.path(analysis.dir, "homer", "hyper"), recursive=TRUE)
    #all
    write.table(data.frame(chr=dmrs_final_df[[i]]$seqnames, start=dmrs_final_df[[i]]$start, end=dmrs_final_df[[i]]$end, id=dmrs_final_df[[i]]$peak_id, 
        notUsed=NA, strand="+"),
        file.path(analysis.dir,  "homer", "all", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #hypo
    write.table(data.frame(chr=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]$seqnames, 
        start=dmrs_final_df[[i]][dmrs_final[[i]]$direction=="hypo",]$start, end=dmrs_final_df[[i]][dmrs_final[[i]]$direction=="hypo",]$end, 
        id=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer", "hypo", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #hyper
    write.table(data.frame(chr=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$seqnames, 
        start=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$start, end=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$end, 
        id=dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer", "hyper", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}

#run in command line
cd /home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200512_DMR_HM_vs_LM
conda activate homer2

for file in `ls homer/*/DMRs.bed`#
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparsedDir ${path}/ -p 6
    echo ${path}
done
