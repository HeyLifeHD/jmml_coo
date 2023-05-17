##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Prepare Homer input and run in comman line

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
output.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"

#load data
bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#assign direction
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
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
cd /home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged
chmod 777 -R ./
conda activate homer2

for file in `ls homer/*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparsedDir ${path}/ -p 6
    echo ${path}
done


#with custom bg
cd /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged
chmod 777 -R ./homer
conda activate homer2

for file in `ls homer/*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path}_customBG -size given -preparsedDir ${path}/ -p 6 -bg /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/EpigenotypeHM/random_BG/random_bg_merged.sorted.bed
    echo ${path}
done
