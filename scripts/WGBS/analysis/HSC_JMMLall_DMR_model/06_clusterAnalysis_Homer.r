
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
library(LOLA)

#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final <- readRDS(dmrs_final,file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_clust6_annotattion.rds"))


#export lists
dmrs_final_df<-list()
for(i in names(dmrs_final)){
    dmrs_final_df<- split(dmrs_final[[i]], dmrs_final[[i]]$cluster)
    names(dmrs_final_df) <- paste0("cluster_", names(dmrs_final_df))
    dmrs_final_df <- lapply(dmrs_final_df, function(x){
        x <- as.data.frame(x)
        x$peak_id <- 1:nrow(x)
        x
    })
    for(j in names(dmrs_final_df)){
        dir.create(file.path(analysis.dir,i, "homer_cluster", j), recursive=TRUE)
        write.table(data.frame(chr=dmrs_final_df[[j]]$seqnames, 
        start=dmrs_final_df[[j]]$start, end=dmrs_final_df[[j]]$end, 
        id=dmrs_final_df[[j]]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,i, "homer_cluster", j, paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    }    
}

#run in command line
cd /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC
chmod 777 -R ./
conda activate homer2

for file in `ls */homer_cluster/*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparse -preparsedDir ${path}/ -p 3
    echo ${path}
done

