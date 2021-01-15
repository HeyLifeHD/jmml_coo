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
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

#load seperated bsseq data
#bsseq<- readRDS(file.path(output.dir, "bsseq.rds"))
#dmrs
dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))
dmrs_red_anal <- readRDS(file.path(analysis.dir, "medecom", "dmrs_red_anal.rds"))
#load medecom results 
medecom.result <- readRDS(file.path(analysis.dir, "medecom", "allDMRred_Medecom_results.rds"))
dmr_LMC_list_unlist <- readRDS(file.path(analysis.dir, "medecom","downstream_combined_K8_lambda0.001", "dmr_LMC_list_unlist.rds"))

#create bg
dir.create(file.path(analysis.dir,"medecom", "downstream_combined_K8_lambda0.001", "homer"))
dmrs_red_anal_df <- as.data.frame(dmrs_red_anal)
dmrs_red_anal_df$peak_id <- 1:nrow(dmrs_red_anal_df)
write.table(data.frame(chr=dmrs_red_anal_df$seqnames, 
        start=dmrs_red_anal_df$start, end=dmrs_red_anal_df$end, 
        id=dmrs_red_anal_df$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,"medecom", "downstream_combined_K8_lambda0.001", "homer", "dmrs_red_anal_df.txt"),
        ,sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

#export lists for hypo 
dmrs_final_df<-list()
dmr_LMC_list_unlist_hypo <- dmr_LMC_list_unlist[grep("hypo", names(dmr_LMC_list_unlist))]
for(i in names(dmr_LMC_list_unlist_hypo)){
    dmrs_final_df[[i]] <- as.data.frame(dmr_LMC_list_unlist_hypo[[i]])
    dmrs_final_df[[i]]$peak_id <- 1:nrow(dmrs_final_df[[i]])
    name <- strsplit(i, ".", fixed=TRUE)
    name <- sapply(name, function(x)paste0(x[1], x[2]))
    dir.create(file.path(analysis.dir,"medecom", "downstream_combined_K8_lambda0.001", "homer",name, "hypo"), recursive=TRUE)
    #hypo
    write.table(data.frame(chr=dmrs_final_df[[i]]$seqnames, 
        start=dmrs_final_df[[i]]$start, end=dmrs_final_df[[i]]$end, 
        id=dmrs_final_df[[i]]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,"medecom", "downstream_combined_K8_lambda0.001", "homer",name, "hypo", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}

#export lists for hyper
dmrs_final_df<-list()
dmr_LMC_list_unlist_hyper <- dmr_LMC_list_unlist[grep("hyper", names(dmr_LMC_list_unlist))]
for(i in names(dmr_LMC_list_unlist_hyper)){
    dmrs_final_df[[i]] <- as.data.frame(dmr_LMC_list_unlist_hyper[[i]])
    dmrs_final_df[[i]]$peak_id <- 1:nrow(dmrs_final_df[[i]])
    name <- strsplit(i, ".", fixed=TRUE)
    name <- sapply(name, function(x)paste0(x[1], x[2]))
    dir.create(file.path(analysis.dir,"medecom", "downstream_combined_K8_lambda0.001", "homer",name, "hyper"), recursive=TRUE)
    #hypo
    write.table(data.frame(chr=dmrs_final_df[[i]]$seqnames, 
        start=dmrs_final_df[[i]]$start, end=dmrs_final_df[[i]]$end, 
        id=dmrs_final_df[[i]]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,"medecom", "downstream_combined_K8_lambda0.001", "homer",name, "hyper", paste0("DMRs.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}

#run in command line
cd /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM/medecom/downstream_combined_K8_lambda0.001/homer
conda activate homer2

for file in `ls */*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparsedDir ${path}/ -p 6
    echo ${path}
done

for file in `ls */*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path}_bg -size given -preparsedDir ${path}/ -p 6 \
    -bg dmrs_red_anal_df.txt
done

#get motif locations
cd /icgc/dkfzlsdf/analysis/C010/cancer_microenvironment/TAM/WGBS/analysis/BMDM_TAM_MG_onlyPBAT/DMR/medecom/downstream_combined_K8_lambda0.001/homer
#source activate homer
conda activate homer2

for file in `ls -1 | grep "dmr_LMC"`
do
    mkdir ${file}/hypo_bg/motif_location
    mkdir ${file}/hyper_bg/motif_location

    for motif in `ls -1 ${file}/hypo_bg/knownResults/*.motif`
    do 
        motif_name=`basename $motif`
        echo ${motif_name}
        echo ${file}
        findMotifsGenome.pl ${file}/hypo/DMRs.bed \
        mm10 ${file}/hypo_bg  -find ${file}/hypo_bg/knownResults/${motif_name} > ${file}/hypo_bg/motif_location/${motif_name}_outputfile.txt -p 40
    done

    for motif in `ls -1 ${file}/hyper_bg/knownResults/*.motif`
    do 
        motif_name=`basename $motif`
        echo ${motif_name}
        echo ${file}
        findMotifsGenome.pl ${file}/hyper/DMRs.bed \
        mm10 ${file}/hyper_bg  -find ${file}/hyper_bg/knownResults/${motif_name} > ${file}/hyper_bg/motif_location/${motif_name}_outputfile.txt -p 40
    done


done

