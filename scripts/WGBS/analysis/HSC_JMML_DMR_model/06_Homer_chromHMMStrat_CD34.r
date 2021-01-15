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
library(data.table)

#load my own data: Epigenotype vs HSC_cb
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

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
dmrs_final$all <- dmrs_red
mcols(dmrs_final$all)$direction <- "hypo"


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
CD34 <- bed_files$BLD.CD34.PC
names(CD34)<- gsub("/","-", names(CD34))

#export lists
dmrs_final_df<-list()
for(i in names(dmrs_final)){
    dmrs_final_df[[i]]<- dmrs_final[[i]]
    dmrs_final_df[[i]]$peak_id <- 1:length(dmrs_final_df[[i]])
    dmrs_final_df[[i]]$direction <- ifelse(dmrs_final_df[[i]]$diff.Methy >0, "hyper", "hypo")
    dir.create(file.path(analysis.dir,i, "homer_CD34_Strat", "hypo"), recursive=TRUE)
    dir.create(file.path(analysis.dir,i,  "homer_CD34_Strat", "hyper"), recursive=TRUE)

    for(j in names(CD34)){
        hypo <- dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]
        hypo_sub <- hypo[hypo %over% CD34[[j]],]
        hyper <- dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]
        hyper_sub <- hyper[hyper %over%CD34[[j]],]
        
        dir.create(file.path(analysis.dir,i,  "homer_CD34_Strat", "hypo", j))
        dir.create(file.path(analysis.dir,i,  "homer_CD34_Strat", "hyper", j))
        #hypo
        write.table(data.frame(chr=seqnames(hypo_sub), 
        start=start(hypo_sub), end=end(hypo_sub), 
        id=hypo_sub$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,i,  "homer_CD34_Strat", "hypo", j, paste0("DMRs.bed")),
            sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
        #hyper
        write.table(data.frame(chr=seqnames(hyper_sub), 
        start=start(hyper_sub), end=end(hyper_sub), 
        id=hyper_sub$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, i, "homer_CD34_Strat", "hyper", j, paste0("DMRs.bed")),
            sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    }
 }

#run in command line
cd /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC
chmod 777 -R ./
conda activate homer2

for file in `ls -r */homer_CD34_Strat/*/*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparsedDir ${path}/ -p 6
    echo ${path}
done

