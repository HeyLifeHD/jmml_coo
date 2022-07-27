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
#Directories
#Directories
odcf.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))
dmrs_final<- list(HM_vs_LM=dmrs)

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
ESC <- bed_files$ESC.H1
names(ESC)<- gsub("/","-", names(ESC))

#export lists
dmrs_final_df<-list()
for(i in names(dmrs_final)){
    dmrs_final_df[[i]]<- dmrs_final[[i]]
    dmrs_final_df[[i]]$peak_id <- 1:length(dmrs_final_df[[i]])

    dir.create(file.path(analysis.dir, "homer_ESC_Strat", "hypo"), recursive=TRUE)
    dir.create(file.path(analysis.dir, "homer_ESC_Strat", "hyper"), recursive=TRUE)

    for(j in names(ESC)){
        hypo <- dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hypo",]
        hypo_sub <- hypo[hypo %over% ESC[[j]],]
        hyper <- dmrs_final_df[[i]][dmrs_final_df[[i]]$direction=="hyper",]
        hyper_sub <- hyper[hyper %over%ESC[[j]],]
        
        dir.create(file.path(analysis.dir, "homer_ESC_Strat", "hypo", j))
        dir.create(file.path(analysis.dir, "homer_ESC_Strat", "hyper", j))
        #hypo
        write.table(data.frame(chr=seqnames(hypo_sub), 
        start=start(hypo_sub), end=end(hypo_sub), 
        id=hypo_sub$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer_ESC_Strat", "hypo", j, paste0("DMRs.bed")),
            sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
        #hyper
        write.table(data.frame(chr=seqnames(hyper_sub), 
        start=start(hyper_sub), end=end(hyper_sub), 
        id=hyper_sub$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer_ESC_Strat", "hyper", j, paste0("DMRs.bed")),
            sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    }
 }

#run in command line
cd /home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM
conda activate homer2

for file in `ls homer_ESC_Strat/*/*/DMRs.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparsedDir ${path}/ -p 6
    echo ${path}
done
