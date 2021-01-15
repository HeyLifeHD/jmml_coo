##### Joschka Hey 
##### 
#### BSseq and DSS
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
library(limma)
library(ChIPseeker)

#directories
odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/DMR"

dir.create(analysis.dir)
#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))

#subset outliers
bsseq_all_sub <- bsseq_all[]
#get pheno
pheno <- pData(bsseq_all)
#create model
modelmatrix = model.matrix(~patient + tumor + Epigenotype  + patient:tumor + patient:Epigenotype + tumor:Epigenotype , pheno)
is.fullrank(modelmatrix)
modelmatrix = model.matrix(~Genotype + tumor + Epigenotype, pheno)
is.fullrank(modelmatrix)
modelmatrix = model.matrix(~patient , pheno)
is.fullrank(modelmatrix)
modelmatrix = model.matrix(~Genotype + Epigenotype + tumor + Epigenotype:tumor, pheno)
is.fullrank(modelmatrix)
#run linear dss model
dml_fit_complete<-  DMLfit.multiFactor(bsseq_all, design=pheno, formula=~Genotype + Epigenotype + tumor,  smoothing = TRUE)
dir.create(file.path(analysis.dir,"medecom_dss"))
saveRDS(dml_fit_complete, file.path(analysis.dir,"medecom_dss" ,"dml_fit_complete.rds"))

#test all different coefficients
DMLtest_complete <- list()
for(i in colnames(dml_fit_complete$X)){
    DMLtest_complete[[i]]<-  DMLtest.multiFactor(dml_fit_complete, coef=i)
    print(paste0("coefficient ", i, " tested"))
}
saveRDS(DMLtest_complete, file.path(analysis.dir,"medecom_dss" ,"dml_test_complete.rds"))

#call dmr from dml test
DMRtest_complete <- list()
for(i in colnames(dml_fit_complete$X)){
    DMRtest_complete[[i]]<-  callDMR(DMLtest_complete[[i]],
        p.threshold=0.05,#delta=0.1, 
        minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
    print(paste0("DMR ", i, " called"))
}
saveRDS(DMRtest_complete, file.path(analysis.dir,"medecom_dss" ,"dmr_test_complete.rds"))
lapply(DMRtest_complete   , function(x)dim(x))

DMRtest_complete <- readRDS(file.path(analysis.dir,"medecom_dss" ,"dmr_test_complete.rds"))
#make granges
dmrs_gr<- list()
dir.create(file.path(analysis.dir,"medecom_dss", "tracks"))
for(i in names(DMRtest_complete)){
dmrs_gr[[i]]<- GRanges(
  seqnames = DMRtest_complete[[i]]$chr,
  ranges = IRanges(start =  DMRtest_complete[[i]]$start,
                   end =  DMRtest_complete[[i]]$end
  )
)
mcols(dmrs_gr[[i]])<- DMRtest_complete[[i]][,5:6]
export.bed(dmrs_gr[[i]],file.path(analysis.dir,"medecom_dss", "tracks", paste0("dmrs_",i, ".bed")))
}
saveRDS(dmrs_gr, file.path(analysis.dir,"medecom_dss" ,"dmrs_gr.rds"))


#Annotate and plot dmrs 
#Load TxDb file
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

dmrs_anno<- list()
dmrs_anno_df<- list()
#plot files and annotate
for(i in names(dmrs_gr)){
#annotate dmrs
dmrs_anno[[i]] <- annotatePeak(peak= dmrs_gr[[i]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

dir.create(file.path(analysis.dir, "medecom_dss",i))
#visualize dmrs annotation
dir.create(file.path(analysis.dir, "medecom_dss",i,  "anno"))
#pdf(file.path(analysis.dir, "medecom_dss",i,  "anno", paste0(i, "_AnnoDist_Pie.pdf")))
#print(plotAnnoPie(dmrs_anno[[i]] ))
#dev.off()
#pdf(file.path(analysis.dir, "medecom_dss",i,"anno", paste0(i, "_AnnoDist_Upset.pdf")))
#print(upsetplot(dmrs_anno[[i]] , vennpie=TRUE))
#dev.off()
#pdf(file.path(analysis.dir, "medecom_dss",i,"anno", paste0(i, "_AnnoDist_Upset_noVenn.pdf")))
#print(upsetplot(dmrs_anno[[i]] , vennpie=FALSE))
#dev.off()
#pdf(file.path(analysis.dir, "medecom_dss", i, "anno",paste0(i, "_DistanceToTSS.pdf")))
#print(plotDistToTSS(dmrs_anno[[i]], 
#              title="Distribution of DMR relative to TSS"))
#dev.off()  

#extract annotated granges
dmrs_anno_df[[i]] <- as.data.frame(dmrs_anno[[i]])
mcols(dmrs_gr[[i]]) <- dmrs_anno_df[[i]][6:19]
}
saveRDS(dmrs_gr, file.path(analysis.dir,"medecom_dss" ,"dmrs_gr_anno.rds"))
dmrs_gr <- readRDS(file.path(analysis.dir,"medecom_dss" ,"dmrs_gr_anno.rds"))

temp <- dmrs_gr$GenotypePTPN11
length(temp)
temp <- temp[grep("Promoter", temp$annotation),]
paste0(head(unique(temp$SYMBOL), 500), collapse=", ")