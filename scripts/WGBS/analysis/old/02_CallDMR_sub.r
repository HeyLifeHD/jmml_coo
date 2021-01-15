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
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/DMR_sub"

dir.create(analysis.dir)
#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))
#subset outliers
bsseq_all$tumor11_JMMLC_D129 <- NULL
#get pheno
pheno <- as.data.frame(pData(bsseq_all))
#relevel pheno so it fits
pheno$Epigenotype <- as.factor(pheno$Epigenotype )
pheno$Epigenotype <-factor(pheno$Epigenotype, levels = c("LM", "HM"))
pheno$Genotype <- as.factor(pheno$Genotype)
pheno$Genotype <-factor(pheno$Genotype, levels = c("neg", "KRAS", "PTPN11"))
pData(bsseq_all) <- pheno
#create model
modelmatrix = model.matrix(~Genotype + Epigenotype + tumor + Epigenotype:tumor, pheno)
is.fullrank(modelmatrix)
#run linear dss model
dml_fit_complete<-  DMLfit.multiFactor(bsseq_all, design=pheno, formula=~Genotype + Epigenotype + tumor + Epigenotype:tumor,  smoothing = TRUE)
dir.create(file.path(analysis.dir))
saveRDS(dml_fit_complete, file.path(analysis.dir ,"dml_fit_complete_inter.rds"))

#test all different coefficients
DMLtest_complete <- list()
for(i in colnames(dml_fit_complete$X)){
    DMLtest_complete[[i]]<-  DMLtest.multiFactor(dml_fit_complete, coef=i)
    print(paste0("coefficient ", i, " tested"))
}
saveRDS(DMLtest_complete, file.path(analysis.dir ,"dml_test_complete_inter.rds"))

#call dmr from dml test
DMRtest_complete <- list()
for(i in colnames(dml_fit_complete$X)){
    DMRtest_complete[[i]]<-  callDMR(DMLtest_complete[[i]],
        p.threshold=0.05,#delta=0.1, 
        minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
    print(paste0("DMR ", i, " called"))
}
saveRDS(DMRtest_complete, file.path(analysis.dir ,"dmr_test_complete_inter.rds"))
lapply(DMRtest_complete   , function(x)dim(x)) 
#make granges
dmrs_gr<- list()
dir.create(file.path(analysis.dir, "tracks"))
for(i in names(DMRtest_complete)){
dmrs_gr[[i]]<- GRanges(
  seqnames = DMRtest_complete[[i]]$chr,
  ranges = IRanges(start =  DMRtest_complete[[i]]$start,
                   end =  DMRtest_complete[[i]]$end
  )
)
mcols(dmrs_gr[[i]])<- DMRtest_complete[[i]][,5:6]
export.bed(dmrs_gr[[i]],file.path(analysis.dir, "tracks", paste0("dmrs_",i, "_inter.bed")))
}
saveRDS(dmrs_gr, file.path(analysis.dir ,"dmrs_gr_inter.rds"))
dmrs_gr <- readRDS(file.path(analysis.dir ,"dmrs_gr_inter.rds"))

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
saveRDS(dmrs_gr, file.path(analysis.dir ,"dmrs_gr_anno_inter.rds"))


#reduce dmrs
dmrs_gr_red <- GRangesList( dmrs_gr[[2]], dmrs_gr[[3]], dmrs_gr[[4]], dmrs_gr[[5]], dmrs_gr[[6]], dmrs_gr[[7]], dmrs_gr[[8]], dmrs_gr[[9]], dmrs_gr[[10]])
dmrs_gr_red <- unlist(dmrs_gr_red)
dmrs_gr_red <- reduce(dmrs_gr_red)
length(dmrs_gr_red)

#subset unique dmrs
dmr_all_meth <- bsseq::getMeth(bsseq_all, regions=dmrs_gr_red, what="perRegion", type="raw")
mcols(dmrs_gr_red)<- dmr_all_meth

#annotate all dmrs
dmrs_red_anno <- annotatePeak(peak= dmrs_gr_red, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
dir.create(file.path(analysis.dir,"all"))

#visualize dmrs annotation
dir.create(file.path(analysis.dir,"all",  "anno"))
pdf(file.path(analysis.dir,"all",  "anno", paste0("all", "_AnnoDist_Pie.pdf")))
plotAnnoPie(dmrs_red_anno)
dev.off()

pdf(file.path(analysis.dir,"all","anno", paste0("all", "_AnnoDist_Upset.pdf")))
upsetplot(dmrs_red_anno, vennpie=TRUE)
dev.off()

pdf(file.path(analysis.dir,"all","anno", paste0("all", "_AnnoDist_Upset_noVenn.pdf")))
upsetplot(dmrs_red_anno , vennpie=FALSE)
dev.off()

pdf(file.path(analysis.dir, "all", "anno",paste0("all", "_DistanceToTSS.pdf")))
plotDistToTSS(dmrs_red_anno, 
              title="Distribution of DMR relative to TSS")
dev.off()  

#extract annotated granges
dmrs_red_anno_df <- as.data.frame(dmrs_red_anno)
mcols(dmrs_gr_red) <- dmrs_red_anno_df[,-c(1:5)]
saveRDS(dmrs_gr_red, file.path(analysis.dir, "dmrs_gr_reduced_anno_inter.rds"))
