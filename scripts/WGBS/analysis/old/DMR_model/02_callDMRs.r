##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##call HM vs LM dmrs

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
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200608_DMR_model"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all_snpfil_sub_cov.rds"))

#get pheno
pheno <- as.data.frame(pData(bsseq_all))
#design matrix
modelmatrix = model.matrix(~Genotype + Epigenotype + tumor , pheno)
modelmatrix = model.matrix(~ Epigenotype + tumor , pheno)

is.fullrank(modelmatrix)
#run linear dss model
dml_fit_complete<-  DMLfit.multiFactor(bsseq_all, design=pheno, formula=~Genotype + Epigenotype + tumor,  smoothing = TRUE)
dir.create(file.path(analysis.dir))
saveRDS(dml_fit_complete, file.path(analysis.dir ,"dml_fit_complete.rds"))

#test all different coefficients
DMLtest_complete <- list()
for(i in colnames(dml_fit_complete$X)){
    DMLtest_complete[[i]]<-  DMLtest.multiFactor(dml_fit_complete, coef=i)
    print(paste0("coefficient ", i, " tested"))
}
saveRDS(DMLtest_complete, file.path(analysis.dir ,"dml_test_complete.rds"))
DMLtest_complete <-readRDS(file.path(analysis.dir ,"dml_test_complete.rds"))

#call dmr from dml test
DMRtest_complete <- list()
for(i in colnames(dml_fit_complete$X)){
    DMRtest_complete[[i]]<-  callDMR(DMLtest_complete[[i]],
        p.threshold=0.05,#delta=0.1, 
        minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
    print(paste0("DMR ", i, " called"))
}
saveRDS(DMRtest_complete, file.path(analysis.dir ,"dmr_test_complete.rds"))
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
export.bed(dmrs_gr[[i]],file.path(analysis.dir, "tracks", paste0("dmrs_",i, ".bed")))
}
saveRDS(dmrs_gr, file.path(analysis.dir ,"dmrs_gr.rds"))


#subset dmrs based on coverage
cov_drms <- list()
for(comp in names(dmrs_gr)){
  cov_drms[[comp]] <- bsseq::getCoverage(bsseq_all, regions= dmrs_gr[[comp]], type = "Cov", what=c("perRegionAverage"))
  print(comp)
}

#define groups
contrasts <- list(GenotypeKRAS=list(KRAS=rownames(pheno[pheno$Genotype=="KRAS",]), remains=rownames(pheno[pheno$Genotype=="neg",])),
    GenotypePTPN11=list(PTPN11=rownames(pheno[pheno$Genotype=="PTPN11",]), remains=rownames(pheno[pheno$Genotype=="neg",])),
    EpigenotypeHM=list(HM=rownames(pheno[pheno$Epigenotype=="HM",]), LM=rownames(pheno[pheno$Epigenotype=="LM",])),
    tumortumor01=list(tumor01=rownames(pheno[pheno$tumor=="tumor01",]), remains=rownames(pheno[pheno$tumor=="tumor00",])),
    tumortumor10=list(tumor10=rownames(pheno[pheno$tumor=="tumor10",]), remains=rownames(pheno[pheno$tumor=="tumor00",])),
    tumortumor11=list(tumor11=rownames(pheno[pheno$tumor=="tumor11",]), remains=rownames(pheno[pheno$tumor=="tumor00",])))

#find out which dmrs have a coverage of at least 3 in at least two of the groups which are being compared
keepLoci <- list()
for(i in names(contrasts)){
    keepLoci[[i]] <- which(rowSums(cov_drms[[i]][, contrasts[[i]][[1]]] >= 5) >= 3 &
                     rowSums(cov_drms[[i]][, contrasts[[i]][[2]]] >= 5) >= 3)
}
temp <- as.vector(c(0,as.numeric(lapply(keepLoci, function(x)length(x)))))/as.vector(as.numeric(lapply(dmrs_gr, function(x)length(x))))
names(temp)<-names(cov_drms)
temp
temp <- as.vector(as.numeric(lapply(dmrs_gr, function(x)length(x))))
names(temp)<-names(cov_drms)
temp

#subset dmrs
dmrs_gr_sub <- list()
for(i in names(keepLoci)){
    dmrs_gr_sub[[i]]<- dmrs_gr[[i]][keepLoci[[i]],]
    dmrs_gr_sub[[i]] <- dmrs_gr_sub[[i]]
}
lapply(dmrs_gr_sub, function(x)length(x))
saveRDS(dmrs_gr_sub, file.path(analysis.dir, "sig_dmrs_5in3_sub.rds"))

#add methylation difference information
meth_dmrs<- list()
for(comp in names(dmrs_gr)){
  meth_dmrs[[comp]] <- bsseq::getMeth(bsseq_all, regions= dmrs_gr[[comp]], type = "raw", what=c("perRegion"))
  print(comp)
}

#take difference of rowmean of different groups defined in contrasts
meth_dmrs_diff <- list()
for(i in names(contrasts)){
    meth_dmrs_diff[[i]] <-rowMeans(meth_dmrs[[i]][, contrasts[[i]][[1]]], na.rm=TRUE) - rowMeans(meth_dmrs[[i]][, contrasts[[i]][[2]]], na.rm=TRUE)
    dmrs_gr[[i]]$diff.Methy <- meth_dmrs_diff[[i]]
    mcols(dmrs_gr[[i]]) <- cbind(mcols(dmrs_gr[[i]]), meth_dmrs[[i]])
}

#subset dmrs with a abs(difference) >0.1 
dmrs_gr_sub <- dmrs_gr
#remove intercept
dmrs_gr_sub[[1]]<- NULL
dmrs_gr_sub <- lapply(dmrs_gr_sub, function(x){
    x <- x[abs(x$diff.Methy)>0.1,]
    x
})
lapply(dmrs_gr_sub, length)

saveRDS(dmrs_gr_sub, file.path(analysis.dir ,"dmrs_gr_sub_MethDiff.rds"))


#Annotate and plot dmrs 
#Load TxDb file
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

dmrs_gr<- list()
dmrs_anno<- list()
dmrs_anno_df<- list()
#plot files and annotate
for(i in names(dmrs_gr_sub)){
#annotate dmrs
dmrs_anno[[i]] <- annotatePeak(peak= dmrs_gr_sub[[i]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
dir.create(file.path(analysis.dir,i))
#visualize dmrs annotation
dir.create(file.path(analysis.dir,i,  "anno"))
pdf(file.path(analysis.dir,i,  "anno", paste0(i, "_AnnoDist_Pie.pdf")))
print(plotAnnoPie(dmrs_anno[[i]] ))
dev.off()

pdf(file.path(analysis.dir,i,"anno", paste0(i, "_AnnoDist_Upset.pdf")))
print(upsetplot(dmrs_anno[[i]] , vennpie=TRUE))
dev.off()

pdf(file.path(analysis.dir,i,"anno", paste0(i, "_AnnoDist_Upset_noVenn.pdf")))
print(upsetplot(dmrs_anno[[i]] , vennpie=FALSE))
dev.off()

pdf(file.path(analysis.dir, i, "anno",paste0(i, "_DistanceToTSS.pdf")))
print(plotDistToTSS(dmrs_anno[[i]], 
              title="Distribution of DMR relative to TSS"))
dev.off()  

#extract annotated granges
dmrs_anno_df[[i]] <- as.data.frame(dmrs_anno[[i]])
mcols(dmrs_gr_sub[[i]]) <- dmrs_anno_df[[i]][6:ncol(dmrs_anno_df[[i]])]
}

#assign direction
dmrs_gr_sub <- lapply(dmrs_gr_sub, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})

saveRDS(dmrs_gr_sub, file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno.rds"))

#get numbers of dmrs
direction_all_numbers <- lapply(dmrs_gr_sub, function(x){table(x$direction)})

#create "new tracks"
#make granges
for(i in names(dmrs_gr_sub)){
export.bed(dmrs_gr_sub[[i]],file.path(analysis.dir, "tracks", paste0("dmrs_sub_MethDiff_",i, ".bed")))
}

#reduce dmrs
dmrs_gr_red <- GRangesList(dmrs_gr_sub[[1]], dmrs_gr_sub[[2]], dmrs_gr_sub[[3]], dmrs_gr_sub[[4]], dmrs_gr_sub[[5]], dmrs_gr_sub[[6]])
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
mcols(dmrs_gr_red) <- dmrs_red_anno_df[, 6:ncol(dmrs_red_anno_df)]
saveRDS(dmrs_gr_red, file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))
