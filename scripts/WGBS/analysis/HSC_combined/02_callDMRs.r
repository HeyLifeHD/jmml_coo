##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##call hierachy developmental dmrs

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
library(BiocParallel)
#directories
odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved.rds"))

#define groups for dmr calling
pheno <- pData(bsseq_all)
pheno$group <- paste0(pheno$Tissue_Source, "_", pheno$Celltype)
pheno_normal <- pheno[grep("normal",pheno$Sample_Type, invert=FALSE),]
groups <- pheno_normal$group
groups_normal <- unique(groups)
contrasts <-as.data.frame(combn(as.character(groups_normal), 2))
contrasts <- sapply(contrasts, as.character)
colnames(contrasts)<- paste0(contrasts[1,],"_vs_", contrasts[2,])
#run dmr tests
dml_test <- list()
mParam = MulticoreParam(workers=4, progressbar=TRUE)

#for(i in colnames(contrasts)){
for(i in colnames(contrasts)[2:length(colnames(contrasts))]){
#for(i in colnames(contrasts)[15]){
    gr1 <- rownames(pheno_normal[pheno_normal$group==contrasts[1,i],])
    gr2 <- rownames(pheno_normal[pheno_normal$group==contrasts[2,i],])

    dml_test[[i]] <- DMLtest(bsseq_all, group1 =gr1, group2 = gr2, smoothing = TRUE, BPPARAM=mParam)
    #dml_test[[i]] <- DMLtest(bsseq_all, group1 =gr1, group2 = gr2, smoothing = TRUE)

    #saveRDS(dml_test[[i]], file.path(analysis.dir ,"dml_fits", paste0(i, "_dml_fit.rds")))
    print(i)
}
saveRDS(dml_test, file.path(analysis.dir, "dml_test.rds"))

#Call Sig DML and DMR
sig_dmls <- list()
sig_dmrs <- list()
for(i in names(dml_test)){
    sig_dmls[[i]]<- callDML(dml_test[[i]] ,delta=0.1, p.threshold=0.05)
    sig_dmrs[[i]] <- callDMR(dml_test[[i]],delta=0.1, p.threshold=0.05,
             minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
    print(i)
    print(paste0(as.character(length( sig_dmrs[[i]])), " of ",as.character(length(dml_test))))
}
saveRDS(sig_dmls, file.path(analysis.dir, "sig_dmls.rds"))
saveRDS(sig_dmrs, file.path(analysis.dir, "sig_dmrs.rds"))

#assign direction
dmrs <- sig_dmrs
dmrs <- lapply(dmrs, function(x){
    x$direction = ifelse(x$diff>0, "hyper", "hypo")
    x
})
dmrs_gr <- list()
for(i in names(dmrs)){
#make granges
dmrs_gr[[i]]<- GRanges(
  seqnames = dmrs[[i]]$chr,
  ranges = IRanges(start =  dmrs[[i]]$start,
                   end =  dmrs[[i]]$end
  )
)
mcols(dmrs_gr[[i]])<- dmrs[[i]][,4:10]
}

saveRDS(dmrs_gr, file.path(analysis.dir, "dmrs_gr.rds"))
dmrs_gr <- readRDS(file.path(analysis.dir, "dmrs_gr.rds"))
#subset dmrs based on coverage
cov_drms <- list()
for(comp in names(dmrs_gr)){
  cov_drms[[comp]] <- bsseq::getCoverage(bsseq_all, regions= dmrs_gr[[comp]], type = "Cov", what=c("perRegionAverage"))
  print(comp)
}

#find out which dmrs have a coverage of at least 3 in at least two of the groups which are being compared
keepLoci <- list()
for(i in names(cov_drms)){
    gr1 <- rownames(pheno_normal[pheno_normal$group==contrasts[1,i],])
    gr2 <- rownames(pheno_normal[pheno_normal$group==contrasts[2,i],])

    keepLoci[[i]] <- which(rowSums(cov_drms[[i]][, gr1] >= 5) >= round(length( gr1)/2) &
                     rowSums(cov_drms[[i]][, gr2] >= 5) >= round(length( gr2)/2))
}
lapply(keepLoci, function(x)length(x))
temp <- as.vector(c(as.numeric(lapply(keepLoci, function(x)length(x)))))/as.vector(as.numeric(lapply(dmrs_gr, function(x)length(x))))
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
saveRDS(dmrs_gr_sub, file.path(analysis.dir, "sig_dmrs_5inHalf_sub.rds"))
dmrs_gr_sub <- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub.rds"))
#add methylation difference information
meth_dmrs<- list()
for(comp in names(dmrs_gr_sub)){
  meth_dmrs[[comp]] <- bsseq::getMeth(bsseq_all, regions= dmrs_gr_sub[[comp]], type = "raw", what=c("perRegion"))
  mcols(dmrs_gr_sub[[comp]])<- cbind( mcols(dmrs_gr_sub[[comp]]), meth_dmrs[[comp]]  )
  print(comp)
}


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

#pdf(file.path(analysis.dir,i,"anno", paste0(i, "_AnnoDist_Upset.pdf")))
#print(upsetplot(dmrs_anno[[i]] , vennpie=TRUE))
#dev.off()

#pdf(file.path(analysis.dir,i,"anno", paste0(i, "_AnnoDist_Upset_noVenn.pdf")))
#print(upsetplot(dmrs_anno[[i]] , vennpie=FALSE))
#dev.off()

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

saveRDS(dmrs_gr_sub, file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno.rds"))
dmrs_gr_sub <- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno.rds"))


#get numbers of dmrs
direction_all_numbers <- lapply(dmrs_gr_sub, function(x){table(x$direction)})
direction_all_numbers
#create "new tracks"
#make granges
dir.create(file.path(analysis.dir, "tracks"))
for(i in names(dmrs_gr_sub)){
temp <- dmrs_gr_sub[[i]]
mcols(temp)<- NULL
score(temp)<-  dmrs_gr_sub[[i]]$diff.Methy
export.bed(temp,file.path(analysis.dir, "tracks", paste0("sig_dmrs_5inHalf_sub_",i, ".bed")))
export.bedGraph(temp,file.path(analysis.dir, "tracks", paste0("sig_dmrs_5inHalf_sub_",i, ".bedGraph")))
}

#reduce dmrs
dmrs_gr_red <- GRangesList(dmrs_gr_sub[[1]], dmrs_gr_sub[[2]], dmrs_gr_sub[[3]], dmrs_gr_sub[[4]], dmrs_gr_sub[[5]],
                            dmrs_gr_sub[[6]], dmrs_gr_sub[[7]], dmrs_gr_sub[[8]], dmrs_gr_sub[[9]], dmrs_gr_sub[[10]],
                            dmrs_gr_sub[[11]], dmrs_gr_sub[[12]], dmrs_gr_sub[[13]], dmrs_gr_sub[[14]], dmrs_gr_sub[[15]])
#dmrs_gr_red <- GRangesList(dmrs_gr_sub[[1]])
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

#pdf(file.path(analysis.dir,"all","anno", paste0("all", "_AnnoDist_Upset.pdf")))
#upsetplot(dmrs_red_anno, vennpie=TRUE)
#dev.off()

#pdf(file.path(analysis.dir,"all","anno", paste0("all", "_AnnoDist_Upset_noVenn.pdf")))
#upsetplot(dmrs_red_anno , vennpie=FALSE)
#dev.off()

pdf(file.path(analysis.dir, "all", "anno",paste0("all", "_DistanceToTSS.pdf")))
plotDistToTSS(dmrs_red_anno, 
              title="Distribution of DMR relative to TSS")
dev.off()  

#extract annotated granges
dmrs_red_anno_df <- as.data.frame(dmrs_red_anno)
mcols(dmrs_gr_red) <- dmrs_red_anno_df[, 6:ncol(dmrs_red_anno_df)]
saveRDS(dmrs_gr_red, file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno_reduced.rds"))
