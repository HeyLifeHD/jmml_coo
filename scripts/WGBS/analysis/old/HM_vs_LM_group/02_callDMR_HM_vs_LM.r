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
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

dir.create(analysis.dir)
#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all_snpfil_sub_cov.rds"))

#Create subsetted bsseq
#bsseq <- list()
#dir.create(file.path(input.dir ,"bsseq", "bsseq_byChr"))
#for (chr in paste0("chr", c(1:22, "X", "Y"))) {
#bsseq[[chr]]<-  chrSelectBSseq(bsseq_all, seqnames =chr, order = TRUE)
#saveRDS(bsseq[[chr]],  paste0(input.dir ,"bsseq", "/bsseq_byChr/bsseq_", chr))
#print(chr)
#}
#saveRDS(bsseq, file.path(input.dir ,"bsseq","bsseq.rds"))
#bsseq<- readRDS(file.path(input.dir ,"bsseq", "bsseq.rds"))

#get pheno
pheno <- as.data.frame(pData(bsseq_all))

#do dml test for HM vs LM ciomparison
pheno <- pData(bsseq_all)
group1 <- rownames(pheno[pheno$Epigenotype=="HM",])
group2 <- rownames(pheno[pheno$Epigenotype=="LM",])

#dml_list <- list()
#for(chr in paste0("chr", c(1:22, "X", "Y"))) {
#    dml_list[[chr]] <- DMLtest(bsseq[[chr]], group1 = group1, group2 = group2, smoothing = TRUE)
#    print(paste0("chromosome ", chr, " processed"))
#}

mParam = MulticoreParam(workers=5, progressbar=TRUE)
dml_list <- DMLtest(bsseq_all, group1 = group1, group2 = group2, smoothing = TRUE,BPPARAM=mParam)

saveRDS(dml_list, file.path(analysis.dir, "dml_list.rds"))

#Call Sig DML and DMR
sig_dmls<- callDML(dml_list,delta=0.1, p.threshold=0.05)
sig_dmrs <- callDMR(dml_list,delta=0.1, p.threshold=0.05,
    minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
dim(sig_dmls)
dim(sig_dmrs)
saveRDS(sig_dmls, file.path(analysis.dir, "sig_dmls.rds"))
saveRDS(sig_dmrs, file.path(analysis.dir, "sig_dmrs.rds") )

#assign direction
sig_dmrs$direction = ifelse(sig_dmrs$diff.Methy>0, "hyper", "hypo")
table(sig_dmrs$direction)

#make granges
dmrs_gr<- GRanges(
  seqnames = sig_dmrs$chr,
  ranges = IRanges(start =  sig_dmrs$start,
                   end =  sig_dmrs$end
  )
)
mcols(dmrs_gr)<- sig_dmrs[,4:10]
saveRDS(dmrs_gr, file.path(analysis.dir, "dmrs_gr.rds") )
export.bed(dmrs_gr, file.path(analysis.dir, "dmrs_gr.bed") )

#subset dmrs based on coverage
cov_drms <- bsseq::getCoverage(bsseq_all, regions= dmrs_gr, type = "Cov", what=c("perRegionAverage"))

#find out which dmrs have a coverage of at least 3 in at least two of the groups which are being compared
keepLoci <- which(rowSums(cov_drms[, bsseq_all$Epigenotype == "HM"] >= 5) >= 4 &
                       rowSums(cov_drms[, bsseq_all$Epigenotype == "LM"] >= 5) >= 4)


#subset dmrs
dmrs_gr_sub  <- dmrs_gr[keepLoci,]

length(dmrs_gr_sub)
saveRDS(dmrs_gr_sub, file.path(analysis.dir, "sig_dmrs_5in4_sub.rds"))
export.bed(dmrs_gr_sub,  file.path(analysis.dir, "sig_dmrs_5in4_sub.bed"))

#Annotate and plot dmrs 
#Load TxDb file
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#annotate dmrs
dmrs_anno <- annotatePeak(peak= dmrs_gr_sub, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
dir.create(file.path(analysis.dir))
#visualize dmrs annotation
dir.create(file.path(analysis.dir,  "anno"))
pdf(file.path(analysis.dir, "anno", paste0( "AnnoDist_Pie.pdf")))
print(plotAnnoPie(dmrs_anno ))
dev.off()

pdf(file.path(analysis.dir,"anno", paste0( "AnnoDist_Upset.pdf")))
print(upsetplot(dmrs_anno , vennpie=TRUE))
dev.off()

pdf(file.path(analysis.dir,"anno", paste0( "AnnoDist_Upset_noVenn.pdf")))
print(upsetplot(dmrs_anno , vennpie=FALSE))
dev.off()


pdf(file.path(analysis.dir,  "anno",paste0( "DistanceToTSS.pdf")))
print(plotDistToTSS(dmrs_anno, 
              title="Distribution of DMR relative to TSS"))
dev.off()  

#extract annotated granges
dmrs_anno_df<- as.data.frame(dmrs_anno)
mcols(dmrs_gr_sub) <- dmrs_anno_df[,6:24]


#assign direction
dmrs_gr_sub$direction = ifelse(dmrs_gr_sub$diff.Methy>0, "hyper", "hypo")
table(dmrs_gr_sub$direction)

saveRDS(dmrs_gr_sub, file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

#subset unique dmrs
dmr_all_meth <- bsseq::getMeth(bsseq_all, regions=dmrs_gr_sub, what="perRegion", type="raw")

idx <- complete.cases(dmr_all_meth)
anno <- dmrs_all_gr[idx,]$SYMBOL
anno<-as.data.frame(anno)
plot <- dmr_all_meth[idx,]
pbat_col = list(
  patient = c(D117 = "#002e89", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
              D213 = "#DE8511", D360 = "#BD5D0C", D124 = "#9C3507", D123 = "#700000"), 
  Genotype = c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458"),
  Epigenotype = c(HM = "#b5241c", LM = "#0041a5"), 
  tumor = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#aeb9d4", tumor11 = "#d9c3ac")
  )
annovst <- as.data.frame(colData(bsseq_all))[, c("Epigenotype", "Patient", "tumor", "Genotype")] 
pheatmap(plot, annotation_col=annovst,scale="row",annotation_color=anno_colors , fontsize_row=5, show_rownames = FALSE, file=file.path(analysis.dir,"Heatmap_allDMRs.pdf"))
pheatmap(plot, annotation_col=annovst,scale="none",annotation_color=anno_colors,       fontsize_row=5, show_rownames = FALSE, file=file.path(analysis.dir,"Heatmap_allDMRs_noscale.pdf"))













