#differential transcript usage analysis
#libraries
library(DRIMSeq)
library(ggpubr)
#folder
base.dir<- "/omics/groups/OE0219/internal/jmmlc_rnaseq/220805_rnaseq_DTU_analysis"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)

#load data
counts <- read.table("/omics/groups/OE0219/internal/jmmlc_rnaseq/220802_rnaseq_processing_known_hg19_starRSEM/star_salmon/salmon.merged.transcript_counts.tsv", header=TRUE)
rownames(counts) <- counts$tx
counts <- counts[,!colnames(counts) %in%c("tx", "gene_id")]
colnames(counts)<- sapply(strsplit(colnames(counts), "_", fixed=TRUE), "[",2)
#create annotation sheet
sample_anno <- data.table::fread("/omics/groups/OE0219/internal/jmmlc_rnaseq/data/2020-09-29_Patient_Characteristics_Full_Meta.csv")
sample_anno <- sample_anno[sample_anno$consensusCluster3 %in% c("LM", "IM", "HM"),]
sample_anno$Epigenotype <- ifelse(sample_anno$consensusCluster3 =="HM", "HM", "non_HM")
sample_anno <- as.data.frame(sample_anno[sample_anno$EWOG_ID %in% colnames(counts),])
rownames(sample_anno)<- sample_anno$EWOG_ID 

counts <- counts[,colnames(counts)%in% sample_anno$EWOG_ID ]
cts<- counts[,sample_anno$EWOG_ID ]
cts <- cts[rowSums(cts) > 0,]

#get tx to gene mapping
txdf <- read.table("/omics/groups/OE0219/internal/jmmlc_rnaseq/220802_rnaseq_processing_known_hg19_starRSEM/star_salmon/salmon_tx2gene.tsv", header=FALSE)
txdf <- txdf[,1:2]
colnames(txdf)<- c( "TXNAME","GENEID")
txdf$ntx <- table(txdf$GENEID)[txdf$GENEID]
txdf <- txdf[,c("GENEID", "TXNAME", "ntx")]
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

#build input
counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)

#create drim seq object
sample_anno$sample_id <- sample_anno$EWOG_ID
d <- dmDSdata(counts=counts, samples=sample_anno)
range(colSums(cts)/1e6)

#filter transcripts
n <- 20
n.small <- 6
d_sub <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)

table(table(counts(d)$gene_id))


#create design formula
design_full <- model.matrix(~Epigenotype, data=DRIMSeq::samples(d))
colnames(design_full)

#run analysis
set.seed(1)

d <- dmPrecision(d, design=design_full)
d <- dmFit(d, design=design_full)
d <- dmTest(d, coef="Epigenotypenon_HM")

saveRDS(d, file.path(results.dir, "DRIMSeq_Object.rds"))
d <- readRDS(file.path(results.dir, "DRIMSeq_Object.rds"))

#extract results
res <- DRIMSeq::results(d)
head(res)
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
nrow(res[which(res$adj_pvalue <0.05),])
write.table(res, file.path(PostDE.dir, "DRIMSeq_results_geneLevel.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
nrow(res.txp[which(res.txp$adj_pvalue <0.05),])

write.table(res.txp, file.path(PostDE.dir, "DRIMSeq_results_transcriptLevel.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)

#plot precision
ggp <- plotPrecision(d)
pdf(file.path(PostDE.dir, "DrimSeq_Precision.pdf"))
ggp + geom_point(size = 4)
dev.off()

#fill up na
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

#stage pvalues
nrow(res)
nrow(res.txp)
pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)
pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)
tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

library(stageR)
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.1)
drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)

head(drim.padj)
nrow(drim.padj[which(drim.padj$gene <0.05),])
write.table(drim.padj, file.path(PostDE.dir, "DRIMSeq_results_postStageR_transcriptLevel.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
#post hoc filtering based on sd
res.txp.filt <- DRIMSeq::results(d, level="feature")
smallProportionSD <- function(d, filter=0.1) {
  cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  propSD <- sqrt(rowVars(props))
  propSD < filter
}
filt <- smallProportionSD(d)
table(filt)
res.txp.filt$pvalue[filt] <- 1 
res.txp.filt$adj_pvalue[filt] <- 1

#plot sign targets
idx <- which(res.txp.filt$adj_pvalue < 0.05)
length(idx)
dir.create(file.path(PostDE.dir, "sign_DTU_proportions"))
for(i in 1:length(idx)){
    print(res.txp.filt[idx[i],])
    pdf(file.path(PostDE.dir, "sign_DTU_proportions",paste0("proportions_", res.txp.filt$gene_id[idx[i]],"_Epigenotype.pdf")))
    print(plotProportions(d, res.txp.filt$gene_id[idx[i]], "Epigenotype"))
    dev.off()
    pdf(file.path(PostDE.dir, "sign_DTU_proportions",paste0("proportions_", res.txp.filt$gene_id[idx[i]],"_consensusCluster3.pdf")))
    print(plotProportions(d, res.txp.filt$gene_id[idx[i]], "consensusCluster3"))
    dev.off()
}

#integrate with deseq2 gene expression results
#load deseq 2 results
DEG_results_list<- readRDS("/omics/groups/OE0219/internal/jmmlc_rnaseq/220805_rnaseq_knownDEG/results/PostDE/DEG_results_group_list.rds")
dres <- DEG_results_list$HM_vs_non_HM

dres_comb <- dplyr::left_join(dres,res, by="gene_id" )
dres_comb$adj_pvalue <- ifelse(is.na(dres_comb$adj_pvalue), 1, dres_comb$adj_pvalue)
dres_comb$col <- "not_significant"
dres_comb$col <- ifelse(dres_comb$padj < 0.05, "DGE", dres_comb$col)
dres_comb$col <- ifelse(dres_comb$adj_pvalue < 0.05, "DTU", dres_comb$col)
dres_comb$col <- ifelse(dres_comb$adj_pvalue < 0.05 & dres_comb$padj < 0.05, "DTE", dres_comb$col)
dres_comb$dge_pval <- -log10(dres_comb$padj)
dres_comb$dte_pval <- -log10(dres_comb$adj_pvalue)

temp <- c(head(dres_comb[dres_comb$adj_pvalue < 0.05,]$gene_id,Inf), 
  head(dres_comb[dres_comb$padj < 0.05,]$gene_id,15))
dres_comb$label <- ifelse(dres_comb$gene_id %in% temp, dres_comb$gene_id, NA)

pdf(file.path(PostDE.dir,  "Drim_seq_transcript_vs_geneUsage.pdf"))
ggscatter(dres_comb, x="dge_pval", y="dte_pval", 
  col="col", palette=c("darkred", "darkblue", "darkgray"),
  label= "label", repel=TRUE,
  xlab="Differential gene expression\n-log10(adjusted P value)", 
  ylab="Differential transcript usage\n-log10(adjusted P value)")
dev.off()

all(res$gene %in% rownames(dres))
dres <- dres[res$gene_id,]

col <- rep(8, nrow(dres))
dte.genes <-  rownames(dres[which(dres$padj < 0.05),])
dte.genes <-dte.genes[dte.genes %in% res[which(res$adj_pvalue < 0.05),]$gene_id]
col[rownames(dres) %in% rownames(dres[which(dres$padj < 0.05),])] <- 1
col[rownames(dres) %in% dte.genes] <- 2
col[rownames(dres) %in% rownames(res[which(res$adj_pvalue < 0.05),])] <- 3
pdf(file.path(PostDE.dir,  "Drim_seq_transcript_vs_geneUsage.pdf"))

bigpar()
# here cap the smallest DESeq2 adj p-value
cap.padj <- pmin(-log10(dres$padj), 100)
# this vector only used for plotting
res[is.na(res$adj_pvalue),]$adj_pvalue <- 1
jitter.padj <- -log10(res$adj_pvalue+ 1e-20)
jp.idx <- jitter.padj == 20
jitter.padj[jp.idx] <- rnorm(sum(jp.idx),20,.25)

pdf(file.path(PostDE.dir,  "Drim_seq_transcript_vs_geneUsage.pdf"))
plot(cap.padj, jitter.padj, col=col,
     xlab="Gene expression",
     ylab="Transcript usage")
legend("topright",
       c("DGE","DTE","DTU","null"),
        pch=20, bty="n")
dev.off()



#look at hm vs non-hm dmrs overlapping with dtudwer
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
paste0(res[which(res$adj_pvalue < 0.05),]$gene_id, collapse=", ")

res_sub <-res[which(res$gene_id %in% dmrs_final$EpigenotypeHM$SYMBOL),]
paste0(res_sub[which(res_sub$adj_pvalue < 0.05),]$gene_id, collapse=", ")
dmrs_final$EpigenotypeHM[dmrs_final$EpigenotypeHM$SYMBOL %in% res_sub[which(res_sub$adj_pvalue < 0.05),]$gene_id,]


res.txp.filt_sub <-res.txp.filt [which(res.txp.filt $gene_id %in% dmrs_final$EpigenotypeHM$SYMBOL),]
paste0(res.txp.filt_sub [which(res.txp.filt_sub $adj_pvalue < 0.05),]$gene_id, collapse=", ")
dmrs_final$EpigenotypeHM[dmrs_final$EpigenotypeHM$SYMBOL %in% rres.txp.filt_sub [which(res.txp.filt_sub $adj_pvalue < 0.05),]$gene_id,]

#continue with dexseq analysis
library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + Epigenotype:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

#run desxeq
system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})
saveRDS(dxd, file.path(results.dir, "DEXSeq_Object.rds"))
dxd <- readRDS(file.path(results.dir, "DEXSeq_Object.rds"))

#extract results
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

#stageR following DEXseq
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
pConfirmation[is.na(pConfirmation),]<- 1
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

#load deseq 2 results
DEG_results_list<- readRDS("/omics/groups/OE0219/internal/jmmlc_rnaseq/220805_rnaseq_knownDEG/results/PostDE/DEG_results_group_list.rds")
dres <- DEG_results_list$HM_vs_non_HM

all(dxr.g$gene %in% rownames(dres))
dres <- dres[dxr.g$gene,]

col <- rep(8, nrow(dres))
dte.genes <-  rownames(dres[which(dres$padj < 0.05),])
dte.genes <-dte.genes[dte.genes %in% rownames(dxr.g[which(dxr.g$qval < 0.05),])]
col[rownames(dres) %in% rownames(dres[which(dres$padj < 0.05),])] <- 1
col[rownames(dres) %in% dte.genes] <- 2
col[rownames(dres) %in% rownames(dxr.g[which(dxr.g$qval < 0.05),])] <- 3

bigpar()
# here cap the smallest DESeq2 adj p-value
cap.padj <- pmin(-log10(dres$padj), 100)
# this vector only used for plotting
jitter.padj <- -log10(dxr.g$qval + 1e-20)
jp.idx <- jitter.padj == 20
jitter.padj[jp.idx] <- rnorm(sum(jp.idx),20,.25)

pdf(file.path(PostDE.dir,  "transcript_vs_geneUsage.pdf"))
plot(cap.padj, jitter.padj, col=col,
     xlab="Gene expression",
     ylab="Transcript usage")
legend("topright",
       c("DGE","DTE","DTU","null"),
       col=c(1:3,8), pch=20, bty="n")
dev.off()



