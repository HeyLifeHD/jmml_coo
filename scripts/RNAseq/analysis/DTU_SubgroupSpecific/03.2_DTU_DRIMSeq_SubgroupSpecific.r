#differential transcript usage analysis
#libraries
library(DRIMSeq)
library(ggpubr)
#folder
base.dir<- "/omics/groups/OE0219/internal/jmmlc_rnaseq/220815_rnaseq_DTU_DRIMSeq_Subgroup_analysis"
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
#relevel consensuscluster 3
sample_anno$consensusCluster3 <- factor(sample_anno$consensusCluster3 ,  levels=c("LM", "IM", "HM"))
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
design_full <- model.matrix(~consensusCluster3, data=DRIMSeq::samples(d))
colnames(design_full)

#run analysis
set.seed(1)
d <- dmPrecision(d, design=design_full, BPPARAM = BiocParallel::MulticoreParam(5))
d <- dmFit(d, design=design_full, BPPARAM = BiocParallel::MulticoreParam(5))
d_IM <- dmTest(d, coef="consensusCluster3IM", BPPARAM = BiocParallel::MulticoreParam(5))
d_HM <- dmTest(d, coef="consensusCluster3HM", BPPARAM = BiocParallel::MulticoreParam(5))

saveRDS(d_IM, file.path(results.dir, "DRIMSeq_Object_IM.rds"))
d_IM <- readRDS(file.path(results.dir, "DRIMSeq_Object_IM.rds"))
saveRDS(d_HM, file.path(results.dir, "DRIMSeq_Object_HM.rds"))
d_HM <- readRDS(file.path(results.dir, "DRIMSeq_Object_HM.rds"))

d_list <- list(IM=d_IM, HM=d_HM)
#save proportions and norm counts
for ( i in names(d_list)){
prop <- DRIMSeq::proportions(d_list[[i]])
dir.create(file.path(results.dir,i))
saveRDS(prop, file.path(results.dir, i, "proportions.rds"))
coun <- DRIMSeq::counts(d_list[[i]])
saveRDS(coun, file.path(results.dir,i,  "counts.rds"))
}

padj <- list()
for ( i in names(d_list)){
    #extract results
    res <- DRIMSeq::results(d_list[[i]])
    res.txp <- DRIMSeq::results(d_list[[i]], level="feature")
    print(i)
    print(nrow(res[which(res$adj_pvalue <0.1),]))
    print(nrow(res.txp[which(res.txp$adj_pvalue <0.1),]))

    #plot precision
    ggp <- plotPrecision(d_list[[i]])
    dir.create(file.path(PreDE.dir, i))
    pdf(file.path(PreDE.dir, i, "DrimSeq_Precision.pdf"))
    print(ggp + geom_point(size = 4))
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
    for (n in 1:2) tx2gene[,n] <- strp(tx2gene[,n])

    library(stageR)
    stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=FALSE, tx2gene=tx2gene)
    stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.1)
    drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                    onlySignificantGenes=TRUE)

    stageR_genes <- getSignificantGenes(stageRObj)
    stageR_transcripts <- getSignificantTx(stageRObj)
    padj[[i]] <- getAdjustedPValues(stageRObj, order = TRUE,
        onlySignificantGenes = FALSE)
    print(head(padj[[i]]))
    dir.create(file.path(PostDE.dir,i))
    write.table(padj[[i]], file.path(PostDE.dir,i,  "StageR_adj_results.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
    saveRDS(padj[[i]], file.path(PostDE.dir,i,  "StageR_adj_results.rds"))

    #plot sign genes
    idx <- unique(padj[[i]][which(padj[[i]]$gene < 0.1),]$geneID)
    print(length(idx))
    dir.create(file.path(PostDE.dir,i,  "sign_DTU_proportions"))
    for(j in idx){
        print(padj[[i]][padj[[i]]$geneID ==j,])
        pdf(file.path(PostDE.dir, i, "sign_DTU_proportions",paste0("proportions_", j,"_Epigenotype.pdf")))
        print(plotProportions(d_list[[i]], j, "Epigenotype"))
        dev.off()
        pdf(file.path(PostDE.dir, i, "sign_DTU_proportions",paste0("proportions_", j,"_Epigenotype_Ribbon.pdf")))
        print(plotProportions(d_list[[i]], j, "Epigenotype",plot_type = "ribbonplot"))
        dev.off()
        pdf(file.path(PostDE.dir,i, "sign_DTU_proportions",paste0("proportions_", j,"_consensusCluster3.pdf")))
        print(plotProportions(d_list[[i]],  j, "consensusCluster3"))
        dev.off()
        pdf(file.path(PostDE.dir,i, "sign_DTU_proportions",paste0("proportions_", j,"_consensusCluster3_Ribbon.pdf")))
        print(plotProportions(d_list[[i]],  j, "consensusCluster3",plot_type = "ribbonplot"))
        dev.off()
    }
}
paste0(unique(padj[["HM"]][which(padj[[i]]$gene < 0.1),]$geneID), collapse=", ")
paste0(unique(padj[["IM"]][which(padj[[i]]$gene < 0.1),]$geneID), collapse=", ")

#look at overlap
overl <- list(HM_DTUs= unique(padj[["HM"]][which(padj[[i]]$gene < 0.1),]$geneID), IM_DTUs=padj[["IM"]][which(padj[[i]]$gene < 0.1),]$geneID)
pdf(file.path(PostDE.dir, "Venn_DTUs.pdf"))
ggVennDiagram::ggVennDiagram(overl, label_alpha = 0, category.names=names(overl))+
    rremove("legend") + 
    theme(text = element_text(size = 8))
dev.off()



#integrate with deseq2 gene expression results
#load deseq 2 results
DEG_results_list<- readRDS("/omics/groups/OE0219/internal/jmmlc_rnaseq/220819_rnaseq_knownDEG_subgroupSpecific/results/PostDE/DEG_results_group_list.rds")
names(DEG_results_list)<- c("IM", "HM")

pal <- c("not_significant"="darkgray", "DGE"="darkred", "DTU"="darkblue", "DTE"="darkgreen")
for(i in names(DEG_results_list)){
  dres <- DEG_results_list[[i]]
  padj_gene <- padj[[i]]
  padj_gene$gene_id <- padj_gene$geneID
  padj_gene<- unique(padj_gene[,c("geneID","gene", "gene_id")])
  dres_comb <- dplyr::left_join(dres,padj_gene, by="gene_id" )
  dres_comb$adj_pvalue <- ifelse(is.na(dres_comb$gene), 1, dres_comb$gene)
  dres_comb$col <- "not_significant"
  dres_comb$col <- ifelse(dres_comb$padj < 0.05, "DGE", dres_comb$col)
  dres_comb$col <- ifelse(dres_comb$adj_pvalue < 0.1, "DTU", dres_comb$col)
  dres_comb$col <- ifelse(dres_comb$adj_pvalue < 0.1 & dres_comb$padj < 0.05, "DTE", dres_comb$col)
  dres_comb$dge_pval <- -log10(dres_comb$padj)
  dres_comb$dte_pval <- -log10(dres_comb$adj_pvalue)
  temp <- c(head(dres_comb[dres_comb$adj_pvalue < 0.1,]$gene_id,Inf), 
    head(dres_comb[dres_comb$padj < 0.05,]$gene_id,15))
  dres_comb$label <- ifelse(dres_comb$gene_id %in% temp, dres_comb$gene_id, NA)
  pdf(file.path(PostDE.dir, i,  "Drim_seq_transcript_vs_geneUsage.pdf"))
  print(ggscatter(dres_comb, x="dge_pval", y="dte_pval", 
    col="col", palette=pal, 
    label= "label", repel=TRUE,
    xlab="Differential gene expression\n-log10(adjusted P value)", 
    ylab="Differential transcript usage\n-log10(adjusted P value)"))
  dev.off()
  pdf(file.path(PostDE.dir,i,   "Drim_seq_transcript_vs_geneUsage_noLabel.pdf"))
  print(ggscatter(dres_comb, x="dge_pval", y="dte_pval", 
    col="col", palette=pal,
    #label= "label", repel=TRUE,
    xlab="Differential gene expression\n-log10(adjusted P value)", 
    ylab="Differential transcript usage\n-log10(adjusted P value)"))
  dev.off()
}

#look at hm vs non-hm dmrs overlapping with dtudwer
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/220512_DMR_hierachy_HSC_TumorNormal/"
dmrs_final<- readRDS(file.path(analysis.dir, "sig_dmrs_5inHalf_sub_anno.rds"))
dmrs_final_sub <- dmrs_final[c("LM_vs_HM", "LM_vs_IM")]
names(dmrs_final_sub)<- c("HM", "IM")
dmrs_final_sub <- lapply(dmrs_final_sub, function(x) {
  x$diff.Methy <- -x$diff.Methy 
  x})
for ( i in names(DEG_results_list)){
  padj_gene <- padj[[i]]
  padj_gene$gene_id <- padj_gene$geneID
  padj_gene<- unique(padj_gene[,c("geneID","gene", "gene_id")])
  res_sub <-padj_gene[which(padj_gene$gene_id %in% dmrs_final_sub[[i]]$SYMBOL),]
  dmrs_final$EpigenotypeHM[dmrs_final$EpigenotypeHM$SYMBOL %in% padj_gene[which(padj_gene$gene < 0.1),]$gene_id,]
  write.table(res_sub, file.path(PostDE.dir,i,  "StageR_adj_results_DMRoverlap.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
  #plot overlap
  overl <- list(DTUs_genes= as.character(padj_gene[which(padj_gene$gene < 0.1),]$gene_id), DMR_genes=unique(dmrs_final_sub[[i]]$SYMBOL))
  pdf(file.path(PostDE.dir, i,  "Venn_DMR_genes_DTU.pdf"))
  print(ggVennDiagram::ggVennDiagram(overl, label_alpha = 0, category.names=names(overl))+
      rremove("legend") + 
      theme(text = element_text(size = 8)))
  dev.off()
}