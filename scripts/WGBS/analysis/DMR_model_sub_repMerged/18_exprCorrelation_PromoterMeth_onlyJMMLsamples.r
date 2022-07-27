
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
library(dplyr)

#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))
#change annotation
pheno <- pData(bsseq_all)
pheno$Tissue<- c( rep("tumor", 15))
pheno$Donor <- c(as.character(pheno$Patient[1:15])) 
pheno$Epigenotype <- as.character(pheno$Epigenotype)
pheno[pheno$Patient %in% c("I217", "D217"),]$Epigenotype <- "IM"
pheno$Epigenotype <- as.factor(pheno$Epigenotype)
pheno$Celltype <- c("MPP","MPP","MPP","MPP" ,"MPP","HSC", "HSC", "HSC","LMPP","LMPP","LMPP","CD45RACD90","CD45RACD90","CD45RACD90","CD45RACD90")
pheno$Sample_Type <- "tumor"
pheno$Tumor<- pheno$tumor
pData(bsseq_all) <- pheno

#extract methylation levels for each promoter per gene
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ENS<-EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
genes <- genes(ENS)
prom <- promoters(genes, upstream=1500, downstream=1500)
prom <- prom[seqnames(prom) %in%  c(paste0("chr", 1:21),"chrY", "chrX")]
meth_prom <- bsseq::getMeth(bsseq_all, regions= prom, type = "raw", what=c("perRegion"))
mcols(prom)<- cbind(mcols(prom), meth_prom)
prom$promID <- paste0("Promoter_", 1:length(prom))
prom$ENSEMBL <- prom$gene_id
prom_list <- list(promoter_meth= prom)


#load expression data
expr <- readRDS("/omics/groups/OE0219/internal/jmmlc_pbat/data/expression/201007_DEG/Statistics/HSC.PseudoBulk.rds")

#loop over each dmr subset and run
corr_res_sub <- list()
expr_final<-list()
prom_final<- list()
for(i in names(prom_list)){
    #select expressed genes of called dmrs and get same ordering of rows
    expr_sub <- expr[rownames(expr) %in%  prom_list[[i]]$ENSEMBL,]
    prom_sub <- as.data.frame(mcols(prom_list[[i]]))
    prom_sub <- prom_sub[prom_sub$ENSEMBL %in% rownames(expr_sub),]
    expr_sub$ENSEMBL <- rownames(expr_sub)
    merged <- inner_join(prom_sub, expr_sub, by="ENSEMBL")
    expr_sub <- merged[,colnames(expr_sub)]
    prom_sub <- merged[,colnames(prom_sub)]
    rownames(prom_sub)<- paste0(merged$symbol, "_",merged$ENSEMBL)
    rownames(expr_sub) <- rownames(prom_sub)
    #get average promoter methylation per patient and cordblood hsc samples
    pheno <-pData(bsseq_all)
    group <- unique(pheno$Donor)
    group <- lapply(group, function(x){
    x <-rownames(pheno[pheno$Donor %in% x,])
    x
    })
    names(group)<-  unique(pheno$Donor)
    prom_sub_mean <- lapply(group, function(x){
        if(length(x) >1){
        x <-  rowMeans(prom_sub[,x], na.rm=TRUE)
        }else{
        x <- prom_sub[,x, drop=TRUE]
        }
    x
    })
    prom_sub_mean <- do.call("cbind", prom_sub_mean)
    #get sample annotation the same
    expr_sub <- expr_sub[,c("D129", "D213", "I217", "D217", "D360", "D123", "D124", "D117")]
    colnames(expr_sub)<- c("D129", "D213", "I217", "D217", "D360", "D123", "D124", "D117")
    expr_sub <- expr_sub[, colnames(prom_sub_mean)]
    dim(expr_sub)
    dim(prom_sub_mean)
    table(colnames(prom_sub_mean)==colnames(expr_sub))
    table(rownames(prom_sub_mean)==rownames(expr_sub))

    #subset only complete observations
    idx <- complete.cases(prom_sub_mean)
    prom_sub_mean <- prom_sub_mean[idx,]
    expr_sub <- expr_sub[idx,]
    merged_sub <- merged[idx,]
    #remove cordblood
    expr_final[[i]]<- expr_sub
    prom_final[[i]]<- prom_sub_mean
    expr_sub <- expr_sub
    prom_sub_mean <- prom_sub_mean
 

    #perform correlation test
    corr <- list()
        for (j in rownames(expr_sub)){
            corr[[j]]<- cor.test(as.numeric(expr_sub[j,]), as.numeric(prom_sub_mean[j,]), method="pearson", alternative = "two.sided", exact=F)
        }
    #extract results
    corr_res <- data.frame(p.value=unlist(lapply(corr, function(x)x$p.value)), estimate=unlist(lapply(corr, function(x)x$estimate)), symbol=merged_sub$symbol, peakID= merged_sub$ENSEMBL)
    #corr_res$symbol <- prom[rownames(corr_res),]$SYMBOL
    #subset meaningful results and adjust for multiple testing
    corr_res_sub[[i]] <- corr_res[complete.cases(corr_res),]
    corr_res_sub[[i]]$p.adjust <- p.adjust(corr_res_sub[[i]]$p.value, "hochberg")
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$p.value, decreasing=FALSE),]
    print(head(corr_res_sub[[i]], 20))
    
    #plot results as histogramm
    corr_res_sub[[i]]$sign <- ifelse(corr_res_sub[[i]]$p.value<0.01, "P value < 0.01", "P value > 0.01")
    dir.create(file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth"), recursive=TRUE)
    write.table(corr_res_sub[[i]], file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "correlation_res.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    pdf(file=file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "correlation_res_histogram.pdf"), height=3.5, width=5)
    print(ggpubr::gghistogram(corr_res_sub[[i]], x="estimate", bins=30, add="mean",
        xlab="Pearson correlation estimate\nAverage DMR methylation vs closest gene expression",ylab="Number of observations",
        fill="sign", palette=c("red", "black")))
    dev.off()

     #plot correlation of significant genes
    corr_res_sub_sig <- corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value<0.01), ]
    #plot correlation of top 50
    #corr_res_sub_sig <- head(corr_res_sub[[i]],50)

    anno_col <- c(D117 = "#0058b4", D129 = "#2188c9", 
                D217 = "#fbbb25", I217 = "#fca349", 
                D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
                cordblood ="#ababab") 

    dir.create(file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth","scatter_plots"), recursive=TRUE)
    for(j in rownames(corr_res_sub_sig)){
        temp <- data.frame(Expression=as.numeric(expr_sub[j,]), Methylation=as.numeric(prom_sub_mean[j,]), Donor=colnames(prom_sub_mean))
        pdf(file=file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth","scatter_plots", paste0("scatter_plot_", j, ".pdf")), height=5, width=5)
            print(ggpubr::ggscatter(temp, x="Expression",y="Methylation",
                xlab="Normalized gene epression",ylab="Average DMR methylation",color="Donor", palette=anno_col,
                shape = 16, size = 2, title=paste0(j, "\n", "estimate = ", round(corr_res_sub_sig[j,]$estimate , 3),"; p.adjust = ", round(corr_res_sub_sig[j,]$p.adjust , 5), ";\np.value = ", round(corr_res_sub_sig[j,]$p.value , 5)   ),
                add = "reg.line",  # Add regressin line
                conf.int = T, # Add confidence interval
                add.params = list(color = "grey", fill = "lightgray"),
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")))
        dev.off()
    }
    print(paste0(i, " done"))
}

saveRDS(corr_res_sub,file.path(analysis.dir,"corr_res_sub_PromoterMeth.rds"))
saveRDS(prom_final,file.path(analysis.dir,"corr_res_sub_methylation_values_PromoterMeth.rds"))
saveRDS(expr_final,file.path(analysis.dir,"corr_res_sub_expression_values_PromoterMeth.rds"))

corr_res_sub <- readRDS(file.path(analysis.dir,"corr_res_sub_PromoterMeth.rds"))
prom_final <- readRDS(prom_final,file.path(analysis.dir,"corr_res_sub_methylation_values_PromoterMeth.rds"))
expr_final <- readRDS(expr_final,file.path(analysis.dir,"corr_res_sub_expression_values_PromoterMeth.rds"))

#Plot histogramm of estimates
for(i in names(corr_res_sub)){
    corr_res_sub[[i]]$id <- rownames(corr_res_sub[[i]])
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$estimate, decreasing=TRUE),]

    pdf(file=file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "waterfall_estimate_all.pdf"), height=3.5, width=5)
    print(ggbarplot(corr_res_sub[[i]], x="id", y="estimate",ylab="correlation estimate\n[r]", fill="gray", color="gray") +rremove("x.axis")+rremove("x.grid")+rremove("xlab")+rremove("x.text")+rremove("x.ticks"))
    dev.off()
}
for(i in names(corr_res_sub)){
    corr_res_sub[[i]]$id <- rownames(corr_res_sub[[i]])
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$estimate, decreasing=TRUE),]
    corr_res_sub_sig <- corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value<0.01), ]

    pdf(file=file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "waterfall_estimate_pvalue0.01.pdf"), height=3.5, width=5)
    print(ggbarplot(corr_res_sub_sig, x="id", y="estimate",ylab="correlation estimate\n[r]", fill="gray", color="gray") +rremove("x.axis")+rremove("x.grid")+rremove("xlab")+rremove("x.text")+rremove("x.ticks"))
    dev.off()
}

#make heatmap of hits
anno_colors = list(
  Patient = c(D117 = "#0058b4", D129 = "#2188c9", 
            D217 = "#fbbb25", I217 = "#fca349", 
            D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"), 
  Epigenotype = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4")
  )
anno_heat <- data.frame(Patient=c("D117", "D129", "D213", "D217", "I217", "D124", "D123", "D360"),
    Epigenotype=c("LM", "LM", "HM", "IM", "IM", "HM","HM",  "HM"))
rownames(anno_heat) <-  anno_heat$Patient

for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
rna_sig <- expr_final[[i]][rownames(corr_res_sub_sig_ord),]
wgbs_sig <-  prom_final[[i]][rownames(corr_res_sub_sig_ord),]

nrow(wgbs_sig)==nrow(rna_sig)
#plot heatmap #estimate order
dir.create(file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps"))
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_RNA_signCorr0.01_estimateOrder_colClustered.pdf")
) 
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_WGBS_signCorr0.01_estimateOrder_colClustered.pdf")
) 

#plot heatmap #clustering order independent
dir.create(file.path(analysis.dir, "heatmaps"))
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_RNA_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
)  
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_WGBS_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
) 

#plot heatmap #clustering order based on rna
dir.create(file.path(analysis.dir, "heatmaps"))
p <- pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps","heatmap_RNA_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
)  
pheatmap::pheatmap(wgbs_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
) 
p <- pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
) 
pheatmap::pheatmap(rna_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
)  
}


#make heatmap of hits that overlap with degs
#load deg data
files <- list.files(path="/omics/groups/OE0219/internal/jmmlc_pbat/data/expression/200825_DEG/HM_vs_LM/", pattern=".rds", full.names=TRUE)
names(files)<- paste0("promoter_meth")
deg_sig <- lapply(files, function(x){
    x <-readRDS(x)
    x$SYMBOL <-x$gene 
    x <- x[which(x$p_val_adj < 0.05),]
    x
    }
)

for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
corr_res_sub_sig_ord <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% unique(deg_sig[[i]]$gene),]
rna_sig <- expr_final[[i]][rownames(corr_res_sub_sig_ord),]
wgbs_sig <-  prom_final[[i]][rownames(corr_res_sub_sig_ord),]

nrow(wgbs_sig)==nrow(rna_sig)
#plot heatmap #estimate order
dir.create(file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs"))
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_estimateOrder_colClustered.pdf")
) 
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_estimateOrder_colClustered.pdf")
) 

#plot heatmap #clustering order independent
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
)  
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
) 


#plot heatmap #clustering order based on rna
p <- pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs","heatmap_RNA_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
)  
pheatmap::pheatmap(wgbs_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
) 
#clustering based on wgbs
p <- pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
) 

pheatmap::pheatmap(rna_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
) 
#row clustering based on wgbs; columns based on methylation level
p <- pheatmap::pheatmap(wgbs_sig[,c(  "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colMeanMeth.pdf")
) 

pheatmap::pheatmap(rna_sig[p$tree_row$order,c(  "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colMeanMeth.pdf")
) 


#row clustering based on rna; columns based on methylation level
p <- pheatmap::pheatmap(rna_sig[,c(  "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredRNA_colMeanMeth.pdf")
) 

pheatmap::pheatmap(wgbs_sig[p$tree_row$order,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colMeanMeth.pdf")
) 
}

#seperate positive and negative estimate
for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
corr_res_sub_sig_ord <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% unique(deg_sig[[i]]$gene),]

rna_sig_neg <- expr_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate < 0,]),, drop=FALSE]
wgbs_sig_neg <-  prom_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate < 0,]),, drop=FALSE]

rna_sig_pos <- expr_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate > 0,]),, drop=FALSE]
wgbs_sig_pos <-  prom_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate > 0,]),, drop=FALSE]

#row clustering ordered; columns based on methylation level
#negative estimate
if(nrow(wgbs_sig_neg)> 1){
pheatmap::pheatmap(wgbs_sig_neg[,c(  "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="none", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowOrdEstimate_colMeanMethh_NegativeEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_neg[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="row", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_rowOrdEstimate_rowClusteredWGBS_colMeanMeth_NegativeEstimate.pdf")
) 
}
if(nrow(wgbs_sig_pos)>1){
#positive estimate
pheatmap::pheatmap(wgbs_sig_pos[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="none", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_rowOrdEstimate_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_pos[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="row", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA__rowOrdEstimate_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
)  
}
if(nrow(wgbs_sig_neg)> 1){
#row clustering based on wgbs; columns based on methylation level
#negative estimate
p <- pheatmap::pheatmap(wgbs_sig_neg[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE], 
    scale="none", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colMeanMethh_NegativeEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_neg[p$tree_row$order,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="row", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colMeanMeth_NegativeEstimate.pdf")
) 
}
if(nrow(wgbs_sig_pos)>1){
#positive estimate
p <- pheatmap::pheatmap(wgbs_sig_pos[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE], 
    scale="none", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_pos[p$tree_row$order,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="row", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
) 
}    
if(nrow(wgbs_sig_neg)> 1){
#row clustering based on rna; columns based on methylation level
#negative estimate
p <- pheatmap::pheatmap(rna_sig_neg[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE], 
    scale="row", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredRNA_colMeanMeth_NegativeEstimate.pdf")
) 

pheatmap::pheatmap(wgbs_sig_neg[p$tree_row$order,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="none", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colMeanMeth_NegativeEstimate.pdf")
) 
}
if(nrow(wgbs_sig_pos)>1){
#positive estimate
p <- pheatmap::pheatmap(rna_sig_pos[,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="row", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredRNA_colMeanMeth_PositiveEstimate.pdf")
) 

pheatmap::pheatmap(wgbs_sig_pos[p$tree_row$order,c(   "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],
    scale="none", show_colnames=F,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "PromoterMeth", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colMeanMeth_PositiveEstimate.pdf")
) 
}
}

#save results
deg_up <- lapply(deg_sig, function(x){
    x <- x[x$avg_logFC>0,]
    x
})
deg_down <- lapply(deg_sig, function(x){
    x <- x[x$avg_logFC<0,]
    x
})
deg_overlap_up <- list()
deg_overlap_down <- list()
for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
deg_overlap_up[[i]] <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% unique(deg_up[[i]]$gene),]
deg_overlap_down[[i]] <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% unique(deg_down[[i]]$gene),]
}
deg_overlap <- list(upregulated=deg_overlap_up, downregulated=deg_overlap_down)
saveRDS(deg_overlap,file.path(analysis.dir,"corr_res_sign0.05_overlapDEGs_up_downStrat_PromoterMeth.rds") )

paste0(unique(deg_overlap$upregulated$promoter_meth$symbol), collapse=", ")
paste0(unique(deg_overlap$downregulated$promoter_meth$symbol), collapse=", ")
paste0(unique(c(as.character(deg_overlap$downregulated$EpigenotypeHM$symbol), as.character(deg_overlap$upregulated$EpigenotypeHM$symbol))), collapse=", ")