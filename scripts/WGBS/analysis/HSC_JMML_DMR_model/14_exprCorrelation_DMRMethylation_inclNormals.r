
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

#load my own data: Epigenotype vs HSC_cb
#Directories
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))
#add all
dmrs_final$all <- dmrs_red
mcols(dmrs_final$all)$direction <- "hypo"

#annotate dmrs
dmrs_anno<- list()
dmrs_anno_df<- list()
#plot files and annotate
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
for(i in names(dmrs_final)){
#annotate dmrs
dmrs_anno[[i]] <- annotatePeak(peak= dmrs_final[[i]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

#extract annotated granges
dmrs_anno_df[[i]] <- as.data.frame(dmrs_anno[[i]])
mcols(dmrs_final[[i]]) <- dmrs_anno_df[[i]]
}

dmrs_final <- lapply(dmrs_final, function(x){
    mcols(x)$peakID <- 1:length(x) 
    x
})

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

#load expression data
expr <- readRDS("/icgc/dkfzlsdf/analysis/C010/JMMLC/scRNA_Data/New/DEG/Statistics/HSC.PseudoBulk.rds")

#loop over each dmr subset and run
corr_res_sub <- list()
expr_final<-list()
prom_final<- list()
for(i in names(dmrs_final)){
    #select expressed genes of called dmrs and get same ordering of rows
    expr_sub <- expr[rownames(expr) %in%  dmrs_final[[i]]$ENSEMBL,]
    prom_sub <- as.data.frame(mcols(dmrs_final[[i]]))
    prom_sub <- prom_sub[prom_sub$ENSEMBL %in% rownames(expr_sub),]
    expr_sub$ENSEMBL <- rownames(expr_sub)
    merged <- inner_join(prom_sub, expr_sub, by="ENSEMBL")
    expr_sub <- merged[,colnames(expr_sub)]
    prom_sub <- merged[,colnames(prom_sub)]
    rownames(prom_sub)<- paste0(merged$SYMBOL, "_distanceToTSS:",merged$distanceToTSS, "_peakID:", merged$peakID)
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
    expr_sub <- expr_sub[,c("D129", "D213", "I217", "D217", "D360", "D123", "D124", "D117", "HCA_CB_Sampled")]
    colnames(expr_sub)<- c("D129", "D213", "I217", "D217", "D360", "D123", "D124", "D117", "cordblood")
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
    corr_res <- data.frame(p.value=unlist(lapply(corr, function(x)x$p.value)), estimate=unlist(lapply(corr, function(x)x$estimate)), symbol=merged_sub$SYMBOL, peakID= merged_sub$peakID)
    #corr_res$symbol <- prom[rownames(corr_res),]$SYMBOL
    #subset meaningful results and adjust for multiple testing
    corr_res_sub[[i]] <- corr_res[complete.cases(corr_res),]
    corr_res_sub[[i]]$p.adjust <- p.adjust(corr_res_sub[[i]]$p.value, "hochberg")
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$p.value, decreasing=FALSE),]
    print(head(corr_res_sub[[i]], 20))
    
    #plot results as histogramm
    corr_res_sub[[i]]$sign <- ifelse(corr_res_sub[[i]]$p.value<0.01, "P value < 0.01", "P value > 0.01")
    dir.create(file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs"), recursive=TRUE)
    write.table(corr_res_sub[[i]], file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "correlation_res.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    pdf(file=file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "correlation_res_histogram.pdf"), height=3.5, width=5)
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

    dir.create(file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs","scatter_plots"), recursive=TRUE)
    for(j in rownames(corr_res_sub_sig)){
        temp <- data.frame(Expression=as.numeric(expr_sub[j,]), Methylation=as.numeric(prom_sub_mean[j,]), Donor=colnames(prom_sub_mean))
        pdf(file=file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs","scatter_plots", paste0("scatter_plot_", j, ".pdf")), height=5, width=5)
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

saveRDS(corr_res_sub,file.path(analysis.dir,"corr_res_sub_inclNorm.rds"))
saveRDS(prom_final,file.path(analysis.dir,"corr_res_sub_methylation_values_inclNorm.rds"))
saveRDS(expr_final,file.path(analysis.dir,"corr_res_sub_expression_values_inclNorm.rds"))


corr_res_sub <- readRDS(file.path(analysis.dir,"corr_res_sub.rds"))
prom_final <- readRDS(prom_final,file.path(analysis.dir,"corr_res_sub_methylation_values.rds"))
expr_final <- readRDS(expr_final,file.path(analysis.dir,"corr_res_sub_expression_values.rds"))

#Plot histogramm of estimates
for(i in names(corr_res_sub)){
    corr_res_sub[[i]]$id <- rownames(corr_res_sub[[i]])
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$estimate, decreasing=TRUE),]

    pdf(file=file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "waterfall_estimate_all.pdf"), height=3.5, width=5)
    print(ggbarplot(corr_res_sub[[i]], x="id", y="estimate",ylab="correlation estimate\n[r]", fill="gray", color="gray") +rremove("x.axis")+rremove("x.grid")+rremove("xlab")+rremove("x.text")+rremove("x.ticks"))
    dev.off()
}
for(i in names(corr_res_sub)){
    corr_res_sub[[i]]$id <- rownames(corr_res_sub[[i]])
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$estimate, decreasing=TRUE),]
    corr_res_sub_sig <- corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value<0.01), ]

    pdf(file=file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "waterfall_estimate_pvalue0.01.pdf"), height=3.5, width=5)
    print(ggbarplot(corr_res_sub_sig, x="id", y="estimate",ylab="correlation estimate\n[r]", fill="gray", color="gray") +rremove("x.axis")+rremove("x.grid")+rremove("xlab")+rremove("x.text")+rremove("x.ticks"))
    dev.off()
}

#make heatmap of hits
anno_colors = list(
  Patient = c(D117 = "#0058b4", D129 = "#2188c9", 
            D217 = "#fbbb25", I217 = "#fca349", 
            D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", cordblood="#737373"), 
  Epigenotype = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4", normal ="#ababab")
  )
anno_heat <- data.frame(Patient=c("D117", "D129", "D213", "D217", "I217", "D124", "D123", "D360", "cordblood"),
    Epigenotype=c("LM", "LM", "HM", "IM", "IM", "HM","HM",  "HM", "normal"))
rownames(anno_heat) <-  anno_heat$Patient

for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
rna_sig <- expr_final[[i]][rownames(corr_res_sub_sig_ord),]
wgbs_sig <-  prom_final[[i]][rownames(corr_res_sub_sig_ord),]

nrow(wgbs_sig)==nrow(rna_sig)
#plot heatmap #estimate order
dir.create(file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps"))
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_RNA_signCorr0.01_estimateOrder_colClustered.pdf")
) 
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.01_estimateOrder_colClustered.pdf")
) 

#plot heatmap #clustering order independent
dir.create(file.path(analysis.dir, "heatmaps"))
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_RNA_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
)  
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
) 

#plot heatmap #clustering order based on rna
dir.create(file.path(analysis.dir, "heatmaps"))
p <- pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps","heatmap_RNA_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
)  
pheatmap::pheatmap(wgbs_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
) 
p <- pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
) 
pheatmap::pheatmap(rna_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
)  
}


#make heatmap of hits that overlap with degs
#load deg data
files <- list.files(path="/icgc/dkfzlsdf/analysis/C010/JMMLC/scRNA_Data/New/DEG/JMML-HCA_CB_Sampled-Epigenotype", pattern=".rds", full.names=TRUE)
files <- files[grep("HSC", files)]
temp <- sapply(strsplit(files, "-", ,fixed=TRUE), "[", 7)
temp <- sapply(strsplit(temp, ".", ,fixed=TRUE), "[", 1)
names(files)<- paste0(temp, "_vs_","cordblood_HSC")
deg_sig <- lapply(files, function(x){
    x <-readRDS(x)
    x$SYMBOL <-x$gene 
    x <- x[which(x$p_val_adj < 0.05),]
    x
    }
)
deg_sig <- list(groupHM=deg_sig$HM_vs_cordblood_HSC,groupLM=deg_sig$nonHM_vs_cordblood_HSC, all=rbind(deg_sig$HM_vs_cordblood_HSC, deg_sig$nonHM_vs_cordblood_HSC),
    groupHM_hierachy_DMRs_removed=deg_sig$HM_vs_cordblood_HSC,groupLM_hierachy_DMRs_removed=deg_sig$nonHM_vs_cordblood_HSC, all_hierachy_DMRs_removed= rbind(deg_sig$HM_vs_cordblood_HSC, deg_sig$nonHM_vs_cordblood_HSC))

for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
corr_res_sub_sig_ord <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% unique(deg_sig[[i]]$gene),]
rna_sig <- expr_final[[i]][rownames(corr_res_sub_sig_ord),]
wgbs_sig <-  prom_final[[i]][rownames(corr_res_sub_sig_ord),]

nrow(wgbs_sig)==nrow(rna_sig)
#plot heatmap #estimate order
dir.create(file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs"))
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_estimateOrder_colClustered.pdf")
) 
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_estimateOrder_colClustered.pdf")
) 

#plot heatmap #clustering order independent
pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
)  
pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredIndependent_colClustered.pdf")
) 


#plot heatmap #clustering order based on rna
p <- pheatmap::pheatmap(rna_sig,color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs","heatmap_RNA_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
)  
pheatmap::pheatmap(wgbs_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colClusteredRNA.pdf")
) 
#clustering based on wgbs
p <- pheatmap::pheatmap(wgbs_sig,color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
) 

pheatmap::pheatmap(rna_sig[p$tree_row$order, p$tree_col$order],color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colClusteredWGBS.pdf")
) 
#row clustering based on wgbs; columns based on methylation level
p <- pheatmap::pheatmap(wgbs_sig[,c( "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colMeanMeth.pdf")
) 

pheatmap::pheatmap(rna_sig[p$tree_row$order,c(  "cordblood","D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colMeanMeth.pdf")
) 


#row clustering based on rna; columns based on methylation level
p <- pheatmap::pheatmap(rna_sig[,c(  "cordblood","D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredRNA_colMeanMeth.pdf")
) 

pheatmap::pheatmap(wgbs_sig[p$tree_row$order,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360")],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colMeanMeth.pdf")
) 
}

#seperate positive and negative estimate
deg_overlap <- list()
for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.01),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
corr_res_sub_sig_ord <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% unique(deg_sig[[i]]$gene),]
deg_overlap[[i]] <- corr_res_sub_sig_ord

rna_sig_neg <- expr_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate < 0,]),, drop=FALSE]
wgbs_sig_neg <-  prom_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate < 0,]),, drop=FALSE]

rna_sig_pos <- expr_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate > 0,]),, drop=FALSE]
wgbs_sig_pos <-  prom_final[[i]][rownames(corr_res_sub_sig_ord[corr_res_sub_sig_ord$estimate > 0,]),, drop=FALSE]

#row clustering ordered; columns based on methylation level
#negative estimate
if(nrow(wgbs_sig_neg)> 1){
pheatmap::pheatmap(wgbs_sig_neg[,c( "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowOrdEstimate_colMeanMethh_NegativeEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_neg[,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_rowOrdEstimate_rowClusteredWGBS_colMeanMeth_NegativeEstimate.pdf")
) 
}
if(nrow(wgbs_sig_pos)>1){
#positive estimate
pheatmap::pheatmap(wgbs_sig_pos[,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_rowOrdEstimate_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_pos[,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA__rowOrdEstimate_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
)  
}
if(nrow(wgbs_sig_neg)> 1){
#row clustering based on wgbs; columns based on methylation level
#negative estimate
p <- pheatmap::pheatmap(wgbs_sig_neg[,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colMeanMethh_NegativeEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_neg[p$tree_row$order,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colMeanMeth_NegativeEstimate.pdf")
) 
}
if(nrow(wgbs_sig_pos)>1){
#positive estimate
p <- pheatmap::pheatmap(wgbs_sig_pos[,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
) 
pheatmap::pheatmap(rna_sig_pos[p$tree_row$order,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredWGBS_colMeanMeth_PositiveEstimate.pdf")
) 
}    
if(nrow(wgbs_sig_neg)> 1){
#row clustering based on rna; columns based on methylation level
#negative estimate
p <- pheatmap::pheatmap(rna_sig_neg[,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredRNA_colMeanMeth_NegativeEstimate.pdf")
) 

pheatmap::pheatmap(wgbs_sig_neg[p$tree_row$order,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colMeanMeth_NegativeEstimate.pdf")
) 
}
if(nrow(wgbs_sig_pos)>1){
#positive estimate
p <- pheatmap::pheatmap(rna_sig_pos[,c( "cordblood",  "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Greens"), cellwidth=25,cellheight=10,
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.01_rowClusteredRNA_colMeanMeth_PositiveEstimate.pdf")
) 

pheatmap::pheatmap(wgbs_sig_pos[p$tree_row$order,c(  "cordblood", "D117", "D129", "D217", "I217", "D213", "D124","D123", "D360"), drop=FALSE],color= RColorBrewer::brewer.pal(9,"Blues"),cellwidth=25,cellheight=10,
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=5,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis_inc_normals", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.01_rowClusteredRNA_colMeanMeth_PositiveEstimate.pdf")
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
saveRDS(deg_overlap,file.path(analysis.dir,"corr_res_sign0.05_overlapDEGs_up_downStrat_inclNorm.rds") )

paste0(unique(deg_overlap$upregulated$groupHM_hierachy_DMRs_removed$symbol), collapse=", ")
paste0(unique(deg_overlap$downregulated$groupHM_hierachy_DMRs_removed$symbol), collapse=", ")
paste0(unique(c(as.character(deg_overlap$downregulated$groupHM_hierachy_DMRs_removed$symbol), as.character(deg_overlap$upregulated$groupHM_hierachy_DMRs_removed$symbol))), collapse=", ")