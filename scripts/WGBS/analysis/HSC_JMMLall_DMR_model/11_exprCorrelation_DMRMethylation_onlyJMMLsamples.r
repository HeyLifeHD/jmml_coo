
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
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

dmrs_final <- lapply(dmrs_final, function(x){
    x$peakID <- 1:length(x) 
    x
})
#create further subset
dmrs_final_sub <- lapply(dmrs_final, function(x){
    x <- x[x %outside% dmrs_HSC_red,]
    x
})
dmrs_final_sub2 <- lapply(dmrs_final, function(x){
    x <- x[x %over% dmrs_HSC_red,]
    x
})
names(dmrs_final_sub)<- paste0(names(dmrs_final_sub), "_hierachy_DMRs_removed")
names(dmrs_final_sub2)<- paste0(names(dmrs_final), "_only_hierachy_DMRs")

lapply(dmrs_final, length)
lapply(dmrs_final_sub, length)
lapply(dmrs_final_sub2, length)
dmrs_final <- c(dmrs_final, dmrs_final_sub,dmrs_final_sub2)


#load expression data
expr <- readRDS("/omics/groups/OE0219/internal/JMMLC/scRNA_Data/New/DEG/Statistics/HSC.PseudoBulk.rds")

#loop over each dmr subset and run
corr_res_sub <- list()
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
    expr_sub <- expr_sub[, -9]
    prom_sub_mean <- prom_sub_mean[, -9]

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
    corr_res_sub[[i]]<- corr_res_sub[[i]][order(corr_res_sub[[i]]$p.adjust, decreasing=FALSE),]
    head(corr_res_sub[[i]], 20)
    
    #plot results as histogramm
    corr_res_sub[[i]]$sign <- ifelse(corr_res_sub[[i]]$p.adjust<0.1, "adj. p-value < 0.1", "adj. p-value > 0.1")
    dir.create(file.path(analysis.dir,i, "correlation_analysis", "promoter_meth_500bp"), recursive=TRUE)
    write.table(corr_res_sub[[i]], file.path(analysis.dir,i, "correlation_analysis", "promoter_meth_500bp", "correlation_res.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    pdf(file=file.path(analysis.dir,i, "correlation_analysis", "promoter_meth_500bp", "correlation_res_histogram.pdf"), height=3.5, width=5)
    print(ggpubr::gghistogram(corr_res_sub[[i]], x="estimate", bins=30, add="mean",
        xlab="Pearson correlation estimate\nAverage DMR methylation vs closest gene expression",ylab="Number of observations",
        fill="sign", palette=c("red", "black")))
    dev.off()

     #plot correlation of significant genes
    #corr_res_sub_sig <- corr_res_sub[[i]][which(corr_res_sub[[i]]$p.adjust<0.1), ]
    #plot correlation of top 50
    corr_res_sub_sig <- head(corr_res_sub[[i]],50)

    anno_col <- c(D117 = "#0058b4", D129 = "#2188c9", 
                D217 = "#fbbb25", I217 = "#fca349", 
                D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
                cordblood ="#ababab") 

    dir.create(file.path(analysis.dir,i, "correlation_analysis", "promoter_meth_500bp","scatter_plots"), recursive=TRUE)
    for(j in rownames(corr_res_sub_sig)){
        temp <- data.frame(Expression=as.numeric(expr_sub[j,]), Methylation=as.numeric(prom_sub_mean[j,]), Donor=colnames(prom_sub_mean))
        pdf(file=file.path(analysis.dir,i, "correlation_analysis", "promoter_meth_500bp","scatter_plots", paste0("scatter_plot_", j, ".pdf")), height=5, width=5)
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
    expr_sub <- expr_sub[, -9]
    prom_sub_mean <- prom_sub_mean[, -9]
    expr_final[[i]]<- expr_sub
    prom_final[[i]]<- prom_sub_mean

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
    head(corr_res_sub[[i]], 20)
    
    #plot results as histogramm
    corr_res_sub[[i]]$sign <- ifelse(corr_res_sub[[i]]$p.adjust<0.1, "adj. p-value < 0.1", "adj. p-value > 0.1")
    dir.create(file.path(analysis.dir,i, "correlation_analysis", "DMRs"), recursive=TRUE)
    write.table(corr_res_sub[[i]], file.path(analysis.dir,i, "correlation_analysis", "DMRs", "correlation_res.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    pdf(file=file.path(analysis.dir,i, "correlation_analysis", "DMRs", "correlation_res_histogram.pdf"), height=3.5, width=5)
    print(ggpubr::gghistogram(corr_res_sub[[i]], x="estimate", bins=30, add="mean",
        xlab="Pearson correlation estimate\nAverage DMR methylation vs closest gene expression",ylab="Number of observations",
        fill="sign", palette=c("red", "black")))
    dev.off()

     #plot correlation of significant genes
    corr_res_sub_sig <- corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value<0.05), ]
    #plot correlation of top 50
    #corr_res_sub_sig <- head(corr_res_sub[[i]],50)

    anno_col <- c(D117 = "#0058b4", D129 = "#2188c9", 
                D217 = "#fbbb25", I217 = "#fca349", 
                D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
                cordblood ="#ababab") 

    dir.create(file.path(analysis.dir,i, "correlation_analysis", "DMRs","scatter_plots"), recursive=TRUE)
    for(j in rownames(corr_res_sub_sig)){
        temp <- data.frame(Expression=as.numeric(expr_sub[j,]), Methylation=as.numeric(prom_sub_mean[j,]), Donor=colnames(prom_sub_mean))
        pdf(file=file.path(analysis.dir,i, "correlation_analysis", "DMRs","scatter_plots", paste0("scatter_plot_", j, ".pdf")), height=5, width=5)
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
saveRDS(corr_res_sub,file.path(analysis.dir,"corr_res_sub.rds") )



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
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.05),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
rna_sig <- expr_final[[i]][rownames(corr_res_sub_sig_ord),]
wgbs_sig <-  prom_final[[i]][rownames(corr_res_sub_sig_ord),]

nrow(wgbs_sig)==nrow(rna_sig)
#plot heatmap #estimate order
dir.create(file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps"))
pheatmap::pheatmap(rna_sig,#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_RNA_signCorr0.05_estimateOrder_colClustered.pdf")
) 
pheatmap::pheatmap(wgbs_sig,#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.05_estimateOrder_colClustered.pdf")
) 

#plot heatmap #clustering order independent
dir.create(file.path(analysis.dir, "heatmaps"))
pheatmap::pheatmap(rna_sig,#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_RNA_signCorr0.05_rowClusteredIndependent_colClustered.pdf")
)  
pheatmap::pheatmap(wgbs_sig,#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.05_rowClusteredIndependent_colClustered.pdf")
) 


#plot heatmap #clustering order based on rna
dir.create(file.path(analysis.dir, "heatmaps"))
p <- pheatmap::pheatmap(rna_sig,#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps","heatmap_RNA_signCorr0.05_rowClusteredRNA_colClusteredRNA.pdf")
)  
pheatmap::pheatmap(wgbs_sig[p$tree_row$order, p$tree_col$order],#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.05_rowClusteredRNA_colClusteredRNA.pdf")
) 

p <- pheatmap::pheatmap(wgbs_sig,#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_WGBS_signCorr0.05_rowClusteredWGBS_colClusteredWGBS.pdf")
) 

pheatmap::pheatmap(rna_sig[p$tree_row$order, p$tree_col$order],#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps", "heatmap_RNA_signCorr0.05_rowClusteredWGBS_colClusteredWGBS.pdf")
)  
}


#make heatmap of hits that overlap with degs
files <- list.files(path="/omics/groups/OE0219/internal/jmmlc_pbat/data/expression/200825_DEG/HM_vs_LM/", pattern=".csv", full.names=TRUE)
deg <- as.data.frame(data.table::fread(files))
deg$SYMBOL <- deg$gene
deg_sig <- deg[deg$p_val_adj < 0.05, ]

for(i in names(corr_res_sub)){
corr_res_sub_sig <-corr_res_sub[[i]][which(corr_res_sub[[i]]$p.value < 0.05),]
corr_res_sub_sig_ord <- corr_res_sub_sig[order(corr_res_sub_sig$estimate),]
corr_res_sub_sig_ord <- corr_res_sub_sig_ord[corr_res_sub_sig_ord$symbol %in% deg_sig$gene,]
rna_sig <- expr_final[[i]][rownames(corr_res_sub_sig_ord),]
wgbs_sig <-  prom_final[[i]][rownames(corr_res_sub_sig_ord),]

nrow(wgbs_sig)==nrow(rna_sig)
#plot heatmap #estimate order
dir.create(file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs"))
pheatmap::pheatmap(rna_sig,#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.05_estimateOrder_colClustered.pdf")
) 
pheatmap::pheatmap(wgbs_sig,#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.05_estimateOrder_colClustered.pdf")
) 

#plot heatmap #clustering order independent
pheatmap::pheatmap(rna_sig,#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.05_rowClusteredIndependent_colClustered.pdf")
)  
pheatmap::pheatmap(wgbs_sig,#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.05_rowClusteredIndependent_colClustered.pdf")
) 


#plot heatmap #clustering order based on rna
p <- pheatmap::pheatmap(rna_sig,#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs","heatmap_RNA_signCorr0.05_rowClusteredRNA_colClusteredRNA.pdf")
)  
pheatmap::pheatmap(wgbs_sig[p$tree_row$order, p$tree_col$order],#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.05_rowClusteredRNA_colClusteredRNA.pdf")
) 

p <- pheatmap::pheatmap(wgbs_sig,#color= RColorBrewer::brewer.pal(9,"Greens"), 
    scale="none", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=TRUE,cluster_cols=TRUE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_WGBS_signCorr0.05_rowClusteredWGBS_colClusteredWGBS.pdf")
) 

pheatmap::pheatmap(rna_sig[p$tree_row$order, p$tree_col$order],#color= RColorBrewer::brewer.pal(9,"Blues"),
    scale="row", show_colnames=F,
    annotation_col=anno_heat,
    fontsize_row=1,
    annotation_colors=anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,# clustering_distance_rows="correlation",
    show_rownames=T,
    filename= file.path(analysis.dir,i, "correlation_analysis", "DMRs", "heatmaps_degs", "heatmap_RNA_signCorr0.05_rowClusteredWGBS_colClusteredWGBS.pdf")
)  
}
