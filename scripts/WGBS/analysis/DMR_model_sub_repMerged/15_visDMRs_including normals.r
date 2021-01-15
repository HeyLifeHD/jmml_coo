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

#Directories
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
input_DMR.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"
analysis.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_model_sub_repMerged_inclNormals"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(input_DMR.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(input_DMR.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#change annotation


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
names(dmrs_final)

#add methylation levels
#add methylation difference information
meth_dmrs<- list()
for(comp in names(dmrs_final)){
    #mcols(dmrs_final[[comp]])[,rownames(pData(bsseq_all)[rownames(pData(bsseq_all)) %in%colnames( mcols(dmrs_final[[comp]])) ]]
  mcols(dmrs_final[[comp]])[,4:ncol(mcols(dmrs_final[[comp]]))]<-NULL
  mcols(dmrs_final[[comp]]) <- bsseq::getMeth(bsseq_all, regions= dmrs_final[[comp]], type = "raw", what=c("perRegion"))
  idx <- complete.cases( mcols(dmrs_final[[comp]]) )
  dmrs_final[[comp]] <- dmrs_final[[comp]][idx,]
  print(comp)
}
#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
dmrs_final$all$direction <- "hypo"

#export text files
for(i in names(dmrs_final)){
    dir.create(file.path(analysis.dir,i))
    write.table(as.data.frame(dmrs_final[[i]]),file.path(analysis.dir,i, paste0(i,"_","dmrs_final.txt")),row.names = FALSE, quote=FALSE, sep="\t")
}


#pca of all original dmrs
#PCA
for(i in names(dmrs_final)){
    dir.create(file.path(analysis.dir, i, "visualization"))
     meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    
    meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
    ir.pca <- prcomp(t(as.matrix(meth_dmr)),
                 center = T,
                 scale. = F) 
    
    print(i)
    print(summary(ir.pca))
    pc <- ir.pca$x
    pc<- as.data.frame(cbind(pc , pData(bsseq_all)))
    dir.create(file.path(analysis.dir, i, "visualization"), recursive=TRUE)
    #Donor and Epigenotype
    pdf(file.path(analysis.dir, i, "visualization","PC12_DonorEpigenotype.pdf"),height = 5, width = 5)
    print(ggscatter(pc,x="PC1", y="PC2",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
    dev.off()
    pdf(file.path(analysis.dir, i, "visualization","PC23_DonorEpigenotype.pdf"),height = 5, width = 5)
    print(ggscatter(pc,x="PC2", y="PC3",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
            star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
            "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    pdf(file.path(analysis.dir, i, "visualization","PC34_DonorEpigenotype.pdf"),height = 5, width = 5)
    print(ggscatter(pc, x="PC3", y="PC4",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
            star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
            "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    #Donor and Genotype
    pdf(file.path(analysis.dir, i, "visualization","PC12_DonorGenotype.pdf"),height = 5, width = 5)
    print(ggscatter(pc,x="PC1", y="PC2",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
            star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
            "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    pdf(file.path(analysis.dir, i, "visualization","PC23_DonorGenotype.pdf"),height = 5, width = 5)
    print(ggscatter(pc, x="PC2", y="PC3",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
            star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
            "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    pdf(file.path(analysis.dir, i, "visualization","PC34_DonorGenotype.pdf"),height = 5, width = 5)
    print(ggscatter(pc, x="PC3", y="PC4",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
            star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
            "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    #Celltype and Tissue
    pdf(file.path(analysis.dir, i, "visualization","PC12_CelltypeTissue.pdf"),height = 5, width = 5)
    print(ggscatter(pc,  x="PC1", y="PC2",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
            star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
            "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    pdf(file.path(analysis.dir, i, "visualization","PC23_CelltypeTissue.pdf"),height = 5, width = 5)
    print(ggscatter(pc,  x="PC2", y="PC3",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HHSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
            star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
            "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    pdf(file.path(analysis.dir, i, "visualization","PC34_CelltypeTissue.pdf"),height = 5, width = 5)
    print(ggscatter(pc, x="PC3", y="PC4",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
            star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
            "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
}


#n of all samples
for(i in names(dmrs_final)){
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
    n <- n::n(t(as.matrix(meth_dmr)))
    n_layout <- n$layout
    colnames(n_layout)<- c("UMAP1", "UMAP2")
    n_layout <- as.data.frame(cbind(n_layout, colData(bsseq_all)))

    #patient and tumor
    pdf(file.path(analysis.dir, i, "visualization","UMAP_PatientTumor_onlyNormal.pdf"),height = 5, width = 5)
    print(ggscatter(n_layout, x="UMAP1", y="UMAP2",
            color = "Patient", shape = "Tumor",#size="Protocol",
            ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                        D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000", 
                        normal ="#e9e9e9"),
            star.plot = F) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()

    #Donor and Epigenotype
    pdf(file.path(analysis.dir, i, "visualization","UMAP_DonorEpigenotype.pdf"),height = 5, width = 5)
    print(ggscatter(n_layout,x="UMAP1", y="UMAP2",
          color = "Donor", shape = "Epigenotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
          star.plot = F) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold")))
    dev.off()

    #Donor and Genotype
    pdf(file.path(analysis.dir, i, "visualization","UMAP_DonorGenotype.pdf"),height = 5, width = 5)
    print(ggscatter(n_layout,x="UMAP1", y="UMAP2",
          color = "Donor", shape = "Genotype",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),
            star.plot = F) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()

    #Celltype and Tissue
    pdf(file.path(analysis.dir, i, "visualization","UMAP_CelltypeTissue.pdf"),height = 5, width = 5)
    print(ggscatter(n_layout,  x="UMAP1", y="UMAP2",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
            star.plot = F) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
           
}



#Sample Clustering
for( i in names(dmrs_final)){
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    drld <- dist(t(as.matrix(meth_dmr)))
    hrld<- hclust(drld)
    dend<- hrld%>% as.dendrogram 

    anno <- pData(bsseq_all)
    #col1 <- hrld$labels
    #names(col1)<- anno[col1,]$group
    #col1 <- col1[order.dendrogram(dend)]
    #col1 <- ifelse(names(col1)=="TAM_tumor", "#EFC00099", ifelse(names(col1)=="BMDM_tumor","#0073C299",  ifelse(names(col1)=="TAM_healthy","#CD534CFF","#86868699")))
    dend <- dend %>% 
    #set("branches_k_color", k=4, c("#EFC00099", "#0073C299", "#CD534CFF","#86868699" )) %>%
    set("branches_lwd", 2) %>%
    #set("labels_colors",col1) %>% 
    set("labels_cex", .6 )%>%
    set("leaves_pch", 19)%>% 
    set("leaves_cex", 1.5)#%>% 
    #set("leaves_col", col1)
    
    pdf(file.path(analysis.dir, i, "visualization", "Clustering_DMR_meth.pdf"), height = 6, width = 5)
    print(dend %>% plot)
    dev.off()

}


#Heatmap of DMRs
#get contrasts
seed(42)
pbat_col = list(
  Tissue = c(tumor = "#99a637", prenatal = "#252525", cordblood = "#737373", adult_bonemarrow = "#ababab"), 
  Celltype = c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"), 
  Genotype = c(neg = "#529a51", KRAS = "#99a637", PTPN11 = "#007458", wildtype ="#ababab"),
  Patient = c(D117 = "#0058b4", D129 = "#2188c9", 
            D217 = "#fbbb25", I217 = "#fca349", 
            D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
            normal ="#ababab"), 
  Epigenotype = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4", wildtype ="#ababab")
  )
for(i in names(dmrs_final)){
    #all dmrs in each comparison with all samples
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    rownames(meth_dmr)<- paste0("Dmr.",1:nrow(meth_dmr),"_", dmrs_final[[i]]$SYMBOL)

    annovst <- as.data.frame(colData(bsseq_all))[, c("Tissue", "Celltype", "Patient", "Genotype", "Epigenotype")] 

            
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
        scale="row",fontsize_row=5,  annotation_color=pbat_col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_rowScale_WarD2_Manhatten.pdf"))
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2",
        scale="none",fontsize_row=5,  annotation_color=pbat_col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_noScale_WarD2_Manhatten.pdf"))
    
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2", clustering_distance_cols= "manhattan",
        scale="row",fontsize_row=5,  annotation_color=pbat_col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_rowScale_WarD2_Manhatten_both.pdf"))
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2",clustering_distance_cols= "manhattan", 
        scale="none",fontsize_row=5,  annotation_color=pbat_col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_noScale_WarD2_Manhatten_both.pdf"))
    print(i)
}


