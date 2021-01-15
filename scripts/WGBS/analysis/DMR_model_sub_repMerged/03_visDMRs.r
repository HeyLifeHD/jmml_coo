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
output.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"

#load data
bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
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

#new output directory
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_model_sub_repMerged"
dir.create(analysis.dir)

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
dmrs_final$all <- dmrs_red
mcols(dmrs_final$all)$direction <- "hypo"

#subset dmrs with a abs(difference) >0.1 
#dmrs_final_sub <- lapply(dmrs_final, function(x){
#    x <- x[abs(x$diff.Methy)>0.3,]
#    x
#})
#lapply(dmrs_final_sub, function(x)table(x$direction))

#export text files
for(i in names(dmrs_final)){
    dir.create(file.path(analysis.dir,i))
    write.table(as.data.frame(dmrs_final[[i]]),file.path(analysis.dir,i, paste0(i,"_","dmrs_final.txt")),row.names = FALSE, quote=FALSE, sep="\t")
}

#width of dmrs
for(i in names(dmrs_final)){
    temp <- width(dmrs_final[[i]])
    dir.create(file.path(analysis.dir, i, "visualization"))
    pdf(file.path(analysis.dir, i, "visualization", paste0("DMR_width_histo",i, ".pdf")), height=4, width=4)
    print(gghistogram(as.data.frame(temp), x="temp", fill="grey", add_density=TRUE,add="mean",bins=50, rug=TRUE, xlab="Width of DMRs [bp]", ylab="# of DMRs")+xscale("log2"))
    dev.off()
    print(mean(temp))
    dmrs_final[[i]]$comparison <- i
}
width <- lapply(dmrs_final, function(x){
    width <- as.data.frame(width(x))
    width$comparison <- x$comparison
    width
})
width$all <- NULL
width <- do.call("rbind", width)
colnames(width)<-c("width", "comparison")
pdf(file.path(analysis.dir,"all","visualization", "DMr_width_allComparison.pdf"), height=4, width=7)
gghistogram(width, x="width", color="comparison", add="mean",add_density=TRUE,bins=50, rug=TRUE, xlab="Width of DMRs [bp]", ylab="# of DMRs", palette="jco",legend="right") +xscale("log2")
dev.off()


#pca of all original dmrs
#PCA
for(i in names(dmrs_final)){
     meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    
    meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
    ir.pca <- prcomp(t(as.matrix(meth_dmr)),
                 center = T,
                 scale. = F) 
    
    print(i)
    print(summary(ir.pca))
    pc <- ir.pca$x
    pc<- as.data.frame(cbind(pc , pData(bsseq_all)))

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


#umap of all samples
for(i in names(dmrs_final)){
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
    umap <- umap::umap(t(as.matrix(meth_dmr)))
    umap_layout <- umap$layout
    colnames(umap_layout)<- c("UMAP1", "UMAP2")
    umap_layout <- as.data.frame(cbind(umap_layout, colData(bsseq_all)))

    #patient and tumor
    pdf(file.path(analysis.dir, i, "visualization","UMAP_PatientTumor_onlyNormal.pdf"),height = 5, width = 5)
    print(ggscatter(umap_layout, x="UMAP1", y="UMAP2",
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
    print(ggscatter(umap_layout,x="UMAP1", y="UMAP2",
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
    print(ggscatter(umap_layout,x="UMAP1", y="UMAP2",
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
    print(ggscatter(umap_layout,  x="UMAP1", y="UMAP2",
          color = "Celltype", shape = "Tissue",label="Patient",repel=TRUE,
          ellipse = F , mean.point = FALSE,palette= c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),
            star.plot = F) +
            theme(legend.position="right",legend.title = element_text(, size=10, 
                                        face="bold")))
    dev.off()
    print(i)
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
seeed(42)
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



#location of dmrs
library(circlize)
for(i in names(dmrs_final)){
    dmrs_final_df <- as.data.frame(dmrs_final[[i]])
    pdf(file.path(analysis.dir,i,"visualization", paste0("DMRlocation_CirclePlot",".pdf")))
    circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", c(1:22, "X", "Y")))
    #bed_list_en = list(dmrs_final_df[[i]][dmrs_final_df[[i]]$direction =="hypo",], dmrs_final_df[[i]][dmrs_final_df[[i]]$direction =="hyper",])
    #circos.genomicRainfall(bed_list_en, pch = 16, cex = 0.8, col = c("#0000FF80", "#FF000080"))
    circos.genomicDensity(dmrs_final_df[dmrs_final_df$direction =="hypo",], col = c("#FF000080"), track.height = 0.1)
    circos.genomicDensity(dmrs_final_df[dmrs_final_df$direction =="hyper",], col = c("#0000FF80"), track.height = 0.1)
    dev.off()
    circos.clear()
}


#pie plot of dmr dsistribution
#define colors
#col <- brewer.pal(3,"RdBu")
col <- c("red","blue")
#For loop
for (comp in names(dmrs_final)){
# calculate Distribution
table <- as.data.frame(table(dmrs_final[[comp]]$direction))
labs <- paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq)
pdf(file.path(analysis.dir,comp,"visualization",paste0("Direction_Pie_DMRs.pdf")), height=3.5, width=3.5)
print(ggpie(table, "Freq", fill="Var1", palette=col, label=labs, lab.pos = "in", main="DMR Analysis",submain=comp,  lab.font = c(5, "bold", "white")) + rremove("legend"))
dev.off()
#Distance to TSS
#dis <- dmrs_final[[comp]]$distanceToTSS
#dis <- as.data.frame(dis)
#dis <- dis/1000
#dis <-dis[dis$dis > (-14) & dis$dis < (14),]
#dis <- as.data.frame(dis)
#histo
#pdf(file.path(analysis.dir,comp,"visualization",paste0("Distance_to_TSS_Histogram.pdf")), height=3.5)
#print(gghistogram(dis, x="dis", fill="lightgray", bins=40, rug = T, xlab="Distance from TSS (kb)", ylab="# of DMRs", xlim=c(-10,10)
#)  + geom_vline(xintercept=0, linetype = 1, color="red", size=2 ))
#dev.off() 
#stratified barplot
dis <- as.data.frame(dmrs_final[[comp]]$diff.Methy)
colnames(dis)<- "dis"
df <- data.frame( ">-100"=sum(dis< (-100)),  "-10 to -100"=sum(dis<(-10) & dis>(-100)), "-1 to -10"=sum(dis<(-1) & dis>(-10)), "0 to -1"=sum(dis<0 & dis>(-1)),
"0 to 1"=sum(dis>0 & dis<1), "1 to 10"=sum(dis>1 & dis<10), "10 to 100"=sum(dis>10 & dis<100), ">100"=sum(dis<100))
colnames(df)<- c( ">-100","-10 to -100",  "-1 to -10","0 to -1", "0 to 1" , "1 to 10", "10 to 100", ">100")
df <- data.frame(distance= colnames(df), quant=t(df)[,1])
df$quant <- df$quant/sum(df$quant)
pdf(file.path(analysis.dir,comp,"visualization", paste0("Distance_to_TSS_barplot.pdf")), height=3.5)
print(ggbarplot(df, x="distance", y="quant",  fill= "gray", ylab="Number of DMRs", main="CpG analysis")+rotate_x_text(angle = 90))
dev.off()
#Difference in Methylation
df <- cut(dis$dis, breaks=c((-0.8), (-0.6) ,(-0.4), (-0.2),0, 0.2,0.4,0.6,0.8), labels=c("-0.8 to -0.6","-0.6 to -0.4", "-0.4 to -0.2","-0.2 to -0","0 to 0.2","0.2 to 0.4","0.4 to 0.6","0.6 to 0.9" ))
df<-as.data.frame(t(table(df)))
df$Var1 <- c(rep("hypo",4), rep("hyper",4))
pdf(file.path(analysis.dir,comp,"visualization",paste0("DifferenceMeth_stratiefied.pdf")), height=3.5, width=3.5)
print(ggbarplot(df ,x="df", y="Freq", fill="Var1", palette=col, xlab="Difference in methylation", ylab="Number of DMRs") +rotate_x_text(angle = 45)+ rremove("legend"))
dev.off()
}