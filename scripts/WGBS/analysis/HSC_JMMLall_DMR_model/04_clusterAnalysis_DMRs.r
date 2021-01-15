
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
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/201005_DMR_tumor_vs_normal_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

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

#add reduced data for common analysis
dmrs_final <- lapply(dmrs_final, function(x){
    x$direction = ifelse(x$diff.Methy>0, "hyper", "hypo")
    x
})
dmrs_final$all <- dmrs_red
mcols(dmrs_final$all)$direction <- "hypo"


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
p <- list()
for(i in names(dmrs_final)){
    #all dmrs in each comparison with all samples
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    rownames(meth_dmr)<- paste0("Dmr.",1:nrow(meth_dmr),"_", dmrs_final[[i]]$SYMBOL)
    
    annovst <- as.data.frame(colData(bsseq_all))[, c("Tissue", "Celltype", "Patient", "Genotype", "Epigenotype")] 

    p[[i]] <- pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        clustering_distance_rows= "manhattan", clustering_method ="ward.D2",clustering_distance_cols= "manhattan", 
        scale="none",fontsize_row=5,  annotation_color=pbat_col#,filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_noScale_WarD2_Manhatten_both.pdf")
        )
    print(i)
}


#plot heatmap with clusters
clust_dmrs <- list()
for(i in names(p)){
    clust_dmrs[[i]]<- list()
    for(j in 3:10){
        meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
        rownames(meth_dmr)<- paste0("Dmr.",1:nrow(meth_dmr),"_", dmrs_final[[i]]$SYMBOL)
        #define k 
        cut <- j
        #cut tree
        clust_dmrs[[i]][[j]] <- cutree(p[[i]]$tree_row, k=cut)
         clust_dmrs[[i]][[j]] <- as.data.frame( clust_dmrs[[i]][[j]])
        colnames( clust_dmrs[[i]][[j]])<- "cluster"
         clust_dmrs[[i]][[j]]$cluster <- as.character( clust_dmrs[[i]][[j]]$cluster)
        #add color row annotation
        color_cluster <- randomcoloR::distinctColorPalette(length(unique( clust_dmrs[[i]][[j]]$cluster)))
        names(color_cluster)<- unique( clust_dmrs[[i]][[j]]$cluster)
        pbat_col = list(
            Tissue = c(tumor = "#99a637", prenatal = "#252525", cordblood = "#737373", adult_bonemarrow = "#ababab"), 
            Celltype = c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"), 
            Genotype = c(neg = "#529a51", KRAS = "#99a637", PTPN11 = "#007458", wildtype ="#ababab"),
            Patient = c(D117 = "#0058b4", D129 = "#2188c9", 
                        D217 = "#fbbb25", I217 = "#fca349", 
                        D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
                        normal ="#ababab"), 
            Epigenotype = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4", wildtype ="#ababab"),
            cluster=color_cluster
            )
        #heatmap
        dir.create(file.path(analysis.dir,i, "visualization","cluster_anal"))
        pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
                clustering_distance_rows= "manhattan", clustering_method ="ward.D2",clustering_distance_cols= "manhattan", 
                scale="none",fontsize_row=5,  annotation_color=pbat_col, cutree_rows =cut, annotation_row=clust_dmrs[[i]][[j]], 
                file=file.path(analysis.dir,i, "visualization","cluster_anal", paste0("Heatmap_DMRs_noScale_WarD2_Manhatten_both_clust_", cut,".pdf")))
    print(j)
    }
    #add cluster annotation to dmr list
    # start with 6 for all dmr sets
    dmrs_final[[i]]$cluster <- clust_dmrs[[i]][[6]]$cluster
    print(i)
}
saveRDS(dmrs_final, )
#seperate dmrs per cluster and plot average methylation per group
dmrs_final_split <- lapply(dmrs_final, function(x)split(x, x$cluster))
#extract average methylaiton per cluster
pheno <- pData(bsseq_all)
meth_median <- list()
for(i in names(dmrs_final_split)){
    meth_median[[i]] <- list()
    for(j in names(dmrs_final_split[[i]])){
        meth_median[[i]][[j]] <-  pheno
        meth_median[[i]][[j]]$cluster <- j
        meth_median[[i]][[j]]$Methylation_level <-  colMedians(as.matrix(as.data.frame(dmrs_final_split[[i]][[j]])[, rownames(pheno)]), na.rm=TRUE)
    }
    meth_median[[i]] <- do.call("rbind", meth_median[[i]])
}

#Plot average Methylation 
for(i in names(meth_median)){
#plotting
dir.create(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6"))
meth_median[[i]] <- as.data.frame(meth_median[[i]])
meth_median[[i]]$cluster <- paste0("cluster_", meth_median[[i]]$cluster)
#Donor
compare_means(Methylation_level ~ Donor ,  data = meth_median[[i]], method = "t.test")
#med <- aggregate(meth_median[[i]][,"Methylation_level",drop=FALSE], list(meth_median[[i]]$Donor), median)
#ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
ord <- c("cordblood", "D117", "D129", "D217", "I217", "D213", "D123", "D360", "D124")

pdf(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6",paste0(i, "_AvMeth_Donor.pdf")), height=7, width=7)
print(ggpubr::ggboxplot(meth_median[[i]], x="Donor", y="Methylation_level", color="Donor", facet.by="cluster",ylim=c(0,1), order=ord,
    palette =c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#Epigenotype
compare_means(Methylation_level ~ Epigenotype ,  data = meth_median[[i]], method = "t.test")
#med <- aggregate(meth_median[[i]][,"Methylation_level",drop=FALSE], list(meth_median[[i]]$Epigenotype), median)
#ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
ord <- c("wildtype", "LM", "IM", "HM")

pdf(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6",paste0(i, "_AvMeth_Epigenotype.pdf")), height=7, width=7)
print(ggpubr::ggboxplot(meth_median[[i]], x="Epigenotype", y="Methylation_level",  facet.by="cluster", ylim=c(0,1), order=ord,
    color="Epigenotype", palette =c(wildtype ="#ababab", LM = "#0058b4", IM = "#fbbb25", HM = "#c33126"),title=i, ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#genotype
compare_means(Methylation_level ~ Genotype,  data = meth_median[[i]], method = "t.test")
med <- aggregate(meth_median[[i]][,"Methylation_level",drop=FALSE], list(meth_median[[i]]$Genotype), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6",paste0(i, "_AvMeth_Genotype.pdf")), height=7, width=7)
print(ggpubr::ggboxplot(meth_median[[i]], x="Genotype", y="Methylation_level", color="Genotype",facet.by="cluster", ylim=c(0,1), order=ord,
    palette =c(wildtype ="#ababab", neg = "#529a51", KRAS = "#99a637", PTPN11 = "#007458"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#tumor
compare_means(Methylation_level ~ Tumor,  data = meth_median[[i]], method = "t.test")
med <- aggregate(meth_median[[i]][,"Methylation_level",drop=FALSE], list(meth_median[[i]]$Tumor), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6",paste0(i, "_AvMeth_tumor.pdf")), height=7, width=7)
print(ggpubr::ggboxplot(meth_median[[i]], x="Tumor", y="Methylation_level", color="Tumor",facet.by="cluster",ylim=c(0,1), order=ord,
    palette =c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#c6a27f", 
                       HSC_adult ="#252525", MPP_adult = "#737373", CMP_adult = "#ffb86f", GMP_adult = "#e27e37", MEP_adult="#d2624a",
                       HSC_CB ="#252525"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#Celltype
compare_means(Methylation_level ~ Celltype,  data = meth_median[[i]], method = "t.test")
med <- aggregate(meth_median[[i]][,"Methylation_level",drop=FALSE], list(meth_median[[i]]$Celltype), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6",paste0(i, "_AvMeth_Celltype.pdf")), height=7, width=7)
print(ggpubr::ggboxplot(meth_median[[i]], x="Celltype", y="Methylation_level", color="Celltype",facet.by="cluster",ylim=c(0,1), order=ord,
    palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
#Sample_Type
compare_means(Methylation_level ~ Sample_Type,  data = meth_median[[i]], method = "t.test")
med <- aggregate(meth_median[[i]][,"Methylation_level",drop=FALSE], list(meth_median[[i]]$Sample_Type), median)
ord <- med[order(med$Methylation_level, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,i, "visualization","cluster_anal","clust_6",paste0(i, "_AvMeth_Sample_Type.pdf")), height=7, width=7)
print(ggpubr::ggboxplot(meth_median[[i]], x="Sample_Type", y="Methylation_level", color="Sample_Type",facet.by="cluster",ylim=c(0,1), order=ord,
    palette =c(normal = "#ababab", tumor = "#99a637"),title=i,ylab="Average methylation levels of DMRs",
    add = "jitter") +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45))
dev.off()
    
}
saveRDS(dmrs_final,file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_clust6_annotattion.rds"))

#prepare clusters for metascape analysis
temp <- dmrs_final_split$Sample_Typetumor_hierachy_DMRs_removed
paste0(unique(temp$"1"$SYMBOL), collapse=", ")
paste0(unique(temp$"2"$SYMBOL), collapse=", ")
paste0(unique(temp$"3"$SYMBOL), collapse=", ")
paste0(unique(temp$"4"$SYMBOL), collapse=", ")
paste0(unique(temp$"5"$SYMBOL), collapse=", ")
paste0(unique(temp$"6"$SYMBOL), collapse=", ")

#continue with plotting
library(readxl)
library(ggpubr)
library(RColorBrewer)
#function
plotMeta <- function(excel_path, n=10, output_path, color="#CD534CFF", height=3.5, width=3.5){
dataset <- read_excel(excel_path, sheet=2)
datasetsub <- dataset[grep("Summary", dataset$GroupID), ]
datasetsub <- head(datasetsub, n)
datasetsub$LogP<- abs(datasetsub$LogP)
datasetsub$Description<- rev(datasetsub$Description)
datasetsub$LogP<- rev(datasetsub$LogP)
datasetsub$Term<- rev(datasetsub$Term)
datasetsub$description <- as.factor(1:nrow(datasetsub))
 p <- ggplot(datasetsub, aes(x=description)) +
geom_bar(aes(x=description, y=LogP, fill = "-log10(p.value)"), stat = "identity", width=0.25)
 p<-p+ annotate("text",x=as.integer(datasetsub$description)+0.35, y=0.01, 
                label=datasetsub$Description, size=4, hjust = 0)
 p <- p + scale_y_continuous( expand = c(0, 0), name = "-log10(p-value)")
 p <- p + scale_fill_manual(values = c( color))
 p<-p+coord_flip() 
 p<-p+ ggpubr:::theme_pubr()
 p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
 p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +rremove(c("legend"))+rremove(c("ylab"))
 p$labels$fill <- ""
 pdf(output_path, height=height, width=width)
 print(p)
 dev.off()
}

#plot results
files <- list.files(file.path(analysis.dir,"Sample_Typetumor_hierachy_DMRs_removed", "visualization","cluster_anal","clust_6"), 
    pattern="metascape_result.xlsx", full.names=TRUE, recursive=TRUE)
paths <- gsub("metascape_result.xlsx", "", files)
files_abs<- files[grep("clust", files)]
paths_abs<- paths[grep("clust", paths)]
#plotting
for(numb in c(3,5,10,15)){
for(i in 1:length(files_abs)){
    plotMeta(files_abs[i], 
        file.path(paths_abs[i], paste0("metascape_vis_abs_",numb,"_1.pdf")), 
        n=numb, width=5, height=5,
        color= "gray")
}
}
