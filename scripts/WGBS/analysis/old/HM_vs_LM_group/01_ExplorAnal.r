##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Exploratory analysis

#libraries
library(bsseq)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggpubr)
library(pheatmap)
library(randomcoloR)

#directories
odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_ExplorAnal"
dir.create(analysis.dir)
#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_all.rds"))
#subset outliers
bsseq_all <- bsseq_all[,!colnames(bsseq_all) %in% c("tumor11_JMMLC_D129_1","tumor11_JMMLC_D117", "tumor00_JMMLC_D123", "tumor10_JMMLC_D123",
    "tumor01_JMMLC_D123", "tumor00_JMMLC_D124", "tumor01_JMMLC_D117")] 
dim(bsseq_all)
pheno<- pData(bsseq_all)
pheno$patient <- droplevels(pheno$patient)
pheno$Epigenotype <- as.factor(pheno$Epigenotype )
pheno$Epigenotype <-factor(pheno$Epigenotype, levels = c("LM", "HM"))
pheno$Genotype <- as.factor(pheno$Genotype)
pheno$Genotype <-factor(pheno$Genotype, levels = c("neg", "KRAS", "PTPN11"))
pData(bsseq_all) <- pheno
saveRDS(bsseq_all, file.path(input.dir ,"bsseq", "bsseq_all_snpfil_sub_cov.rds"))
bsseq_all <- readRDS( file.path(input.dir ,"bsseq", "bsseq_all_snpfil_sub_cov.rds"))

#Plot average Methylation 
## extract the methylation values
pheno <- colData(bsseq_all)
meth_per_cpg <- bsseq::getMeth(bsseq_all, type = "raw")
#saveRDS(meth_per_cpg, file.path(output.dir, "bsseq","filtered_meth_matrix.rds"))#
meth_per_cpg <- as.data.frame(meth_per_cpg)
av.meth_per_cpg<- colMeans(meth_per_cpg, na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Epigenotype=as.character(pheno$Epigenotype), Patient=as.character(pheno$Patient), 
    Genotype=as.character(pheno$Genotype), Tumor= as.character(pheno$tumor))

#plotting
#Epigenotype
compare_means(AverageMethylation ~ Epigenotype,  data = av.meth_per_cpg, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_Epigenotype.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,1), palette =c( HM = "#b5241c", LM = "#0041a5"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()
#Patient
compare_means(AverageMethylation ~ Patient,  data = av.meth_per_cpg, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_Patient.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient", y="AverageMethylation", color="Patient",ylim=c(0,1), 
    palette =c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()

compare_means(AverageMethylation ~ Genotype,  data = av.meth_per_cpg, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_Genotype.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Genotype", y="AverageMethylation", color="Genotype",ylim=c(0,1),
    palette =c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()
#tumor
compare_means(AverageMethylation ~ Tumor,  data = av.meth_per_cpg, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_tumor.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor", y="AverageMethylation", color="Tumor",ylim=c(0,1), 
    palette =c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    #stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()



#Plot average Methylation of CpGislans
#cpg locations
library(RnBeads)
cgi <- unlist(rnb.get.annotation("cpgislands"))

## extract the methylation values
meth_per_cpg_cpgi <- bsseq::getMeth(bsseq_all, type = "raw", what="perRegion", regions=cgi)
colnames(meth_per_cpg_cpgi)<- rownames(pheno)
av.meth_per_cpg_cpgi<- colMeans(meth_per_cpg_cpgi, na.rm = TRUE)
av.meth_per_cpg_cpgi <- data.frame(AverageMethylation=av.meth_per_cpg_cpgi, Epigenotype=as.character(pheno$Epigenotype), Patient=as.character(pheno$Patient), 
    Genotype=as.character(pheno$Genotype), Tumor= as.character(pheno$tumor))

#plotting
#Epigenotype
compare_means(AverageMethylation ~ Epigenotype,  data = av.meth_per_cpg_cpgi, method = "t.test")
my_comparisons <- list(c("HM", "LM"))
pdf(file.path(analysis.dir ,"AvMeth_CpGisland_Epigenotype.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg_cpgi, x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,1),
    palette =c(HM = "#b5241c", LM = "#0041a5"),
    add = "jitter")  +  
    stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c( .91)) + 
    stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()
#Patient
compare_means(AverageMethylation ~ Patient,  data = av.meth_per_cpg_cpgi, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_CpGisland_Patient.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg_cpgi, x="Patient", y="AverageMethylation", color="Patient",ylim=c(0,1), 
    palette =c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()
#Genotype

compare_means(AverageMethylation ~ Genotype,  data = av.meth_per_cpg_cpgi, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_CpGisland_Genotype.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg_cpgi, x="Genotype", y="AverageMethylation", color="Genotype",ylim=c(0,1),
    palette =c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()
#tumor
compare_means(AverageMethylation ~ Tumor,  data = av.meth_per_cpg_cpgi, method = "t.test")
pdf(file.path(analysis.dir ,"AvMeth_CpGisland_tumor.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg_cpgi, x="Tumor", y="AverageMethylation", color="Tumor",ylim=c(0,1),
    palette = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d"),
    add = "jitter")  +  
    #stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.80, .85, .91, .95)) + 
    stat_compare_means(label.y = .25) +
    rremove("legend")+
    rremove("xlab")+
    rotate_x_text(angle = 45)
dev.off()


#Heatmap of most variable CpGs
#meth_per_cpg <- bsseq::getMeth(bsseq_all, type = "raw")
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),20000)
meth_per_cpg<- as.data.frame(meth_per_cpg)
matrld <- meth_per_cpg[topVarCpGs,]
pbat_col = list(
  Patient = c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
              D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"), 
  Genotype = c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458"),
  Epigenotype = c(HM = "#b5241c", LM = "#0041a5"), 
  tumor = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d"))
annovst <- as.data.frame(colData(bsseq_all))[, c("Epigenotype", "Patient", "tumor", "Genotype")] 

pheatmap(matrld,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,
    show_colnames=FALSE, scale="row",fontsize_row=5,  annotation_colors=pbat_col,
    filename=file.path(analysis.dir,"Heatmap20000mvCpG.pdf"))
pheatmap(matrld,  annotation_col=as.data.frame(annovst), annotation_colors=pbat_col,
    show_rownames=FALSE,show_colnames=FALSE, scale="none",fontsize_row=5,  
    filename=file.path(analysis.dir,"Heatmap20000mvCpG_noScale.pdf"))


#Heatmap of most variable average Promoter Methylation
#Heatmap on Prtomoter
library("EnsDb.Hsapiens.v75")
ENS<-EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
GENE = genes(ENS)

#1500 upstream and 500 downstream as standard
#20000mv
Promoter<- promoters(GENE, upstream = 1500, downstream = 500)
PromoterMeth <- bsseq::getMeth(bsseq_all, regions=Promoter, what="perRegion", type="raw")
row.names(PromoterMeth)<-Promoter$symbol
PromoterMeth<- as.matrix(PromoterMeth)
topVarGenesRld<- head(order(rowVars(PromoterMeth), decreasing=TRUE),20000)
anno <- Promoter
anno<-as.data.frame(anno)
annotop <- anno[topVarGenesRld,]
matrld <- PromoterMeth[topVarGenesRld,]
pheatmap(matrld,  annotation_col=as.data.frame(annovst),
    show_rownames=FALSE,show_colnames=FALSE, scale="row",fontsize_row=5, annotation_colors=pbat_col,, 
    filename=file.path(analysis.dir,"Heatmap_20000mvPromoterMeth.pdf"))
pheatmap(matrld,  annotation_col=as.data.frame(annovst),
    show_rownames=FALSE,show_colnames=FALSE, scale="none",fontsize_row=5, annotation_colors=pbat_col,,  
    filename=file.path(analysis.dir,"Heatmap_20000mvPromoterMeth_noScale.pdf"))
#1000 mv
topVarGenesRld<- head(order(rowVars(PromoterMeth), decreasing=TRUE),1000)
anno <- Promoter
anno<-as.data.frame(anno)
annotop <- anno[topVarGenesRld,]
matrld <- PromoterMeth[topVarGenesRld,]
pheatmap(matrld,  annotation_col=as.data.frame(annovst), annotation_colors=pbat_col,,
    show_rownames=FALSE,show_colnames=FALSE, scale="row",fontsize_row=5, 
    filename=file.path(analysis.dir,"Heatmap_1000mvPromoterMeth.pdf"))
pheatmap(matrld,  annotation_col=as.data.frame(annovst), annotation_colors=pbat_col,,
    show_rownames=FALSE,show_colnames=FALSE, scale="none",fontsize_row=5,  
    filename=file.path(analysis.dir,"Heatmap_1000mvPromoterMeth_noScale.pdf"))


#PCA mv Promoter
#meth_per_cpg <- as.data.frame(bsseq::getMeth(bsseq_all, type = "raw"))
topVarGenesRld<- head(order(rowVars(PromoterMeth), decreasing=TRUE),20000)
anno <- Promoter
anno<-as.data.frame(anno)
annotop <- anno[topVarGenesRld,]
matrld <- PromoterMeth[topVarGenesRld,]
ir.pca <- prcomp(t(matrld),
                 center = TRUE,
                 scale. = TRUE) 
summary(ir.pca)
x <- ir.pca$x
x<- as.data.frame(cbind(x , pheno))
pdf(file.path(analysis.dir, "PC12_200000Promoter_subsetted.pdf"),height = 5, width = 5)
ggscatter(x, x="PC1", y="PC2",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC23_200000Promoter_subsetted.pdf"),height = 5, width = 5)
ggscatter(x, x="PC2", y="PC3",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC34_200000Promoter_subsetted.pdf"),height = 5, width = 5)
ggscatter(x, x="PC3", y="PC4",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()


#PCA mv CpG
#meth_per_cpg <- as.data.frame(bsseq::getMeth(bsseq_all, type = "raw"))
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),200000)
ir.pca <- prcomp(t(meth_per_cpg[topVarCpGs,]),
                 center = TRUE,
                 scale. = TRUE) 
summary(ir.pca)
x <- ir.pca$x
x<- as.data.frame(cbind(x , pheno))
pdf(file.path(analysis.dir, "PC12_200000CpGs_subsetted.pdf"),height = 5, width = 5)
ggscatter(x, x="PC1", y="PC2",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), 
          "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC23_200000CpGs_subsetted.pdf"),height = 5, width = 5)
ggscatter(x, x="PC2", y="PC3",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), 
          "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()
pdf(file.path(analysis.dir, "PC34_200000CpGs_subsetted.pdf"),height = 5, width = 5)
ggscatter(x, x="PC3", y="PC4",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
          "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
           theme(legend.position="right",legend.title = element_text(, size=10, 
                                      face="bold"))
dev.off()


#Sample Clustering
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),200000)
dissimilarity <- 1 - cor(meth_per_cpg[topVarCpGs,], use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)

library(dendextend)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
#col1 <- hrld$labels
#names(col1)<- anno[col1,]$group
#col1 <- col1[order.dendrogram(dend)]
#col1 <- ifelse(names(col1)=="TAM_tumor", "#EFC00099", ifelse(names(col1)=="BMDM_tumor","#0073C299",  ifelse(names(col1)=="TAM_healthy","#CD534CFF","#86868699")))
dend <- dend %>% 
#set("branches_k_color", k=4, c("#EFC00099", "#0073C299", "#CD534CFF","#86868699" )) %>%
#set("branches_col",col1) %>%
#set("branches_lwd", 2) %>%
#set("labels_colors",col1) %>% 
#set("labels_cex", .6 )%>%
#set("leaves_pch", 19)%>% 
#set("leaves_cex", 1.5)%>% 
#set("leaves_col", col1)
pdf(file.path(analysis.dir, "Clustering_correlation_200000CpGs.pdf"), height = 5, width = 5)
dend %>% plot
dev.off()

#export bigwig
meth_per_cpg_all <- as.data.frame(bsseq::getMeth(bsseq_all, type = "raw"))
library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
anno <- granges(bsseq_all)
dir.create(file.path(analysis.dir,"bw"))
for(i in rownames(pheno)){
    temp<- anno 
    score(temp) <- meth_per_cpg_all[,i]
    temp<- temp[!is.na(score(temp)),]
    seqlengths(temp)<-seqlengths(hg19)[names(seqlengths(temp))]
    export.bw(temp, file.path(analysis.dir,"bw",paste0(i, ".bw")))
}

cov_per_cpg_all <- as.data.frame(bsseq::getCoverage(bsseq_all))
dir.create(file.path(analysis.dir,"bw_cov"))
for(i in rownames(pheno)){
    temp<- anno 
    score(temp) <- cov_per_cpg_all[,i]
    temp<- temp[!is.na(score(temp)),]
    seqlengths(temp)<-seqlengths(hg19)[names(seqlengths(temp))]
    export.bw(temp, file.path(analysis.dir,"bw_cov",paste0(i, ".bw")))
    print(i)
}







#####not adapted yet
#violin plots of average methylation
## extract the methylation values
library(reshape2)
pheno <- pData(bsseq_all)
pheno$Origin <- as.character(pheno$Origin)
meth_per_cpg <- bsseq::getMeth(bsseq_all, type = "raw")
meth_per_cpg <- as.data.frame(meth_per_cpg)
met_per_cpg_group <- data.frame(BMDM_healthy=rowMeans(meth_per_cpg[, rownames(pheno[pheno$group =="BMDM_healthy",])], na.rm=TRUE),
    BMDM_tumor=rowMeans(meth_per_cpg[, rownames(pheno[pheno$group =="BMDM_tumor",])], na.rm=TRUE),
    TAM_tumor=rowMeans(meth_per_cpg[, rownames(pheno[pheno$group =="TAM_tumor",])], na.rm=TRUE),
    TAM_healthy=rowMeans(meth_per_cpg[, rownames(pheno[pheno$group =="TAM_healthy",])], na.rm=TRUE)
 )
met_per_cpg_group_melt <- reshape2::melt(met_per_cpg_group)
pdf(file.path(dir=analysis.dir, "AvMeth_cpg_violing.pdf"), height=4, width=4)
ggviolin(met_per_cpg_group_melt, x="variable", y="value", fill="variable",ylab="Average Methylation per CpG", ylim=c(0,1), palette =c("#86868699", "#0073C299", "#CD534CFF","#EFC00099" )) +rremove("legend") +rremove("xlab")+
rotate_x_text(angle = 45)
dev.off()

#ß value distribution
# from wide to long format
pdf(file.path(dir=analysis.dir, "ßvalue2.pdf")) 
gghistogram(met_per_cpg_group_melt, "value", color="variable",add_density=TRUE, palette =c("#86868699", "#0073C299","#EFC00099" , "#CD534CFF"))
dev.off()

#Clustering
meth_per_cpg <- as.data.frame(bsseq::getMeth(bsseq_all, type = "raw"))
topVarCpGs<- head(order(rowVars(as.matrix(meth_per_cpg)), decreasing=TRUE),200000)
dissimilarity <- 1 - cor(meth_per_cpg[topVarCpGs,], use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
library(rafalib)
pdf(file.path(analysis.dir, "Correlation_100000MVCPGs_Region_subsetted.pdf"), height = 5, width = 5)
myplclust(hclust(distance), labels=pheno$group, lab.col= as.fumeric(pheno$group))
dev.off()

scale_rows <- function(x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}
dissimilarity <- 1 - cor(scale_rows(meth_per_cpg[topVarCpGs,]), use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
pdf(file.path(analysis.dir, "Correlation_100000MVCPGs_Region_scaled_subsetted.pdf"), height = 5, width = 5)
myplclust(hclust(distance), labels=pheno$group, lab.col= as.fumeric(pheno$group))
dev.off()

#Sample Clustering2
#c("#A7303099", "#7AA6DC99")
dend<- dist(t(meth_per_cpg[topVarCpGs,]), method="manhattan")
pdf(file.path(analysis.dir, "Manhattan_100000MVCPGs_Region_scaled_subsetted.pdf"), height = 5, width = 5)
plot(hclust(dend))
dev.off()
