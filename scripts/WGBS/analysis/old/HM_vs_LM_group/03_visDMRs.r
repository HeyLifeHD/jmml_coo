##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Vis DMRs

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
#Directories
odcf.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub_cov.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

#select 
dmrs_final<- list()
dmrs_final$HM_vs_LM <-dmrs
#add methylation information of each samples
meth <- list()
for(i in names(dmrs_final)){
    meth[[i]] <- bsseq::getMeth(bsseq_all, regions=dmrs_final[[i]], what="perRegion", type="raw")
    mcols(dmrs_final[[i]]) <- cbind(mcols(dmrs_final[[i]]), meth[[i]]) 
}
saveRDS(dmrs_final, file.path(analysis.dir, "dmrs_gr_anno_meth.rds"))


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

    pdf(file.path(analysis.dir,i, "visualization", "PCA12_allDMRs.pdf"), height=4, width=4)
    print(ggscatter(pc, x="PC1", y="PC2",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), "% variance")), 
          ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance")))+theme(legend.position="right"))
    dev.off()

    pdf(file.path(analysis.dir, i, "visualization","PCA23_allDMRs.pdf"), height=4, width=4)
    print(ggscatter(pc, x="PC2", y="PC3",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance")), 
          ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance")))+theme(legend.position="right"))
    dev.off()

    pdf(file.path(analysis.dir, i, "visualization","PCA34_allDMRs.pdf"),height = 3.5, width = 5)
    print(ggscatter(pc, x="PC3", y="PC4",
          color = "Patient", shape = "tumor",#size="Protocol",
          ellipse = F , mean.point = FALSE,palette= c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"),
          star.plot = F, xlab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), 
            "% variance")), ylab=(paste0("PC4: ", round(summary(ir.pca)$importance[2,4]*100,2), "% variance"))) +
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

for(i in names(dmrs_final)){
    #all dmrs in each comparison with all samples
    meth_dmr <- mcols(dmrs_final[[i]])[,rownames(pData(bsseq_all))]
    rownames(meth_dmr)<- paste0("Dmr.",1:nrow(meth_dmr),"_", dmrs_final[[i]]$SYMBOL)
    pbat_col = list(
    Patient = c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000"), 
    Genotype = c(neg = "#319a38", KRAS = "#99a637", PTPN11 = "#007458"),
    Epigenotype = c(HM = "#b5241c", LM = "#0041a5"), 
    tumor = c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d"))
    annovst <- as.data.frame(colData(bsseq_all))[, c("Epigenotype", "Patient", "tumor", "Genotype")] 

    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        scale="row",fontsize_row=5,  annotation_color=pbat_col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_rowScale.pdf"))
    pheatmap(meth_dmr,  annotation_col=as.data.frame(annovst),show_rownames=FALSE,show_colnames=FALSE, 
        scale="none",fontsize_row=5,  annotation_color=pbat_col,
        filename=file.path(analysis.dir,i, "visualization","Heatmap_DMRs_noScale.pdf"))

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
dis <- dmrs_final[[comp]]$distanceToTSS
dis <- as.data.frame(dis)
dis <- dis/1000
dis <-dis[dis$dis > (-14) & dis$dis < (14),]
dis <- as.data.frame(dis)
#histo
pdf(file.path(analysis.dir,comp,"visualization",paste0("Distance_to_TSS_Histogram.pdf")), height=3.5)
print(gghistogram(dis, x="dis", fill="lightgray", bins=40, rug = T, xlab="Distance from TSS (kb)", ylab="# of DMRs", xlim=c(-10,10)
)  + geom_vline(xintercept=0, linetype = 1, color="red", size=2 ))
dev.off() 
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