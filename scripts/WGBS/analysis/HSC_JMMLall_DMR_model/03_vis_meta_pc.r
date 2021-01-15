
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

#pca of main comparison of interestt
meth_dmr <- mcols(dmrs_final[["Sample_Typetumor_hierachy_DMRs_removed"]])[,rownames(pData(bsseq_all))]
rownames(meth_dmr)<- paste0(1:nrow(meth_dmr),"_",dmrs_final[["Sample_Typetumor_hierachy_DMRs_removed"]]$SYMBOL )

meth_dmr <- meth_dmr[complete.cases(meth_dmr),]
ir.pca <- prcomp(t(as.matrix(meth_dmr)),
            center = T,
            scale. = F) 

print(i)
print(summary(ir.pca))
pc <- ir.pca$x
pc<- as.data.frame(cbind(pc , pData(bsseq_all)))



#extract ir.pca1 and ir.pca3
#pc1
rotation <- as.data.frame(ir.pca$rotation)
rot_pc1 <- rotation[order(rotation$PC1,decreasing=TRUE), ]
top100_ir.pc1 <- head(rot_pc1, 100)
paste0(unique(sapply(strsplit(rownames(top100_ir.pc1), "_"),"[",2)), collapse=", ")
last100_ir.pc1  <- tail(rot_pc1, 100)
paste0(unique(sapply(strsplit(rownames(last100_ir.pc1), "_"),"[",2)), collapse=", ")
rot_pc1_abs <- rotation[order(abs(rotation$PC1),decreasing=TRUE), ]
abs200_ir.pc1<- head(rot_pc1_abs, 200)
paste0(unique(sapply(strsplit(rownames(abs200_ir.pc1), "_"),"[",2)), collapse=", ")

rotation <- as.data.frame(ir.pca$rotation)
rot_pc1 <- rotation[order(rotation$PC1,decreasing=TRUE), ]
top100_ir.pc1 <- head(rot_pc1, 1000)
paste0(unique(sapply(strsplit(rownames(top100_ir.pc1), "_"),"[",2)), collapse=", ")
last100_ir.pc1  <- tail(rot_pc1, 1000)
paste0(unique(sapply(strsplit(rownames(last100_ir.pc1), "_"),"[",2)), collapse=", ")
rot_pc1_abs <- rotation[order(abs(rotation$PC1),decreasing=TRUE), ]
abs200_ir.pc1<- head(rot_pc1_abs, 2000)
paste0(unique(sapply(strsplit(rownames(abs200_ir.pc1), "_"),"[",2)), collapse=", ")



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
files <- list.files(file.path(analysis.dir,"Sample_Typetumor_hierachy_DMRs_removed",  "metascape_jmml_vs_normal_dmrs_pc1"), 
    pattern="metascape_result.xlsx", full.names=TRUE, recursive=TRUE)
paths <- gsub("metascape_result.xlsx", "", files)
files_hypo<- files[grep("top", files)]
paths_hypo<- paths[grep("top", paths)]
files_hyper<- files[grep("last", files)]
paths_hyper<- paths[grep("last", paths)]
files_abs<- files[grep("abs", files)]
paths_abs<- paths[grep("abs", paths)]
#plotting
for(numb in c(3,5,10,15)){
for(i in 1:length(files_hypo)){
    plotMeta(files_hypo[i], 
        file.path(paths_hypo[i], paste0("metascape_vis_top_",numb,"_1.pdf")), 
        n=numb, width=5, height=5)
}
#start from 2 because no enrichment for hyper
for(i in 1:length(files_hyper)){
    plotMeta(files_hyper[i], 
        file.path(paths_hyper[i], paste0("metascape_vis_last_",numb,"_1.pdf")), 
        n=numb, width=5, height=5,
        color= "#67A9CF")
}

for(i in 1:length(files_abs)){
    plotMeta(files_abs[i], 
        file.path(paths_abs[i], paste0("metascape_vis_abs_",numb,"_1.pdf")), 
        n=numb, width=5, height=5,
        color= "gray")
}
}
