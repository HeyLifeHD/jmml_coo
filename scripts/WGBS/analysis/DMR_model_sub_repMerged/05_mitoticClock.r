##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
#look at mitotic clock

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
data.dir <- "c010-datasets/Internal/mitotic_clock/data/"
output.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"

#load data
bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#change annotation
pheno <- pData(bsseq_all)
pheno$Tissue<- c( rep("adult_bonemarrow", 15))
pheno$Donor <- c(as.character(pheno$Patient[1:15])) 
pheno$Epigenotype <- as.character(pheno$Epigenotype)
pheno[pheno$Patient %in% c("I217", "D217"),]$Epigenotype <- "IM"
pheno$Epigenotype <- as.factor(pheno$Epigenotype)
pheno$Celltype <- c("MPP","MPP","MPP","MPP" ,"MPP","HSC", "HSC", "HSC","LMPP","LMPP","LMPP","CD45RACD90","CD45RACD90","CD45RACD90","CD45RACD90")
pheno$Sample_Type <- "tumor"
pheno$Tumor<- pheno$tumor
pData(bsseq_all) <- pheno

#new output directory
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_model_sub_repMerged"
dir.create(analysis.dir)

#load cpgs
cpg_human <- read.table(file.path(data.dir, "cpgs_human.csv"), sep=",", header=TRUE)
#cpg_human$Chromosome <- gsub("chr","", as.character(cpg_human$Chromosome))
cpg_human_gr<-makeGRangesFromDataFrame(cpg_human)

jmml_mc <- bsseq::getMeth(bsseq_all, regions=cpg_human_gr, what="perRegion", type="raw")
round(colMeans(jmml_mc, na.rm=TRUE), 3)
M <- jmml_mc


#get average methylation
#per methylation status
av.meth_per_cpg <- as.data.frame(colMeans(M, na.rm = TRUE))
colnames(av.meth_per_cpg)<- "AverageMethylation"
pheno <- pData(bsseq_all)
av.meth_per_cpg <- as.data.frame(cbind(av.meth_per_cpg, pheno))

compare_means(AverageMethylation ~ Epigenotype,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("HM", "LM"), c("HM", "IM"), c("LM", "IM") )
med <- aggregate(av.meth_per_cpg[,"AverageMethylation",drop=FALSE], list(av.meth_per_cpg$Epigenotype), median)
ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
dir.create(file.path(analysis.dir,"mitotic_clock"))
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_MethylationStatus.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,1), palette = c(wildtype ="#ababab", LM = "#0058b4", IM = "#fbbb25", HM = "#c33126")
,add = "jitter", order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(0.5,0.6,.70))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_MethylationStatus2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,.4), palette = c(wildtype ="#ababab", LM = "#0058b4", IM = "#fbbb25", HM = "#c33126")
,add = "jitter", order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.20,.3,.4))+stat_compare_means(label.y = .001) +rremove("legend")
dev.off()

#per patient
med <- aggregate(av.meth_per_cpg[,"AverageMethylation",drop=FALSE], list(av.meth_per_cpg$Patient), median)
ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_patient.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="patient", y="AverageMethylation", color="patient",ylim=c(0,1), palette =c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220")

,add = "jitter", order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_patient2.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="patient", y="AverageMethylation", color="patient",ylim=c(0,.2), palette =c(cordblood = "#737373", adult_bonemarrow ="#ababab", 
                    D117 = "#0058b4", D129 = "#2188c9", 
                    D217 = "#fbbb25", I217 = "#fca349", 
                    D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220")
,add = "jitter", order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.17, .18,.19))+stat_compare_means(label.y = .0) +rremove("legend")
dev.off()

#per tumor
compare_means(AverageMethylation ~ tumor,  data = av.meth_per_cpg, method = "t.test")
med <- aggregate(av.meth_per_cpg[,"AverageMethylation",drop=FALSE], list(av.meth_per_cpg$tumor), median)
ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_tumor.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="tumor", y="AverageMethylation", color="tumor",ylim=c(0,1), palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12", 
                       tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#99a637")
,add = "jitter",order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_tumor2.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="tumor", y="AverageMethylation", color="tumor",ylim=c(0,.2), palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12", 
                       tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#99a637")
,add = "jitter", order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

#per Celltype
compare_means(AverageMethylation ~ Celltype,  data = av.meth_per_cpg, method = "t.test")
med <- aggregate(av.meth_per_cpg[,"AverageMethylation",drop=FALSE], list(av.meth_per_cpg$Celltype), median)
ord <- med[order(med$AverageMethylation, decreasing=FALSE),]$Group.1
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_Celltype.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="Celltype", y="AverageMethylation", color="tumor",ylim=c(0,1), palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12", 
                       tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#99a637")
,add = "jitter",order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_Celltype2.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="Celltype", y="AverageMethylation", color="tumor",ylim=c(0,.2), palette =c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12", 
                       tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#99a637")
,add = "jitter", order=ord
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
