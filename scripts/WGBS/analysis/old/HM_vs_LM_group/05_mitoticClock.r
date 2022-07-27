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
odcf.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

data.dir <- "c010-datasets/Internal/mitotic_clock/data/"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub_cov.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

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
my_comparisons <- list( c("HM", "LM") )

dir.create(file.path(analysis.dir,"mitotic_clock"))
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_MethylationStatus.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,1), palette = c(HM = "#b5241c", LM = "#0041a5")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_MethylationStatus2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Epigenotype", y="AverageMethylation", color="Epigenotype",ylim=c(0,.2), palette = c(HM = "#b5241c", LM = "#0041a5")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.20))+stat_compare_means(label.y = .001) +rremove("legend")
dev.off()



#per patient
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_patient.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="patient", y="AverageMethylation", color="patient",ylim=c(0,1), palette =c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000")

,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_patient2.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="patient", y="AverageMethylation", color="patient",ylim=c(0,.2), palette =c(D117 = "#001885", D129 = "#305382", D217 = "#77694F", I217 = "#B18637", 
                      D213 = "#ea7e23", D360 = "#dd5b1c", D124 = "#e2321b", D123 = "#700000")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.17, .18,.19))+stat_compare_means(label.y = .0) +rremove("legend")
dev.off()
#per tumor
pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_tumor.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="tumor", y="AverageMethylation", color="tumor",ylim=c(0,1), palette =c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path(analysis.dir,"mitotic_clock","AvMethMitoticClock_allCpGs_tumor2.pdf"), height=4, width=7)
ggpubr::ggboxplot(av.meth_per_cpg, x="tumor", y="AverageMethylation", color="tumor",ylim=c(0,.2), palette =c(tumor01 = "#252525", tumor00 = "#737373", tumor10 = "#9babcf", tumor11 = "#d4b38d")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
