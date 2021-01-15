library(methrix)
library(ggpubr)
load("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/JMMLC_methrix_obj.Rda")
data.dir <- "/home/heyj/c010-datasets/Internal/mitotic_clock/data/"

M <- methrix::get_matrix(jmmlc_meth, type="M", add_loci=F)
pheno <- jmmlc_meth@colData
#mitotic clock on jmmml samples

#load cpgs
cpg_human <- read.table(file.path(data.dir, "cpgs_human.csv"), sep=",", header=TRUE)
cpg_human$Chromosome <- gsub("chr","", as.character(cpg_human$Chromosome))
cpg_human_gr<-makeGRangesFromDataFrame(cpg_human)

jmml_mc<- methrix::get_region_summary(jmmlc_meth, regions=cpg_human_gr, type="M",how="mean")
jmml_mc<- as.data.frame(jmml_mc)
jmml_mc <- jmml_mc[,5:ncol(jmml_mc)]
round(colMeans(jmml_mc, na.rm=TRUE), 3)

M <- jmml_mc


#get average methylation
#per methylation status
av.meth_per_cpg<- colMeans(M, na.rm = TRUE)

av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Methylation_Level=as.character(pheno$Methylation_Level))

compare_means(AverageMethylation ~ Methylation_Level,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_MethylationStatus.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Methylation_Level", y="AverageMethylation", color="Methylation_Level",ylim=c(0,1), palette =c("#86868699", "#0073C299" )
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_MethylationStatus2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Methylation_Level", y="AverageMethylation", color="Methylation_Level",ylim=c(0,.2), palette =c("red", "#0073C299" )
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.20))+stat_compare_means(label.y = .001) +rremove("legend")
dev.off()



#per patient
av.meth_per_cpg<- colMeans(M, na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Patient_ID=as.character(pheno$Patient_ID))

compare_means(AverageMethylation ~ Patient_ID,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("D117", "D123"),c("D117", "D124"),c("D123", "D124"))

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_Patient_ID.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient_ID", y="AverageMethylation", color="Patient_ID",ylim=c(0,1), palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_Patient_ID2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient_ID", y="AverageMethylation", color="Patient_ID",ylim=c(0,.2), palette =c( "#0073C299", "red","red" )
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.17, .18,.19))+stat_compare_means(label.y = .0) +rremove("legend")
dev.off()
#per tumor section
av.meth_per_cpg<- colMeans(M, na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Tumor_Section=as.character(pheno$Tumor_Section))

compare_means(AverageMethylation ~ Tumor_Section,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_Tumor_Section.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="Tumor_Section",ylim=c(0,1), palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_Tumor_Section2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="Tumor_Section",ylim=c(.7,.9), palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMethMitoticClock_allCpGs_Tumor_Section2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="gray44",ylim=c(0,.2), palette =c("gray44")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .18) +rremove("legend")
dev.off()




