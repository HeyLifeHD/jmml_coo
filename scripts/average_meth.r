library(methrix)
library(ggpubr)
load("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/JMMLC_methrix_obj.Rda")

M <- methrix::get_matrix(jmmlc_meth, type="M", add_loci=F)
M <- as.data.frame(M)
pheno <- colData(jmmlc_meth)


#get average methylation
#per methylation status
av.meth_per_cpg<- colMeans(M, na.rm = TRUE)

av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Methylation_Level=as.character(pheno$Methylation_Level))

compare_means(AverageMethylation ~ Methylation_Level,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_MethylationStatus.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Methylation_Level", y="AverageMethylation", color="Methylation_Level",ylim=c(0,1), palette =c("#86868699", "#0073C299" )
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_MethylationStatus2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Methylation_Level", y="AverageMethylation", color="Methylation_Level",ylim=c(.7,.9), palette =c("#86868699", "#0073C299" )
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()



#per patient
av.meth_per_cpg<- colMeans(M, na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Patient_ID=as.character(pheno$Patient_ID))

compare_means(AverageMethylation ~ Patient_ID,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_Patient_ID.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient_ID", y="AverageMethylation", color="Patient_ID",ylim=c(0,1), palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_Patient_ID2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient_ID", y="AverageMethylation", color="Patient_ID",ylim=c(.7,.9), palette =c("#0073C299", "red", "red")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
#per tumor section
av.meth_per_cpg<- colMeans(M, na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Tumor_Section=as.character(pheno$Tumor_Section))

compare_means(AverageMethylation ~ Tumor_Section,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_Tumor_Section.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="Tumor_Section",ylim=c(0,1), palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_Tumor_Section2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="Tumor_Section",ylim=c(.7,.9), palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_allCpGs_Tumor_Section2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="gray44",ylim=c(.7,.9), palette =c("gray44")
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()





#Plot average Methylation of CpGs islands
#cpg locations
library(RnBeads)
rnb.options(assembly="hg19")
cgi <- rnb.annotation2data.frame(rnb.get.annotation("cpgislands"))
chrm<-as.vector(gsub("chr", "", cgi$Chromosome))

cgi$Chromosome <-NULL
cgi$Chromosome <- chrm
cgi <-makeGRangesFromDataFrame(cgi)
av.meth_per_cpg_island <- methrix::get_region_summary(jmmlc_meth, regions=cgi, type="M",how="mean")
av.meth_per_cpg_island <- as.data.frame(av.meth_per_cpg_island)
av.meth_per_cpg_island <- av.meth_per_cpg_island [, 5:ncol( av.meth_per_cpg_island )]

av.meth_per_cpg<- colMeans(av.meth_per_cpg_island , na.rm = TRUE)

av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Methylation_Level=as.character(pheno$Methylation_Level), Patient_ID=pheno$Patient_ID)

compare_means(AverageMethylation ~ Methylation_Level,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_MethylationStatus.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Methylation_Level", y="AverageMethylation", color="Methylation_Level",ylim=c(0,1), palette =c("red", "#0073C299" ), legend="right"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.5))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_MethylationStatus2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Methylation_Level", y="AverageMethylation", color="Methylation_Level", palette =c("red", "#0073C299" )
,add = "jitter", ylim=c(0,.4)
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
#per patient

av.meth_per_cpg_island <- methrix::get_region_summary(jmmlc_meth, regions=cgi, type="M",how="mean")
av.meth_per_cpg_island <- as.data.frame(av.meth_per_cpg_island)
av.meth_per_cpg_island <- av.meth_per_cpg_island [, 5:ncol( av.meth_per_cpg_island )]

av.meth_per_cpg<- colMeans(av.meth_per_cpg_island , na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Patient_ID=as.character(pheno$Patient_ID))

compare_means(AverageMethylation ~ Patient_ID,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_Patient_ID.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient_ID", y="AverageMethylation", color="Patient_ID",ylim=c(0,1), palette =c("#0073C299", "red", "red"),
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_Patient_ID2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Patient_ID", y="AverageMethylation", color="Patient_ID",  palette =c("#0073C299", "red", "red"),
,add = "jitter", ylim=c(0,.5)
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()

#per patient
av.meth_per_cpg_island <- methrix::get_region_summary(jmmlc_meth, regions=cgi, type="M",how="mean")
av.meth_per_cpg_island <- as.data.frame(av.meth_per_cpg_island)
av.meth_per_cpg_island <- av.meth_per_cpg_island [, 5:ncol( av.meth_per_cpg_island )]

av.meth_per_cpg<- colMeans(av.meth_per_cpg_island , na.rm = TRUE)
av.meth_per_cpg <- data.frame(AverageMethylation=av.meth_per_cpg, Tumor_Section=as.character(pheno$Tumor_Section), Patient_ID=as.character(pheno$Patient_ID))

compare_means(AverageMethylation ~ Tumor_Section,  data = av.meth_per_cpg, method = "t.test")
my_comparisons <- list( c("Hyper", "Low") )

dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_Tumor_Section.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="Tumor_Section", palette ="jco"
,add = "jitter"
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
dir.create(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka"))
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_Tumor_Section2.pdf"), height=4, width=4)
ggpubr::ggboxplot(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="gray44",
,add = "jitter", ylim=c(0.25,.35)
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) +rremove("legend")
dev.off()
pdf(file.path("c010-datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/05_results/03_Methirx/", "temp_analysis_Joschka","AvMeth_CpGislands_Tumor_Section_scatter.pdf"), height=4, width=4)
ggpubr::ggscatter(av.meth_per_cpg, x="Tumor_Section", y="AverageMethylation", color="Patient_ID",ylim=c(0.2,.4), palette ="jco" 
)  +  stat_compare_means(method="t.test", comparisons = my_comparisons,label.y = c(.90))+stat_compare_means(label.y = .25) 
dev.off()

