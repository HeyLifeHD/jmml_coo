#Combine methrix objects of hierachy and jmml

#libraries
library(methrix)
library(bsseq)
#directories
output.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"

#load jmml dataset
meth_jmml <- readRDS(file.path("/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/" ,"methrix",  "methrix_snpRemoved.rds"))
#do changes that were previously done in bsseq object
meth_jmml <- meth_jmml[,!colnames(meth_jmml) %in% c("tumor11_JMMLC_D129_1","tumor11_JMMLC_D117", "tumor00_JMMLC_D123", "tumor10_JMMLC_D123",
    "tumor01_JMMLC_D123", "tumor00_JMMLC_D124", "tumor01_JMMLC_D117")] 
dim(meth_jmml)
pheno_jmml<- colData(meth_jmml)
pheno_jmml$patient <- droplevels(pheno_jmml$patient)
pheno_jmml$Epigenotype <- as.factor(pheno_jmml$Epigenotype )
pheno_jmml$Epigenotype <-factor(pheno_jmml$Epigenotype, levels = c("LM", "HM"))
pheno_jmml$Genotype <- as.factor(pheno_jmml$Genotype)
pheno_jmml$Genotype <-factor(pheno_jmml$Genotype, levels = c("neg", "KRAS", "PTPN11"))
#add column for cell type
pheno_jmml$CellType <- "JMML"
colData(meth_jmml) <- pheno_jmml

#load hierachy dataset
meth_hierachy <- readRDS(file.path("/home/heyj/c010-projects/Reka/39_scMeth_HSC/JMML_all_HSC_data_scripts/" , "raw_methrix.RDS"))
write.table(colData(meth_hierachy),file.path(output.dir,"methrix", "sample_anno_hierachy_all.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#do subset 
pheno_hierachy<- colData(meth_hierachy)
pheno_hierachy_sub <- pheno_hierachy[pheno_hierachy$Sample_Type=="normal",]
write.table(pheno_hierachy_sub,file.path(output.dir,"methrix", "sample_anno_hierachy_NORMAL.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

meth_hierachy_sub <- meth_hierachy[,rownames(pheno_hierachy_sub)]
#add column for cell type
pheno_hierachy_sub$CellType <- pheno_hierachy_sub$Celltype
pheno_hierachy_sub$X <- NULL
colData(meth_hierachy_sub) <- pheno_hierachy_sub
#filter snps
ranges <-  get_matrix(meth_hierachy_sub, type = "M", add_loci = TRUE, in_granges = FALSE)
ranges_sub <- ranges[ranges$ch %in% paste0("chr",c(1:22, "X","Y")),]
ranges_sub<-ranges_sub[,c(1,2)]
ranges_sub$end <- ranges_sub$start+1
meth_hierachy_sub <- subset_methrix(meth_hierachy_sub, regions=ranges_sub)
meth_hierachy_sub_snpfil <- remove_snps(meth_hierachy_sub)

#get similar pheno pheno objects
sample_anno_combined <- read.table(file.path(output.dir,"methrix", "sample_anno_HSC_comb.txt"), sep="\t", header=TRUE)
row.names(sample_anno_combined) <- sample_anno_combined$Name
sample_anno_hierachy <- sample_anno_combined[sample_anno_combined$Patient=="normal",]
colData(meth_hierachy_sub_snpfil) <- DataFrame(sample_anno_hierachy)
sample_anno_jmml <- sample_anno_combined[!sample_anno_combined$Patient=="normal",]
sample_anno_jmml <- sample_anno_jmml[sample_anno_jmml$Remove ==FALSE,]
sample_anno_jmml <- sample_anno_jmml[colnames(meth_jmml),]
colData(meth_jmml)<- DataFrame(sample_anno_jmml)

#combine
meth_comb <- cbind(meth_jmml, meth_hierachy_sub_snpfil)
write.table(colData(meth_comb),file.path(output.dir,"methrix", "sample_anno_HSC_comb_finalOrder.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#save methrix
saveRDS(meth_comb, file.path(output.dir,"methrix", "methrix_HSC_comb_snpRemoved.rds"))
meth_comb <- readRDS(file.path(output.dir,"methrix", "methrix_HSC_comb_snpRemoved.rds"))
#run QC report
methrix::methrix_report(meth = meth_comb, recal_stats=TRUE,output_dir = file.path(output.dir, "methrix",  "QC_report_HSC_comb_snpRemoved"), n_thr=10)

#export as bsseq
rm(meth_jmml,meth_hierachy_sub, meth_hierachy )
bsseq_all <-methrix::methrix2bsseq(meth_comb)
saveRDS(bsseq_all, file.path(output.dir , "bsseq", "bsseq_HSC_comb_snpRemoved.rds"))
