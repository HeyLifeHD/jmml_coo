#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/jmmlc_rnaseq/220805_rnaseq_knownDEG"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)


#load data
counts <- readRDS("/omics/groups/OE0219/internal/jmmlc_rnaseq/220802_rnaseq_processing_known_hg19_starRSEM/star_salmon/salmon.merged.gene_counts.rds")
counts <- assay(counts)
colnames(counts)<- sapply(strsplit(colnames(counts), "_", fixed=TRUE), "[",2)
#create annotation sheet
sample_anno <- data.table::fread("/omics/groups/OE0219/internal/jmmlc_rnaseq/data/2020-09-29_Patient_Characteristics_Full_Meta.csv")
sample_anno <- sample_anno[sample_anno$consensusCluster3 %in% c("LM", "IM", "HM"),]
sample_anno$Epigenotype <- ifelse(sample_anno$consensusCluster3 =="HM", "HM", "non_HM")
sample_anno <- sample_anno[sample_anno$EWOG_ID %in% colnames(counts),]
rownames(sample_anno)<- sample_anno$EWOG_ID 
counts <- counts[,colnames(counts)%in% sample_anno$EWOG_ID ]
counts<- counts[,sample_anno$EWOG_ID ]

counts_int <- sapply(counts, as.integer)
rownames(counts_int)<- rownames(counts)
#create dds
dds <- DESeqDataSetFromMatrix(countData = counts_int, 
                              colData = sample_anno, 
                              design = ~ Epigenotype )

#subset samples of wrong strandedness
dds <- dds[,!colnames(counts)%in% c("D124", "D123", "D117")]
#Estimating size factors
#for Genotype comparison
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
table(idx)
#dds <- dds[idx,]
#dim(dds)

#Running the DESEQ
dds <- DESeq(dds)
saveRDS(dds, file.path(results.dir,"dds.rds"))