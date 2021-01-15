#libraries
library(methrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
#directories
odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
output.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
dir.create(output.dir,recursive=TRUE)

#load reference CpGs
hg19_cpgs = methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
bdg_PBAT_files = dir(path = odcf.dir , recursive=TRUE,pattern = "chr.bedgraph$", full.names = TRUE)

#prepare sample sheet
#PBAT
temp<-strsplit(bdg_PBAT_files, "_", fixed=TRUE)
patient<- sapply(temp , function(x)x[9])
temp <- sapply(temp , function(x)x[7])
temp <- sapply(strsplit(temp, "/", fixed=TRUE), "[",2)
tumor <-  sapply(strsplit(temp, "-", fixed=TRUE), "[",1)
replicate <-  sapply(strsplit(temp, "-", fixed=TRUE), "[",2)
replicate[!is.na(replicate)]<- paste0("_", replicate[!is.na(replicate)])
SampleID <- paste0(tumor,"_JMMLC_", patient, replicate)
SampleID <- gsub("NA", "", SampleID)
sample_anno_PBAT <- data.frame(SampleID=SampleID, 
    patient=patient, tumor=tumor)
#load marks sample_sheet
sample_anno_mark <- read.table(file.path(output.dir, "sample_anno_mark.txt"), stringsAsFactors=FALSE, header=TRUE)
sample_anno_mark$SampleID <- sample_anno_mark$Name
sample_anno<- left_join(sample_anno_PBAT,sample_anno_mark)
rownames(sample_anno)<- sample_anno$SampleID
#combine sample anno as well as bdg files
#read begraphs
meth = methrix::read_bedgraphs(files = bdg_PBAT_files, ref_cpgs = hg19_cpgs, 
    chr_idx = 1, start_idx = 2, beta_idx=4 ,M_idx = 5, U_idx = 6,
    stranded = TRUE, collapse_strands = TRUE,
    coldata=sample_anno,
    n_threads=3#, h5=TRUE, h5_dir=file.path(output.dir, "190923_methrix")
    )
dir.create(file.path(output.dir ,"methrix"))
saveRDS(meth, file.path(output.dir ,"methrix",  "methrix.rds"))


#QC report
meth <- readRDS( file.path(output.dir ,"methrix",  "methrix.rds"))
dir.create(file.path(output.dir, "methrix",  "QC_report_noFilter"))
methrix::methrix_report(meth = meth, recal_stats=TRUE,output_dir = file.path(output.dir, "methrix",  "QC_report_noFilter"))

#filter snps
ranges <-  get_matrix(meth, type = "M", add_loci = TRUE, in_granges = FALSE)
ranges_sub <- ranges[ranges$ch %in% paste0("chr",c(1:22, "X","Y")),]
ranges_sub<-ranges_sub[,c(1,2)]
ranges_sub$end <- ranges_sub$start+1
meth <- subset_methrix(meth, regions=ranges_sub)
meth_snpfil <- remove_snps(meth)
saveRDS(meth_snpfil, file.path(output.dir ,"methrix",  "methrix_snpRemoved.rds"))
#run QC report
methrix::methrix_report(meth = meth_snpfil, recal_stats=TRUE,output_dir = file.path(output.dir, "methrix",  "QC_report_snpRemoved"))

# 1.2: mask low coverage CpGs
meth_fil <- methrix::mask_methrix(m = meth_fil, low_count=3, high_quantile = 0.99)
saveRDS(meth_fil, file.path(output.dir ,"methrix",  "methrix_masked.rds"))

# 1.3: filter out uncovered CpGs
meth_fil <- methrix::remove_uncovered(m = meth_fil)
saveRDS(meth_fil, file.path(output.dir , "methrix_fil.rds"))

#QC report
dir.create(file.path(output.dir, "methrix", "QC_report_Filter"))
methrix::methrix_report(meth = meth_fil, output_dir = file.path(output.dir, "methrix", "QC_report_Filter"),recal_stats=TRUE)

#export as bsseq
bsseq_all <-methrix::methrix2bsseq(meth)
bsseq_all_sub <-methrix::methrix2bsseq(meth_fil)
bsseq_all_snpfil <-methrix::methrix2bsseq(meth_snpfil)

dir.create(file.path(output.dir, "bsseq"))
saveRDS(bsseq_all_sub, file.path(output.dir , "bsseq", "bsseq_all_fil.rds"))
saveRDS(bsseq_all, file.path(output.dir , "bsseq", "bsseq_all.rds"))
saveRDS(bsseq_all_snpfil, file.path(output.dir , "bsseq", "bsseq_all_snpfil.rds"))

