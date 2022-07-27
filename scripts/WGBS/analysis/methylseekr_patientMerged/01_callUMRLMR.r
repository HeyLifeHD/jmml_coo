#Libraries
library(MethylSeekR)
library(parallel)
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPseeker)
library(ggpubr)
#load my own data
#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
methylseekr.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200811_mergedPatientsHSCcb_methylseekr"
dir.create(analysis.dir)

#load data
#bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved.rds"))#
#collapse patient replicates and normal celltypes
#pheno <- pData(bsseq_all)
#pheno$group <- paste0(pheno$Tissue_Source, "_", as.character(unlist(pheno$Celltype)))
#pheno[pheno$Sample_Type=="tumor", ]$group <- as.character(pheno[pheno$Sample_Type=="tumor", ]$Patient)
#bsseq_merged <- collapseBSseq(bsseq_all, pheno$group)
#saveRDS(bsseq_merged, file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_PatientMerged_CelltypeMerged.rds"))
bsseq_merged <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_PatientMerged_CelltypeMerged.rds"))

#prepare methylation data in bedgraph format
anno <- granges(bsseq_merged)
pheno <- pData(bsseq_merged)
cov <- bsseq::getCoverage(bsseq_merged)
meth <- bsseq::getMeth(bsseq_merged, type="raw")
meth_reads <- meth*cov
meth.gr <- list()
for(i in rownames(pheno)){
    meth.gr[[i]] <- anno
    meth.gr[[i]]$T <- cov[,i]
    meth.gr[[i]]$M <- meth_reads[,i]
}

#calculate the distribution of α-values for one selected chromosome. α characterizes the distribution of methylation levels in sliding windows containing i00 consecutive CpGs along the genome.
for(i in names(meth.gr)){
    dir.create(file.path(methylseekr.dir, i))
    plotAlphaDistributionOneChr(m= meth.gr[[i]], chr.sel="chr19", num.cores=6, pdfFilename=file.path(methylseekr.dir, i, paste0("alphaDistribution_chr19_",i, ".pdf")))
    print(i)
    print(paste0(which(names(meth.gr)==i), " of ", length(names(meth.gr))))
}
for(i in names(meth.gr)){
    dir.create(file.path(methylseekr.dir, i))
    plotAlphaDistributionOneChr(m= meth.gr[[i]], chr.sel="chr2", num.cores=6, pdfFilename=file.path(methylseekr.dir, i, paste0("alphaDistribution_chr2_",i, ".pdf")))
    print(i)
    print(paste0(which(names(meth.gr)==i), " of ", length(names(meth.gr))))
}


#Hidden Markov Model (HMM) to identify PMDs genome-wide 
#--> not needed
#PMDsegments.gr <- list()
#for(i in names(meth.gr)){
#    PMDsegments.gr[[i]] <- segmentPMDs(m=meth.gr[[i]], chr.sel="chr19", seqLengths=sLengths, num.cores=10,  
#        pdfFilename=file.path(methylseekr.dir,i, paste0("PMDsegmentationHisto_chr19_",i, ".pdf")))
#    plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr[[i]], PMDsegments.gr[[i]][values(PMDsegments.gr[[i]])$type=="notPMD"]), 
#        chr.sel="chr19", num.cores=10, 
#        pdfFilename=file.path(methylseekr.dir, i, paste0("PMDsegmentationHisto_PMDremoved_chr19_",i, ".pdf")))
#    plotPMDSegmentation(m=meth.gr[[i]], segs=PMDsegments.gr[[i]], 
#        pdfFilename=file.path(methylseekr.dir, i, paste0("PMDsegmentation_chr19_",i, ".pdf")))
#    print(i)
#    print(paste0(which(names(meth.gr)==i), " of ", length(names(meth.gr))))#
#
#}

#Identification of UMRs and LMRs
library(rtracklayer)
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "cpgIslandExt")
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))

#calculate fdr rates for parameter determination
stats <- list()
for(i in names(meth.gr)){
    stats[[i]] <- calculateFDRs(m=meth.gr[[i]], CGIs=CpGislands.gr, num.cores=5,  pdfFilename=file.path(methylseekr.dir, i, paste0("FDRcalculation_", i, ".pdf")))
    print(i)
    print(paste0(which(names(meth.gr)==i), " of ", length(names(meth.gr))))
}

#select parameters
FDR.cutoff <- 5
m.sel <- 0.5
n.sel <- list()
for(i in names(stats)){
    n.sel[[i]]=as.integer(names(stats[[i]]$FDRs[as.character(m.sel), ][stats[[i]]$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
}
n.sel[[i]]

#identify UMRs and LMRs
#extract seqlengths
sLengths=seqlengths(Hsapiens)

UMRLMRsegments.gr <- list()
for(i in names(meth.gr)){
    UMRLMRsegments.gr[[i]] <- segmentUMRsLMRs(m=meth.gr[[i]], meth.cutoff=m.sel, nCpG.cutoff=n.sel[[i]],
        num.cores=5, myGenomeSeq=Hsapiens, seqLengths=sLengths,
        pdfFilename=file.path(methylseekr.dir, i, paste0("UMRLMRSegmentHeat_",i, ".pdf")))
    saveUMRLMRSegments(segs=UMRLMRsegments.gr[[i]],
        GRangesFilename=file.path(methylseekr.dir,i, paste0("UMRsLMRs",i, ".rds")), 
        TableFilename=file.path(methylseekr.dir,i, paste0("UMRsLMRs",i, ".tab")))
    plotFinalSegmentation(m=meth.gr[[i]], segs=UMRLMRsegments.gr[[i]], meth.cutoff=m.sel, 
        pdfFilename=file.path(methylseekr.dir,i, paste0("FinalSegmentation_",i, ".pdf")))
    print(i)
    print(paste0(which(names(meth.gr)==i), " of ", length(names(meth.gr))))

}
saveRDS(UMRLMRsegments.gr,file.path(methylseekr.dir,  paste0("UMRsLMRs", ".rds")))
head(UMRLMRsegments.gr[[i]])
lapply(UMRLMRsegments.gr, length)
