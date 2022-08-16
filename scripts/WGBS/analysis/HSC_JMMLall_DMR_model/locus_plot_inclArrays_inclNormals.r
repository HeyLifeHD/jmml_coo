#Libraries
library(bsseq)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(Gviz)
#load my own data: Epigenotype vs HSC_cb
#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
input_external.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/DNAmethylationArrays"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"
#create directory
dir.create("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal")

#load meth smooth data
bsseq_all_smooth <- readRDS( "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/bsseq/bsseq_HSC_combTotal_snpRemoved_smooth.rds")
#subset bsseq
phenoTotal <- pData(bsseq_all_smooth)
drop <- phenoTotal$Name %in% c("FL2_HSC_PBAT_2", "JU3_HSC_PBAT_11")
bsseq_all_smooth <- bsseq_all_smooth[, !drop]
pheno <- phenoTotal[, !drop]
#load dmrs
dmrs_final<- readRDS(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged", "dmrs_gr_sub_MethDiff.rds"))
dmrs_normal<- readRDS(file.path( "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC", "dmrs_gr_sub_MethDiff.rds"))

#load anno
anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")

#load 450k
HM_array <- fread(file.path(input_external.dir, "HM_samples_nvCpGs_betas_2016-12-28.igv"),skip=1)
IM_array <- fread(file.path(input_external.dir, "IM_samples_nvCpGs_betas_2016-12-28.igv"),skip=1)
LM_array <- fread(file.path(input_external.dir, "LM_samples_nvCpGs_betas_2016-12-28.igv"),skip=1)

#prepare methylation data of 450k arrays
HM_array_gr <-  GRanges(
   seqnames = HM_array$V1,
   ranges = IRanges(start = HM_array$V2,
                  end =HM_array$V3))
score(HM_array_gr)<- HM_array$V5/100
IM_array_gr <-  GRanges(
   seqnames = IM_array$V1,
   ranges = IRanges(start = IM_array$V2,
                  end =IM_array$V3))
score(IM_array_gr)<- IM_array$V5/100
LM_array_gr <-  GRanges(
   seqnames = LM_array$V1,
   ranges = IRanges(start = LM_array$V2,
                  end =LM_array$V3))
score(LM_array_gr)<- LM_array$V5/100

#subset annos for symbol parsing
anno_original_transcript_symbol <- data.frame(row.names=anno_original[anno_original$type =="transcript",]$transcript_id,
    SYMBOL=  anno_original[anno_original$type =="transcript",]$transcript_name)

#prepare methylation data for plotting
anno_meth <- granges(bsseq_all_smooth)
meth_per_cpg <- bsseq::getMeth(bsseq_all_smooth, type = "smooth")
anno_patient <- pheno

#get different subgroup annotation
meth_data_all <- anno_meth
mcols(meth_data_all) <-  meth_per_cpg

#get memory
lapply(list( HM_array, IM_array, LM_array ), rm)
rm(anno_meth)
rm(meth_per_cpg)
rm(bsseq_all_smooth)
gc()

#color annotation
all_patient<- c(as.character(anno_patient[anno_patient$Sample_Type=="tumor",]$Patient),paste0(as.character(anno_patient[anno_patient$Sample_Type=="normal",]$Tissue_Source),"_",as.character(anno_patient[anno_patient$Sample_Type=="normal",]$Developmental_Stage) ))
names(all_patient)<- c(as.character(rownames(anno_patient[anno_patient$Sample_Type=="tumor",])),as.character(rownames(anno_patient[anno_patient$Sample_Type=="normal",])))
all_patient <- all_patient[colnames(mcols(meth_data_all))]
all_patient <- gsub("cordblood_perinatal", "NEO",all_patient)
all_patient <- gsub("bonemarrow_young_adult", "ADU",all_patient)
all_patient <- gsub("liver_fetal", "FL",all_patient)
all_patient <- gsub("spleen_fetal", "FS",all_patient)
all_patient <- gsub("bonemarrow_juvenile", "JUV",all_patient)


all_patient <- factor(all_patient, levels=c("FS", "FL", "NEO","JUV","ADU",  "D117", "D129", "D217", "I217", "D213", "D124", "D123", "D360"))

all_color <- c(  ADU = "#E5E5E5", JUV = "#A7D9CB", NEO = "#978474", FL = "#6D7A9F", FS = "#595959",
                   D117 = "#0058b4", D129 = "#2188c9", 
                   D217 = "#fbbb25", I217 = "#fca349", 
                   D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220")[levels(all_patient)]
Epigenotype_color = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4", wildtype ="#ababab")

#parameters
#region
fontSize <- 12
#select region
#roi <- anno[grep("^FBP2$", anno$gene_name),][1,]
# roi <- GRanges(
#    seqnames = c("chr5", "chr1", "chr12", "chr6", "chr3", "chr1", "chr17", "chr3", "chr11", "chr2", "chr11"),
#    ranges = IRanges(start = c(149780031, 26642500, 9895867, 109682438, 107805978, 208071775, 80271968, 111253926, 119289296, 9599998, 61716486),
#                   end = c(149794830, 26651000, 9936715, 109710665, 107811844, 208119712, 80276186, 111264698, 119296214, 9711622, 61724653)))
# roi$transcript_id <- c("CD74", "CD52", "CD69", "CD164", "CD47", "CD34", "CD7", 
#     "CD96", "CD90THY1", "ADAM17", "BEST1")
    
roi <- GRanges(
   seqnames = c("chr1","chr1","chr1","chr1","chr11", "chr5", "chr11", "chr21", "chr3", "chr12", "chr2", "chr1", "chr1", "chr1"),
   ranges = IRanges(start = c(26630140, 26727773, 26642080, 26727684, 128554112, 76244805, 119290449, 36252213, 111254513, 66212596, 74198395, 164616841, 8916352, 26639916),
                  end = c(26635920,26732593 , 26648772, 26742939,128563997, 76254968, 119295937, 36266709, 111269948, 66281293, 74233481, 164631308, 8941474, 26648286)))
roi$transcript_id <- c( "CD52_DMR1","CD52_DMR2", "CD52_DMR3","CD52_DMR4","FLI1", "CRHBP", "THY1", "RUNX1", "CD96", "HMGA2", "TET3", "PBX1", "ENO1", "CD52")

#run plotting
for(j in c(0, 500, 1500)){
    ext <- j
#plot regions
for (i in 1:length(roi)){

        print(paste0(i, " starts"))
        #Define Region
        lim <- c(start(roi[i,]), end(roi[i,]))
        Chr<- as.character(seqnames(roi[i,]))
        range<- GRanges(
        seqnames = Chr,
        ranges = IRanges(start = lim[1]-ext,
                    end = lim[2]+ext))

        #get ideogramm tracks
        itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr,fontcolor="black")

        #genome axis track
        getrack <- GenomeAxisTrack(fontcolor="black")

        #gene annotation
        #original
        grtrack_original <- GeneRegionTrack(anno_original_gtf, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, 
            geneSymbol=TRUE, name="Genecode",collapseTrack=TRUE, fill="darkgray",color="black", cex.feature = 0.5,
            fontsize=fontSize-2,fontcolor.title="black", fontcolor.group="black")
        symbol(grtrack_original) <- as.character(anno_original_transcript_symbol[transcript(grtrack_original),])

        #get dmr annotation
        dmrs_HM <- AnnotationTrack(dmrs_final$EpigenotypeHM[seqnames(dmrs_final$EpigenotypeHM)==Chr,], 
            name="HM vs LM\nDMRs", chromosome=Chr,genome = "hg19", fill="darkgray", fontsize=2)
        dmrs_HMvsNormal <- AnnotationTrack(dmrs_normal$groupHM[seqnames(dmrs_normal$groupHM)==Chr,], 
            name="HM vs cordblood\nDMRs", chromosome=Chr,genome = "hg19", fill="darkgray",fontsize=2)
        dmrs_LMvsNormal <- AnnotationTrack(dmrs_normal$groupLM[seqnames(dmrs_normal$groupLM)==Chr,], 
            name="LM vs cordblood\nDMRs", chromosome=Chr,genome = "hg19", fill="darkgray",fontsize=2)

        #Data Tracks
        #WGBS Methylation 
        #beta_genes<-  meth_data_all[meth_data_all %over% range,]
        #span <-20/length(beta_genes)
        meth_data_all_track <- DataTrack(range = meth_data_all, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c("a"), chromosome = Chr, name = "WGBS",  #span=span,
            groups=all_patient,col=all_color, lwd=2, alpha=.8,alpha.title=1)     #confint
        #Array Methylation
        HM_array_track <- DataTrack(range = HM_array_gr, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c("histogram"), chromosome = Chr, name = "HM 450k",  #span=span,
            col.histogram=Epigenotype_color["HM"],fill.histogram=Epigenotype_color["HM"], lwd=2, alpha=.8,alpha.title=1)     #confint
        IM_array_track <- DataTrack(range = IM_array_gr, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c("histogram"), chromosome = Chr, name = "IM 450k",  #span=span,
        col.histogram=Epigenotype_color["IM"],fill.histogram=Epigenotype_color["IM"], lwd=2, alpha=.8,alpha.title=1)     #confint
        LM_array_track <- DataTrack(range = LM_array_gr, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c("histogram"), chromosome = Chr, name = "LM 450k",  #span=span,
            col.histogram=Epigenotype_color["LM"],fill.histogram=Epigenotype_color["LM"],  lwd=2, alpha=.8,alpha.title=1)     #confint
        #Plot track
        #pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext, "_Span",span,"_bsseqSmoot_line.pdf")), width = 10, height = 5)
        pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext,"_bsseqSmoot_Averageline.pdf")), width = 10, height = 10)
        plotTracks(list(itrack , getrack, grtrack_original, 
            dmrs_HM,# dmrs_HMvsNormal,dmrs_LMvsNormal,
            meth_data_all_track,  LM_array_track, IM_array_track, HM_array_track),fontcolor.title="black",
            from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
        dev.off()

        #with dots
        meth_data_all_track <- DataTrack(range = meth_data_all, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c("a", "p"), chromosome = Chr, name = "WGBS",  #span=span,
            groups=all_patient,col=all_color, lwd=2, alpha=.8,alpha.title=1)     #confint

        #Plot track
        #pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext, "_Span",span,"_bsseqSmoot_line.pdf")), width = 10, height = 5)
        pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext,"_bsseqSmoot_line_point.pdf")), width = 10, height = 10)
        plotTracks(list(itrack , getrack, grtrack_original, 
            dmrs_HM, #dmrs_HMvsNormal,dmrs_LMvsNormal,
            meth_data_all_track,  LM_array_track, IM_array_track, HM_array_track),fontcolor.title="black",
            from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
        dev.off()

        #with dots only
        meth_data_all_track <- DataTrack(range = meth_data_all, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c( "p"), chromosome = Chr, name = "WGBS",  #span=span,
            groups=all_patient,col=all_color, lwd=2, alpha=.8,alpha.title=1)     #confint

        #Plot track
        #pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext, "_Span",span,"_bsseqSmoot_line.pdf")), width = 10, height = 5)
        pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext,"_bsseqSmoot_point.pdf")), width = 10, height = 10)
        plotTracks(list(itrack , getrack, grtrack_original, 
            dmrs_HM, #dmrs_HMvsNormal,dmrs_LMvsNormal,
            meth_data_all_track, LM_array_track, IM_array_track, HM_array_track),fontcolor.title="black",
            from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
        dev.off()

            #with dots only
        meth_data_all_track <- DataTrack(range = meth_data_all, genome = "hg19",
            fontcolor.title="black", fontcolor.group="black",
            fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
            type = c( "a", "confint"), chromosome = Chr, name = "WGBS",  #span=span,
            groups=all_patient,col=all_color, lwd=2, alpha=.8,alpha.title=1)     #confint

        #Plot track
        #pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext, "_Span",span,"_bsseqSmoot_line.pdf")), width = 10, height = 5)
        pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots/incl_array_novelNormal",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext,"_bsseqSmoot_line_confint.pdf")), width = 10, height = 10)
        plotTracks(list(itrack , getrack, grtrack_original, 
            dmrs_HM,# dmrs_HMvsNormal,dmrs_LMvsNormal,
            meth_data_all_track, LM_array_track, IM_array_track, HM_array_track),fontcolor.title="black",
            from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
        dev.off()
        print(i) 
    }
    print(j)
}

