#Libraries
library(bsseq)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(Gviz)
#load my own data: Epigenotype vs HSC_cb
#Directories
input.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"
#create directory
dir.create("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots")

#load meth
bsseq_all <- readRDS( file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged.rds"))

#load dmrs
dmrs_final<- readRDS(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged", "dmrs_gr_sub_MethDiff.rds"))
dmrs_normal<- readRDS(file.path( "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC", "dmrs_gr_sub_MethDiff.rds"))

#load anno
anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")

#subset annos for symbol parsing
anno_original_transcript_symbol <- data.frame(row.names=anno_original[anno_original$type =="transcript",]$transcript_id,
    SYMBOL=  anno_original[anno_original$type =="transcript",]$transcript_name)
#prepare methylation data for plotting
anno_meth <- granges(bsseq_all)
meth_per_cpg <- bsseq::getMeth(bsseq_all, type = "raw")
anno_patient <- pData(bsseq_all)

#get different subgroup annotation
meth_data_all <- anno_meth
mcols(meth_data_all) <-  meth_per_cpg

#color annotation
all_patient<- c(as.character(anno_patient[anno_patient$Sample_Type=="tumor",]$Patient),as.character(anno_patient[anno_patient$Sample_Type=="normal",]$Tissue_Source))
names(all_patient)<- c(as.character(rownames(anno_patient[anno_patient$Sample_Type=="tumor",])),as.character(rownames(anno_patient[anno_patient$Sample_Type=="normal",])))
all_patient <- all_patient[colnames(mcols(meth_data_all))]
all_patient <- factor(all_patient, levels=c("bonemarrow", "cordblood", "D117", "D129", "D217", "I217", "D213", "D124", "D123", "D360"))
all_color <- c(prenatal = "#252525", cordblood = "#737373", bonemarrow = "#ababab", 
                   D117 = "#0058b4", D129 = "#2188c9", 
                   D217 = "#fbbb25", I217 = "#fca349", 
                   D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220")[levels(all_patient)]


#parameters
#region
ext <- 5000 
fontSize <- 12
#select region
#roi <- anno[grep("^FBP2$", anno$gene_name),][1,]
roi <- GRanges(
   seqnames = c("chr5", "chr1", "chr12", "chr6", "chr3", "chr1", "chr17", "chr3", "chr11", "chr2", "chr11"),
   ranges = IRanges(start = c(149780031, 26642500, 9895867, 109682438, 107805978, 208071775, 80271968, 111253926, 119289296, 9599998, 61716486),
                  end = c(149794830, 26651000, 9936715, 109710665, 107811844, 208119712, 80276186, 111264698, 119296214, 9711622, 61724653)))
roi$transcript_id <- c("CD74", "CD52", "CD69", "CD164", "CD47", "CD34", "CD7", 
    "CD96", "CD90/THY1", "ADAM17", "BEST1")
    
#plot regions
for (i in 6:length(roi)){
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
    #Methylation 
    beta_genes<-  meth_data_all[meth_data_all %over% range,]
    span <-20/length(beta_genes)

    meth_data_all_track <- DataTrack(range = meth_data_all, genome = "hg19",
        fontcolor.title="black", fontcolor.group="black",
        fontsize=fontSize,col.axis="black",ylim=c(-.05,1.05),
        type = c("smooth"), chromosome = Chr, name = "WGBS",  span=span,
        groups=all_patient,col=all_color, lwd=2, alpha=.8)     #confint


    #Plot track
    pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext, "_Span",span,".pdf")), width = 10, height = 5)
    plotTracks(list(itrack , getrack, grtrack_original, 
        dmrs_HM, dmrs_HMvsNormal,dmrs_LMvsNormal,
        meth_data_all_track),fontcolor.title="black",
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()

    print(i)
}





#get different subgroup annotation
# meth_data_LM <-anno_meth
# mcols(meth_data_LM) <-  meth_per_cpg[, rownames(anno_patient[anno_patient$Epigenotype=="LM",])]
# meth_data_IM <-anno_meth
# mcols(meth_data_IM) <-  meth_per_cpg[, rownames(anno_patient[anno_patient$Epigenotype=="IM",])]
# meth_data_HM <-anno_meth
# mcols(meth_data_HM) <-  meth_per_cpg[, rownames(anno_patient[anno_patient$Epigenotype=="HM",])]
# meth_data_normal <-anno_meth
# mcols(meth_data_normal) <-  meth_per_cpg[, rownames(anno_patient[anno_patient$Sample_Type=="normal",])]

#color annotation
# pbat_col = list(
#   Tissue = c(tumor = "#99a637", prenatal = "#252525", cordblood = "#737373", adult_bonemarrow = "#ababab"), 
#   Celltype = c(HSC ="#252525", MPP = "#737373", LMPP = "#9babcf", CD45RACD90 = "#99a637", MEP = "#e62628", CMP = "#f6be13", GMP = "#f57e12"), 
#   Genotype = c(neg = "#529a51", KRAS = "#99a637", PTPN11 = "#007458", wildtype ="#ababab"),
#   Patient = c(D117 = "#0058b4", D129 = "#2188c9", 
#             D217 = "#fbbb25", I217 = "#fca349", 
#             D213 = "#ff6b36", D124 = "#e34e2e", D123 = "#c33126", D360 = "#a41220", 
#             normal ="#ababab"), 
#   Epigenotype = c(HM = "#c33126", IM = "#fbbb25", LM = "#0058b4", wildtype ="#ababab")
#   )
# LM_patient<- anno_patient[anno_patient$Epigenotype=="LM",]$Patient
# LM_patient <- as.character(LM_patient)
# LM_color <- pbat_col[["Patient"]][LM_patient]
# HM_patient<- anno_patient[anno_patient$Epigenotype=="HM",]$Patient
# HM_patient <- as.character(HM_patient)
# HM_color <- pbat_col[["Patient"]][HM_patient]
# IM_patient<- anno_patient[anno_patient$Epigenotype=="IM",]$Patient
# IM_patient <- as.character(IM_patient)
# IM_color <- pbat_col[["Patient"]][IM_patient]
# normal_patient<- anno_patient[anno_patient$Sample_Type=="normal",]$Patient
# normal_patient <- as.character(normal_patient)
# normal_color <- pbat_col[["Patient"]][normal_patient]

 #Data Tracks
    #Methylation 
    # meth_data_LM_track <- DataTrack(range = meth_data_LM, genome = "hg19",
    #     fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
    #     type = c("smooth"), chromosome = Chr, name = "WGBS\nLM", 
    #     groups=LM_patient,col=LM_color)
    # meth_data_IM_track <- DataTrack(range = meth_data_IM, genome = "hg19",
    #     fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
    #     type = c("smooth"), chromosome = Chr, name = "WGBS\nIM", 
    #     groups=IM_patient,col=IM_color)
    # meth_data_HM_track <- DataTrack(range = meth_data_HM, genome = "hg19",
    #     fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
    #     type = c("smooth"), chromosome = Chr, name = "WGBS\nHM", 
    #     groups=HM_patient,col=HM_color)
    # meth_data_normal_track <- DataTrack(range = meth_data_normal, genome = "hg19",
    #     fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
    #     type = c("smooth"), chromosome = Chr, name = "WGBS\nnormal", 
    #     groups=normal_patient,col=normal_color)
  
    #Plot track
    # pdf(file.path("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/locus_plots",paste0(roi[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 10)
    # plotTracks(list(itrack , getrack, grtrack_original, 
    #     meth_data_LM_track, meth_data_IM_track, meth_data_HM_track,meth_data_normal_track),fontcolor.title="black",
    #     from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    # dev.off()