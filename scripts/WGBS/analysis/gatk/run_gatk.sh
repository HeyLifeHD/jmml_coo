
# GATK segmentation 
#get segmentation in R
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
seq <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
gr <- GRanges(
  seqnames =names(seq),
  ranges = IRanges(start =  1,
                   end =  seq
 )
)
gr<- gr[seqnames(gr) %in% paste0("chr",1:22)]
auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC",
group="auto")
gr <- keepSeqlevels(gr,auto)

newStyle <- mapSeqlevels(seqlevels(gr),"NCBI")
ncbi_gr <- renameSeqlevels(gr, newStyle)

gr_tiled <- unlist(tile(ncbi_gr, width=1000))
dir.create("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input", recursive=TRUE)
export.bed(gr_tiled, "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.bed")

#bed to interval list
cp /icgc/dkfzlsdf/dmg/otp/production/processing/reference_genomes/bwa06_methylCtools_hs37d5_PhiX_Lambda/hs37d5_PhiX_Lambda.fa /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/
cp /icgc/dkfzlsdf/dmg/otp/production/processing/reference_genomes/bwa06_methylCtools_hs37d5_PhiX_Lambda/hs37d5_PhiX_Lambda.fa.fai /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/

conda activate picard
picard CreateSequenceDictionary R=/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.fa \
    O=/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.dict
picard BedToIntervalList \
      I=/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.bed \
      O=/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.interval_list \
      SD=/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.dict



#get docker running
cd /home/heyj/docker/
sshfs heyj@10.128.130.12:/icgc/project /icgc/dkfzlsdf/project/
sudo mkdir -p  /icgc2/dkfzlsdf/project/
sudo chmod 777 -R /icgc2/dkfzlsdf/project/
sudo mkdir -p  /icgc2/dkfzlsdf/analysis/
sudo chmod 777 -R /icgc2/dkfzlsdf/analysis/
sshfs -o allow_other heyj@10.128.130.12:/icgc/project /icgc2/dkfzlsdf/project/
sshfs -o allow_other heyj@10.128.130.12:/icgc/analysis /icgc2/dkfzlsdf/analysis/

## Loop over samples and the folders
mkdir /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output/
chmod 777 -R /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/

docker run -it --mount type=bind,source=/icgc2/,target=/icgc/ broadinstitute/gatk
OUTPUT_DIR=/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output/

for patient in `ls /icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/*/*/paired/merged-alignment/*bam`;
do
echo $patient
name=`basename $patient`
echo $name
# collect read counts for all samples
gatk CollectReadCounts -I ${patient} \
-L /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.interval_list \
--interval-merging-rule OVERLAPPING_ONLY -O $OUTPUT_DIR/${name}.counts.h5
done

## Create a panel for control sample. One can also put all normals together to create a panel!
#gatk CreateReadCountPanelOfNormals \
#-I /icgc/dkfzlsdf/analysis/C010/cancer_microenvironment/TAM/WGBS/analysis/GATK/healthy-rep1_C010_CME_PBAT_TAM_merged.mdup.bam.counts.h5 \
#-I /icgc/dkfzlsdf/analysis/C010/cancer_microenvironment/TAM/WGBS/analysis/GATK/healthy-rep2_C010_CME_PBAT_TAM_merged.mdup.bam.counts.h5 \
#-I /icgc/dkfzlsdf/analysis/C010/cancer_microenvironment/TAM/WGBS/analysis/GATK/healthy-rep3_C010_CME_PBAT_TAM_merged.mdup.bam.counts.h5 \
#--minimum-interval-median-percentile 5.0 --aximum-zeros-in-sample-percentage 20 \
#-O $OUTPUT_DIR/MG_panel.h5 

#annotate intervals for donoiced copy without normals
gatk AnnotateIntervals --intervals /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.interval_list \
--output /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.annotated_intervals \
--reference /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.fa \
--interval-merging-rule OVERLAPPING_ONLY 

## Standardize and denoise case read counts against the PoN with DenoiseReadCounts
mkdir $OUTPUT_DIR/denoised_copy/
mkdir $OUTPUT_DIR/stand_copy/

PID=`ls -l $OUTPUT_DIR | grep  '.h5' | awk '{print $9}'`
echo $PID
for patient in $PID;
do
echo $patient
gatk DenoiseReadCounts -I $OUTPUT_DIR/${patient} \
--annotated-intervals /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.annotated_intervals \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy/${patient}_denoised.tsv 
done
#--count-panel-of-normals $OUTPUT_DIR/MG_panel.h5 \ #should be replaced by --annotated-intervals

## Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios.
mkdir $OUTPUT_DIR/plots/

for patient in $PID;
do
echo $patient
gatk PlotDenoisedCopyRatios \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy/${patient}_denoised.tsv \
--sequence-dictionary /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.dict \
-O $OUTPUT_DIR/plots/ --output-prefix ${patient}
done

#anlyse segmented geneomes and their cnv
#Rscript /home/epicwl/c010-datasets/External/2018-07-20-Coral/3_rep/cnv_calling/segment_denoiseCR.R \
#/home/epicwl/c010-datasets/External/2018-07-20-Coral/3_rep/cnv_calling/odcf_bam/gatk/readCounts/denoised_copy/normal_C010_CME_classII_22_1_merged.mdup.bam.counts.h5_denoised.tsv







#do the same with "normal" pannel
docker run -it --mount type=bind,source=/icgc2/,target=/icgc/ broadinstitute/gatk
OUTPUT_DIR=/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output/

## Create a panel for control sample. One can also put all normals together to create a panel!
gatk CreateReadCountPanelOfNormals \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_D129_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_D213_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_D217_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_I217_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-2_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-2_JMMLC_D217_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00_JMMLC_D123_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00_JMMLC_D124_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor01-1_JMMLC_D213_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor01-2_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor01-2_JMMLC_D124_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor01_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor01_JMMLC_D123_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor01_JMMLC_D124_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor10-1_JMMLC_D129_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor10-1_JMMLC_D217_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor10-2_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor10_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor10_JMMLC_D123_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11-1_JMMLC_D129_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11-1_JMMLC_D217_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11-1_JMMLC_D360_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11-2_JMMLC_D123_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11-2_JMMLC_D124_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11_JMMLC_D123_PBAT_merged.mdup.bam.counts.h5 \
-I /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor11_JMMLC_D124_PBAT_merged.mdup.bam.counts.h5 \
--minimum-interval-median-percentile 5.0 --maximum-zeros-in-sample-percentage 20 \
-O $OUTPUT_DIR/All_panel.h5 

#annotate intervals for donoiced copy without normals
gatk AnnotateIntervals --intervals /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.interval_list \
--output /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.annotated_intervals \
--reference /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.fa \
--interval-merging-rule OVERLAPPING_ONLY 

## Standardize and denoise case read counts against the PoN with DenoiseReadCounts
mkdir $OUTPUT_DIR/denoised_copy_vs_All/
mkdir $OUTPUT_DIR/stand_copy_vs_All/

PID=`ls -l $OUTPUT_DIR | grep  '.h5' | awk '{print $9}'`
echo $PID
for patient in $PID;
do
echo $patient
gatk DenoiseReadCounts -I $OUTPUT_DIR/${patient} \
--count-panel-of-normals $OUTPUT_DIR/All_panel.h5 \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy_vs_All/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy_vs_All/${patient}_denoised.tsv 
done
#--count-panel-of-normals $OUTPUT_DIR/MG_panel.h5 \ #should be replaced by --annotated-intervals
#--annotated-intervals /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.annotated_intervals \

## Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios.
mkdir $OUTPUT_DIR/plots_vs_All/

for patient in $PID;
do
echo $patient
gatk PlotDenoisedCopyRatios \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy_vs_All/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy_vs_All/${patient}_denoised.tsv \
--sequence-dictionary /icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.dict \
-O $OUTPUT_DIR/plots_vs_All/ --output-prefix ${patient}
done

#anlyse segmented geneomes and their cnv
#Rscript /home/epicwl/c010-datasets/External/2018-07-20-Coral/3_rep/cnv_calling/segment_denoiseCR.R \
#/home/epicwl/c010-datasets/External/2018-07-20-Coral/3_rep/cnv_calling/odcf_bam/gatk/readCounts/denoised_copy/normal_C010_CME_classII_22_1_merged.mdup.bam.counts.h5_denoised.tsv
