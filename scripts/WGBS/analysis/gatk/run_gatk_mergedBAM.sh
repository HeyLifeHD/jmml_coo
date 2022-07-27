
#merge bamfiles from the same patients
cd /icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/

#merge sort and index
#conda activate samtools
for pid in `ls /icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/`
do
fname=`basename $pid`
echo $pid
samtools merge -f /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}.bam /icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/$pid/*/paired/merged-alignment/*.bam
#samtools sort  /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}.bam > /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}_sorted.bam
#samtools index /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}_sorted.bam
done


#kick out tumor11_JMMLC_D129_1
samtools merge -f /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/JMMLC_D129_PBAT.bam JMMLC_D129_PBAT/tumor00-1/paired/merged-alignment/tumor00-1_JMMLC_D129_PBAT_merged.mdup.bam JMMLC_D129_PBAT/tumor10-1/paired/merged-alignment/tumor10-1_JMMLC_D129_PBAT_merged.mdup.bam

#sort and index
for pid in `ls /icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/`
do
fname=`basename $pid`
echo $pid
samtools sort  --threads 6 /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}.bam > /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}_sorted.bam
samtools index /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/${pid}_sorted.bam
done

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
mkdir /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output_merged/
chmod 777 -R /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output_merged/

docker run -it --mount type=bind,source=/icgc2/,target=/icgc/ broadinstitute/gatk
OUTPUT_DIR=/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output_merged/

for patient in `ls /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/merged_bam/*_sorted.bam`;
do
echo $patient
name=`basename $patient`
echo $name
# collect read counts for all samples
gatk CollectReadCounts -I ${patient} \
-L /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hg19_1000bp_binning.interval_list \
--interval-merging-rule OVERLAPPING_ONLY -O $OUTPUT_DIR/${name}.counts.h5
done

#######
## Correct this and run
#######
## Create a panel for control sample. One can also put all normals together to create a panel!
gatk CreateReadCountPanelOfNormals \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_D129_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_D213_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_D217_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-1_JMMLC_I217_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-2_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00-2_JMMLC_D217_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00_JMMLC_D117_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00_JMMLC_D123_PBAT_merged.mdup.bam.counts.h5 \
-I /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_output//tumor00_JMMLC_D124_PBAT_merged.mdup.bam.counts.h5 \

--minimum-interval-median-percentile 5.0 --maximum-zeros-in-sample-percentage 20 \
-O $OUTPUT_DIR/All_panel.h5 

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

## Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios.
mkdir $OUTPUT_DIR/plots_vs_All/

for patient in $PID;
do
echo $patient
gatk PlotDenoisedCopyRatios \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy_vs_All/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy_vs_All/${patient}_denoised.tsv \
--sequence-dictionary /omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/gatk/gatk_input/hs37d5_PhiX_Lambda.dict \
-O $OUTPUT_DIR/plots_vs_All/ --output-prefix ${patient}
done

