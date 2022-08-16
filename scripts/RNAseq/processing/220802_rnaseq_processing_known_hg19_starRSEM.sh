#singularity
module load anaconda3/2019.07
#create new nextflow environment
conda create -n nextflow_v22
source activate nextflow_v22
conda install -c bioconda nextflow

#update nextflow to developer version
export NXF_EDGE=1
nextflow self-update

#merge RNAseq data
#copy and merge fastq files
cd /omics/odcf/project/OE0219/JMMLC/sequencing/rna_sequencing/view-by-pid/
TARGET_DIR="/omics/groups/OE0219/internal/jmmlc_rnaseq/data"
mkdir -p $TARGET_DIR
for file in `ls | uniq`; do
cat ./${file}/tumor01/paired/*/sequence/*R1.fastq.gz > $TARGET_DIR/${file}.R1.fastq.gz
cat ./${file}/tumor01/paired/*/sequence/*R2.fastq.gz > $TARGET_DIR/${file}.R2.fastq.gz
echo $file
done

#run processing
OUTPUT_DIR="/omics/groups/OE0219/internal/jmmlc_rnaseq/220802_rnaseq_processing_known_hg19_starRSEM"
mkdir $OUTPUT_DIR
cd $OUTPUT_DIR

#create design + config file
wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py
chmod 777 ./fastq_dir_to_samplesheet.py
./fastq_dir_to_samplesheet.py $TARGET_DIR design.csv --strandedness reverse -r1 .R1.fastq.gz -r2 .R2.fastq.gz 

export NXF_HOME=$OUTPUT_DIR
nextflow run nf-core/rnaseq -profile singularity -c rnaseq.config \
    --input design.csv \
    --outdir $OUTPUT_DIR -resume --aligner star_salmon 
    --fasta '/omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa' --gtf '/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf' 
    
    --star_index ''
    
    
    \
    -resume
#try in parallel wit