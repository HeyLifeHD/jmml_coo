conda activate MethylDackel
mkdir /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/Mbias
#Check for mbias
for file in readlink -f `ls /home/heyj/icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/*/*/paired/merged-alignment/*.bam`; do
    fname=`basename $file`
    MethylDackel mbias -@ 10 /home/heyj/icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/GRCh37/seq/hs37d5_PhiX_Lambda.fa $file \
    /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/Mbias/${fname}
    echo $fname
done 

#Extract methylation values with adjusting for mbias
mkdir /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/
for file in  readlink -f `ls /home/heyj/icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/*/*/paired/merged-alignment/*.bam`; do
    fname=`basename $file`
    dirname=`echo $fname | awk -F\. '{print $1}'`
    mkdir /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md//methylationCalls/${dirname}/    
    MethylDackel extract -@ 5 --OT 6,146,2,144 --OB 7,146,12,150  \
    /home/heyj/icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/GRCh37/seq/hs37d5_PhiX_Lambda.fa\
    $file \
    -o /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/${dirname}/${fname}
    echo $dirname
    echo $fname
done 

#replace chromosome names
for file in /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/*/*.bedGraph; do
awk '{print "chr" $0}' $file > $file.chr.bedgraph
echo $file 
done 

