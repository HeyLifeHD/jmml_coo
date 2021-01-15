#run cellphonedb
conda activate cpdb
#test dataset
cd c010-projects/Mark/JMML/COO/10X/CellPhoneDB/Test_D124
cellphonedb method statistical_analysis cellphonedb_meta_label.txt cellphonedb_count_seurat.txt --threads=10
#run default cellphonedb plotting
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot cellphonedb_meta_label.txt 

#propper dataset
#cd /home/heyj/c010-datasets/Internal/JMMLC/CellphoneDB/D124
#cellphonedb method statistical_analysis cellphonedb_meta_label.txt cellphonedb_count_seurat.txt --threads=5
##run default cellphonedb plotting
#cellphonedb plot dot_plot
#cellphonedb plot heatmap_plot cellphonedb_meta_label.txt

cd /home/heyj/c010-datasets/Internal/JMMLC/CellphoneDB/
for COUNT in `ls */cellphonedb_count_seurat.txt`;do
DIR=`dirname $COUNT`
cellphonedb method statistical_analysis ./${DIR}/cellphonedb_meta_label.txt $COUNT   --output-path ./${DIR}/output/ --threads=5
cellphonedb plot dot_plot --means-path ./${DIR}/output/means.txt --pvalues-path ./${DIR}/output/pvalues.txt \
    --output-path ./${DIR}/output/ --output-name dotplot.pdf
cellphonedb plot heatmap_plot --pvalues-path ./${DIR}/output/pvalues.txt \
    --output-path ./${DIR}/output --count-name heatmap_count.pdf --log-name heatmap_log_count.pdf \
    --count-network-name count_network.txt --interaction-count-name interaction_count.txt ./${DIR}/cellphonedb_meta_label.txt
echo $DIR
done


cd /home/heyj/c010-datasets/Internal/JMMLC/CellphoneDB/
for COUNT in `ls */cellphonedb_count_seurat_2.txt`;do
DIR=`dirname $COUNT`
cellphonedb method statistical_analysis ./${DIR}/cellphonedb_meta_label.txt $COUNT   --output-path ./${DIR}/output/ --threads=5
cellphonedb plot dot_plot --means-path ./${DIR}/output/means.txt --pvalues-path ./${DIR}/output/pvalues.txt \
    --output-path ./${DIR}/output/ --output-name dotplot.pdf
cellphonedb plot heatmap_plot --pvalues-path ./${DIR}/output/pvalues.txt \
    --output-path ./${DIR}/output --count-name heatmap_count.pdf --log-name heatmap_log_count.pdf \
    --count-network-name count_network.txt --interaction-count-name interaction_count.txt ./${DIR}/cellphonedb_meta_label.txt
echo $DIR
done

ls */*count*.txt

cd /home/heyj/c010-datasets/Internal/JMMLC/CellphoneDB/
for COUNT in `ls *M/cellphonedb_count_seurat_2.txt`;do
DIR=`dirname $COUNT`
cellphonedb method statistical_analysis ./${DIR}/cellphonedb_meta_label.txt $COUNT   --output-path ./${DIR}/output/ --threads=5
cellphonedb plot dot_plot --means-path ./${DIR}/output/means.txt --pvalues-path ./${DIR}/output/pvalues.txt \
    --output-path ./${DIR}/output/ --output-name dotplot.pdf
cellphonedb plot heatmap_plot --pvalues-path ./${DIR}/output/pvalues.txt \
    --output-path ./${DIR}/output --count-name heatmap_count.pdf --log-name heatmap_log_count.pdf \
    --count-network-name count_network.txt --interaction-count-name interaction_count.txt ./${DIR}/cellphonedb_meta_label.txt
echo $DIR
done