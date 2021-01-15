library(RnBeads.mm10)
library(RnBeads)
library(data.table)
library(ggplot2)
library(plotly)

chromosomes <- paste0("chr", c(1:21, "X", "Y", "M"))
cpg_anno <- rnb.get.annotation(type="CpG", assembly="hg19")

no_of_cpgs <- sum(unlist(lapply(cpg_anno, length)))/2
  
read_file_regions <- function(file, keyword="", lines=3){
    #browser()
    con <- file(file)
    read_in_lines <- readLines(con=con)  
    read_in_data <- (which(read_in_lines==paste0("### ", keyword))+1):(which(read_in_lines==paste0("### ", keyword))+lines)
    data <- as.data.frame(read.table(textConnection(read_in_lines[read_in_data])))
    data <- apply(data, 2, as.character)
    colnames(data) <- as.character(data[1,])
    data <- data[-1,]
    return(data)
}

transform_sum <- function(DT, sum_cols){
  
    changeCols <- colnames(DT)
  DT [,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols]
  DT [,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]
  DT[, lapply(.SD,sum), by=sum_cols]

}

###sections to read in 
### global methylation -- 3 lines
### methylation vs. position 2*125+1
### methylation vs. baseQ 2*42+1
### coverage 201??+1
FOLDER <- "/icgc/dkfzlsdf/project/OE0219/JMMLC_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/"


odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
output.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
bdg_PBAT_files = dir(path = odcf.dir , recursive=TRUE,pattern = "chr.bedgraph$", full.names = TRUE)
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
sample_anno_mark <- read.table(file.path("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data", "sample_anno_mark.txt"), stringsAsFactors=FALSE, header=TRUE)
sample_anno_mark$SampleID <- sample_anno_mark$Name
library(dplyr)
sample_anno<- left_join(sample_anno_PBAT,sample_anno_mark)
rownames(sample_anno)<- sample_anno$SampleID
annotation <- sample_anno
methylation_qc_folder <- paste0(FOLDER, "/", annotation$"ODCF_Patient_ID", "/", tolower(annotation$"Alignment_SampleType"), "/paired/merged-alignment/methylation/merged/methylationCallingMetrics/")

methylation_qc_files <- lapply(methylation_qc_folder, dir, full.names=T)
names(methylation_qc_files) <- annotation$SampleID

### global methylation -- 3 lines
global_methylation <- list()
global_methylation_sum <- list()
methylation_vs_baseQ <- list()
methylation_vs_position <- list()
coverage <- list()

for (sample in annotation$SampleID){

  #global methylation
  global_methylation[[sample]] <- as.data.table(do.call(rbind.data.frame,lapply(methylation_qc_files[[sample]], read_file_regions, keyword="global methylation", lines=3)))
  changeCols <- c("mC", "C")
 global_methylation[[sample]] [,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols]
 global_methylation[[sample]] [,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]
 #summarizing by sample
 global_methylation_sum[[sample]]  <- global_methylation[[sample]] [,.(context, mC, C)]
 global_methylation_sum[[sample]] <- global_methylation_sum[[sample]][, lapply(.SD,sum), by=context]
 global_methylation_sum[[sample]] <- global_methylation_sum[[sample]][, ratio:=mC/(mC+C)]
 global_methylation_sum[[sample]] <- global_methylation_sum[[sample]][,sample_name:=sample]

 
 # methylation and position for mbias
  methylation_vs_position[[sample]] <- as.data.table(do.call(rbind.data.frame,lapply(methylation_qc_files[[sample]], read_file_regions, keyword="methylation vs. position", lines=251)))
  #summarizing 
   methylation_vs_position[[sample]] <- transform_sum( methylation_vs_position[[sample]], sum_cols = c("mate", "pos"))
   methylation_vs_position[[sample]] <- methylation_vs_position[[sample]][,sample_name:=sample]

 
 #methylation vs. base quality
  methylation_vs_baseQ[[sample]] <-  as.data.table(do.call(rbind.data.frame,lapply(methylation_qc_files[[sample]], read_file_regions, keyword="methylation vs. baseQ", lines=85)))
  #summarizing
  methylation_vs_baseQ[[sample]] <- transform_sum( methylation_vs_baseQ[[sample]], sum_cols = c("mate", "baseQ"))
  methylation_vs_baseQ[[sample]] <- methylation_vs_baseQ[[sample]][,sample_name:=sample]
  #coverage
  coverage[[sample]] <-   as.data.table(do.call(rbind.data.frame,lapply(methylation_qc_files[[sample]], read_file_regions, keyword="coverage", lines=202)))
  coverage[[sample]] <- transform_sum( coverage[[sample]], sum_cols = c("cov"))
  coverage[[sample]] <- coverage[[sample]][,sample_name:=sample]
}



methylation_vs_position <- rbindlist(methylation_vs_position)
methylation_vs_baseQ <- rbindlist(methylation_vs_baseQ)
coverage <- rbindlist(coverage)
global_methylation_sum <- rbindlist(global_methylation_sum)
global_methylation_sum <- global_methylation_sum[,coverage:=(mC+C)/no_of_cpgs]
global_methylation_sum <- global_methylation_sum[context=="CH", coverage:=NA]

global_methylation_sum <- global_methylation_sum[,conversion_rate:=1-ratio]
global_methylation_sum <- global_methylation_sum[context=="CG", conversion_rate:=NA]

output.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/QC"
dir.create(output.dir)
save.image(file.path(output.dir, "QC_data.RData"))

#Quality measurements  {.tabset .tabset-fade} 

## Global methylation rate
#load(file.path(output.dir, "QC_data.RData"))
g <- ggplot(data=global_methylation_sum[context=="CG",], aes(ratio, coverage))+geom_point(aes(color=sample_name))+theme_bw()+xlab("Methylation ratio")+ylab("Mean coverage")
ggplotly(g)

pdf(file.path(output.dir, "global_methylation_sum_CG.pdf"))
print(g)
dev.off()

## Conversion rate
g <- ggplot(data=global_methylation_sum[context=="CH",], aes(sample_name, conversion_rate))+geom_point()+theme_bw()+coord_flip()+xlab("Conversion rate")+ylab("")
ggplotly(g)

pdf(file.path(output.dir, "global_methylation_sum_CH.pdf"))
print(g)
dev.off()

##M-bias plot
methylation_vs_position <- methylation_vs_position[,Mratio:=CG.mC/(CG.mC+CG.C)]
methylation_vs_position <- methylation_vs_position[,mate:=as.factor(mate)]
g <- ggplot(data=methylation_vs_position, aes(pos, Mratio))+geom_line(aes(color=sample_name, linetype= mate))+theme_bw()+ylab("Methylation reatio")+xlab("position in read")
ggplotly(g)

pdf(file.path(output.dir, "Mbias.pdf"))
print(g)
dev.off()

methylation_vs_position_sub <- methylation_vs_position[methylation_vs_position$pos %in% 1:50, ]
g <- ggplot(data=methylation_vs_position_sub, aes(pos, Mratio))+geom_line(aes(color=sample_name, linetype= mate))+theme_bw()+ylab("Methylation reatio")+xlab("position in read")
ggplotly(g)
pdf(file.path(output.dir, "Mbias_0_50.pdf"))
print(g)
dev.off()

methylation_vs_position_sub <- methylation_vs_position[methylation_vs_position$pos %in% 90:160, ]
g <- ggplot(data=methylation_vs_position_sub, aes(pos, Mratio))+geom_line(aes(color=sample_name, linetype= mate))+theme_bw()+ylab("Methylation reatio")+xlab("position in read")
ggplotly(g)
pdf(file.path(output.dir, "Mbias_90_160.pdf"))
print(g)
dev.off()

##Coverage
coverage <- coverage[,log_CG:=log(CG)]
g <- ggplot(data=coverage)+geom_bar(aes(cov, log_CG), stat="identity")+facet_wrap(.~sample_name, nrow=6)+theme_bw()+xlab("Coverage, number of reads")+ylab("Number of CpGs, log10")
ggplotly(g)

pdf(file.path(output.dir, "coverage.pdf"))
print(g)
dev.off()

#Quality measurements, merged replicates  {.tabset .tabset-fade} 
## Global methylation rate
global_methylation_sum <- merge(global_methylation_sum, annotation[, c("sample_name", "CellType","SamplePrep")], by="sample_name")

global_methylation_sum_merge <- global_methylation_sum[,lapply(.SD, sum), .SDcols=c("mC", "C"),by=.(CellType, SamplePrep, context)]
global_methylation_sum_merge<- global_methylation_sum_merge[,c("ratio", "coverage", "sample_name") := list(mC/(mC+C), (mC+C)/no_of_cpgs, paste0(global_methylation_sum_merge$CellType, "_", global_methylation_sum_merge$SamplePrep))]

g <- ggplot(data=global_methylation_sum_merge[context=="CG",], aes(ratio, coverage))+geom_point(aes(color=sample_name))+theme_bw()+xlab("Methylation ratio")+ylab("Mean coverage")
pdf(file.path(output.dir, "MethylationRateCpG_merged.pdf"))
print(g)
dev.off()
ggplotly(g)


## Conversion rate
global_methylation_sum_merge <- global_methylation_sum_merge[,conversion_rate:=1-ratio]

g <- ggplot(data=global_methylation_sum_merge[context=="CH",], aes(sample_name, conversion_rate))+geom_point()+theme_bw()+coord_flip()+xlab("Conversion rate")+ylab("")
pdf(file.path(output.dir, "MethylationRateCH_merged.pdf"))
print(g)
dev.off()
ggplotly(g)

##Coverage
annotation$sample_name <-annotation$SampleID
coverage <- merge(coverage, annotation[, c("sample_name", "SampleID","patient", "tumor")], by="sample_name")


coverage_merge <- coverage[,lapply(.SD, sum), .SDcols=c("CG", "CH"),by=.(SampleID, cov)]


coverage_merge<- coverage_merge[,c("log_CG",  "sample_name") := list(log(CG), paste0(coverage_merge$CellType, "_", coverage_merge$SamplePrep))]

g <- ggplot(data=coverage_merge)+geom_bar(aes(cov, log_CG), stat="identity")+facet_wrap(.~sample_name, nrow=6)+theme_bw()+xlab("Coverage, number of reads")+ylab("Number of CpGs, log10")
ggplotly(g)

pdf(file.path(output.dir, "Coverage_merged.pdf"))
print(g)
dev.off()



