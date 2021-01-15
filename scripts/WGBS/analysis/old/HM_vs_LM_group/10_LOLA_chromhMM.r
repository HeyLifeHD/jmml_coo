##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
#run lola on  DMRs chromHMM

#Libraries
library(DSS)
library(bsseq)
library(doParallel)
library(ChIPseeker)
library(foreach)
library(SummarizedExperiment)
library(doMC)
library(rtracklayer)
library(HDF5Array)
library(org.Mm.eg.db)
library(ggpubr)
library(randomcoloR)
library(RColorBrewer)
library(VennDiagram)
library(ChIPpeakAnno)
library(dendextend)
library(LOLA)
library(data.table)
#Directories
odcf.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

dmrs_final<- list(HM_vs_LM=dmrs)

#prepare data for lola
bed_files <- paste0("/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/roadmap_tracks_jmml/",list.files("/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/roadmap_tracks_jmml"))
bed_files<-lapply(bed_files,function(x){
    x <- fread(x)
    x <- x[,1:4]
    x <- as.data.frame(x)
    colnames(x)<- c("chromosome", "start","end", "state")
    x <-makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
    x
})
names(bed_files) <- c("ESC.H1", "STRM.MRW.MSC", "BLD.CD14.PC", "BLD.CD15.PC", "BLD.CD19.PPC","BLD.CD3.PPC","BLD.CD34.PC", "BLD.CD56.PC","BLD.PER.MONUC.PC", "SPLN", "BLD.GM12878", "BLD.K562.CNCR")
#split them by state
bed_files <-lapply(bed_files, function(x){
    x <- split(x, x$state)
    x
})

#export them
dir.create( file.path(datasets.dir,"jmmlc"), recursive=TRUE)
for(i in names(bed_files)){
    dir.create( file.path(datasets.dir,"jmmlc", i,"regions"), recursive=TRUE)
    names(bed_files[[i]])<- gsub("/","-", names(bed_files[[i]]))
    for(j in names(bed_files[[i]]))
        export.bed(bed_files[[i]][[j]], file.path(file.path(datasets.dir,"jmmlc",i, "regions", paste0(j, ".bed"))))
}
#load lola
regionDB <- loadRegionDB(file.path(datasets.dir,"jmmlc"))


#Run Enrichment
results_regionDB <- list() 
for (i in names(dmrs_final)){
    hypo <- dmrs_final[[i]][which( dmrs_final[[i]]$diff.Meth<0),]
    hyper <- dmrs_final[[i]][which( dmrs_final[[i]]$diff.Methy>0),]
    userSets<- list(hypo, hyper)
    names(userSets)<- c("hypo","hyper")
    dir.create(file.path(analysis.dir, "LOLA"))
    #set  Universe
    userUnisverse <-dmrs_final[[i]]
    #Run analysis
    #results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
    results_regionDB[[i]]= runLOLA(userSets, userUnisverse, regionDB, cores=3)
    print(i)
}


saveRDS(results_regionDB, file.path(analysis.dir, "LOLA_resultTables", "results_chromHMM_jmmlc.rds"))
results_regionDB<- readRDS(file.path(analysis.dir, "LOLA_resultTables", "results_chromHMM_jmmlc.rds"))

#function
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#run plotting
for(i in names(results_regionDB)){
    result <- results_regionDB[[i]]
    for(j in unique(unique(result$collection))){
        result_sub<- result[result$collection == j,]
        temp <- unique(c(head(result_sub[result_sub$userSet=="hypo",]$filename, 20), head(result_sub[result_sub$userSet=="hyper",]$filename, 20)))
        result_sub2 <- result_sub[result_sub$filename %in% temp, ]


        #data preparation
        combined_data <- result_sub2[,c("userSet","dbSet", "pValueLog", "oddsRatio" ,"filename", "qValue")]#[userSet=="closed",]
        #combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
        combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
        #change infinite values
        #combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
        combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
        combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)
        combined_data$filename <-as.factor(combined_data$filename)
        combined_data$filename <- ordered(combined_data$filename,c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF-Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC",  "14_ReprPCWk", "15_Quies"))
        #plot
        g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
                geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
                scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
                scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
                scale_size(name="P-value\n(-log10)", labels = label_func) +
                scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
                theme(text =element_text(size=14, color="black", family = "sans"),
                    axis.ticks = element_blank(), axis.line = element_blank(), 
                    axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
                    axis.text.y=element_text(size=12, family="sans", color="black"))+
                scale_x_discrete(name=NULL)+
                theme(legend.text=element_text(size=12, family="sans"), 
                    legend.title=element_text(size=12, family= "sans"),
                    legend.background = element_rect(fill="white", color="white"),
                    panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
                    legend.key = element_rect(fill="white"))+rremove("ylab")
        pdf(file.path(analysis.dir,"LOLA", paste0("EnrichLOLA_", j, ".pdf")), width=9, height=7)
        print(g)
        dev.off()
    }
} 

