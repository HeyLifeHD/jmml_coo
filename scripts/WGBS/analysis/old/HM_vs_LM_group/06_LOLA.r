##Joschka hey
#20.05.2020
#PBAT analysis JMMLC
##Run Lola enrichment

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
#Directories
odcf.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/methylationCalls/"
input.dir <- "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

#load data
bsseq_all <- readRDS(file.path(input.dir , "bsseq","bsseq_all_snpfil_sub.rds"))
dmrs <- readRDS(file.path(analysis.dir, "sig_dmrs_5in4_sub_anno.rds"))

dmrs_final<- list(HM_vs_LM=dmrs)
#export lists

#load LOLA
datasets.dir <- "/C010-datasets/Internal/COPD/enrichment_databases/"
regionDB = loadRegionDB(file.path(datasets.dir,"scratch/ns5bc/resources/regions/LOLACore/hg19/"))


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

dir.create(file.path(analysis.dir,"LOLA_resultTables"))
saveRDS(results_regionDB, file.path(analysis.dir, "LOLA_resultTables", "results_regionDB.rds"))
results_regionDB<- readRDS(file.path(analysis.dir, "LOLA_resultTables", "results_regionDB.rds"))

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
        pdf(file.path(analysis.dir,"LOLA", paste0("EnrichLOLA_", j, ".pdf")), width=15, height=20)
        print(g)
        dev.off()
    }
}

