#Directories
#meth
output.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"
#lola
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"

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
library(LOLA)
library(dplyr)
library(forcats)
library(ggpubr)

#dmrs_red<- readRDS(file.path(analysis.dir, "sig_dmrs_3in2_sub_anno_reduced.rds"))
dmrs_red_anal <- readRDS(file.path(analysis.dir, "medecom", "dmrs_red_anal.rds"))
#load medecom results 
medecom.result <- readRDS(file.path(analysis.dir, "medecom", "allDMRred_Medecom_results.rds"))
dmr_LMC_list_unlist <- readRDS(file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_list_unlist.rds"))
dmr_LMC_list_promoters_unlist <- readRDS(file.path(analysis.dir, "medecom","downstream_K8_lambda0.001", "dmr_LMC_list_promoters_unlist.rds"))

#load LOLA
#regionDB = loadRegionDB(file.path(datasets.dir,"scratch/ns5bc/resources/regions/LOLACore/hg19/"))
#saveRDS(regionDB, file.path(datasets.dir, "regionDB_hg19.rds"))
regionDB <- readRDS( file.path(datasets.dir, "regionDB_hg19.rds"))

#Create Lists

results_regionDB <- list()
#create lists for analysis 
dmr_LMC_03<- dmr_LMC_list_unlist[grep("LMC_03", names(dmr_LMC_list_unlist))]
dmr_LMC_04<- dmr_LMC_list_unlist[grep("LMC_04", names(dmr_LMC_list_unlist))]
dmr_LMC_05<- dmr_LMC_list_unlist[grep("LMC_05", names(dmr_LMC_list_unlist))]
dmr_LMC_07<- dmr_LMC_list_unlist[grep("LMC_07", names(dmr_LMC_list_unlist))]
dmr_LMC_08<- dmr_LMC_list_unlist[grep("LMC_08", names(dmr_LMC_list_unlist))]

test <- list(dmr_LMC_03, dmr_LMC_04, dmr_LMC_05, dmr_LMC_07, dmr_LMC_08)
names(test)<- c("dmr_LMC_03_combined", "dmr_LMC_04_combined", "dmr_LMC_05_combined","dmr_LMC_07_combined","dmr_LMC_08_combined")
#Run Enrichment
results_Core <- list()
results_genomicRegions <- list()
results_homer <- list()
results_chipSeq <- list()
results_msigdb <- list() 
for (i in names(test)){
dir.create(file.path(analysis.dir,"medecom", "downstream_K8_lambda0.001", "LOLA",i), recursive=TRUE)
#stratify open and closed
userSets<- test[[i]]
#set  Universe
userUnisverse <-dmrs_red_anal
#Run analysis
#results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
results_regionDB[[i]]= runLOLA(userSets, userUnisverse, regionDB, cores=3)
print(i)
}
dir.create(file.path(analysis.dir,"medecom", "downstream_K8_lambda0.001", "LOLA_resultTables"))
saveRDS(results_regionDB, file.path(analysis.dir, "medecom", "downstream_K8_lambda0.001", "LOLA_resultTables", "results_regionDB_combined.rds"))
results_regionDB<- readRDS(file.path(analysis.dir,"medecom", "downstream_K8_lambda0.001",  "LOLA_resultTables", "results_regionDB_combined.rds"))

#Plot genomic regions in bubble plot
#function 
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#plotting
library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]



for(i in names(results_regionDB)){
    result <- results_regionDB[[i]]
    dir.create(file.path(analysis.dir,"medecom", "downstream_K8_lambda0.001", "LOLA",i))
    for(j in unique(unique(result$collection))){
        result_sub<- result[result$collection == j,]
            temp <- list()
        for( us in unique(result_sub$userSet)){
            temp[[us]]<- unique(head(result_sub[result_sub$userSet==us,]$filename, 10))
        }
        #temp <- unique(c(head(result_sub[result_sub$userSet=="hypo",]$filename, 20), head(result_sub[result_sub$userSet=="hyper",]$filename, 20)))
        temp <- unlist(temp)

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
        pdf(file.path(analysis.dir,"medecom", "downstream_K8_lambda0.001", "LOLA",i, paste0("EnrichLOLA_", j, ".pdf")), width=15, height=20)
        print(g)
        dev.off()
    }
}