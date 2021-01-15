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
input.dir <- "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200830_DMR_Model_CelltypeGroup_sub_cbHSC"
dir.create(analysis.dir)

#load data
bsseq_all <- readRDS(file.path(input.dir ,"bsseq", "bsseq_HSC_comb_snpRemoved_repMerged_sub_cbHSC_comparison.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#load HSC DMRs
dmrs_HSC_red<- readRDS(file.path( "/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200801_DMR_hierachy_HSC_comb", "sig_dmrs_5inHalf_sub_anno_reduced.rds"))

#create further subset
dmrs_final_sub <- lapply(dmrs_final, function(x){
    x <- x[x %outside% dmrs_HSC_red,]
    x
})
dmrs_final_sub2 <- lapply(dmrs_final, function(x){
    x <- x[x %over% dmrs_HSC_red,]
    x
})
dmrs_final_sub3 <- list(uniqueLM= dmrs_final$groupLM[dmrs_final$groupLM %outside% dmrs_final$groupHM,],
    uniqueHM=dmrs_final$groupHM[dmrs_final$groupHM %outside% dmrs_final$groupLM,],
    common=dmrs_final$groupLM[dmrs_final$groupLM %over% dmrs_final$groupHM,])

names(dmrs_final_sub)<- paste0(names(dmrs_final_sub), "_hierachy_DMRs_removed")
names(dmrs_final_sub2)<- paste0(names(dmrs_final_sub), "_only_hierachy_DMRs")
lapply(dmrs_final, length)
lapply(dmrs_final_sub, length)
lapply(dmrs_final_sub2, length)
lapply(dmrs_final_sub3, length)

dmrs_final <- c(dmrs_final, dmrs_final_sub, dmrs_final_sub2, dmrs_final_sub3)


#load LOLA
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
regionDB = loadRegionDB(file.path(datasets.dir,"MSigDB_hg19_Hallmarks"))


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
saveRDS(results_regionDB, file.path(analysis.dir, "LOLA_resultTables", "results_MSigDB_hallmarks_hg19.rds"))
results_regionDB<- readRDS(file.path(analysis.dir, "LOLA_resultTables", "results_MSigDB_hallmarks_hg19.rds"))

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
        pdf(file.path(analysis.dir,i, "LOLA", paste0("EnrichLOLA_", j, ".pdf")), width=15, height=15)
        print(g)
        dev.off()
    }
}




#run plotting for all Sign
for(i in names(results_regionDB)){
    result <- results_regionDB[[i]]
    for(j in unique(unique(result$collection))){
        result_sub<- result[result$collection == j,]
        temp <- unique(c(result_sub[result_sub$userSet=="hypo"& result_sub$pValueLog > 2,]$filename, result_sub[result_sub$userSet=="hyper"& result_sub$pValueLog > 2,]$filename))
        if(length(temp)>1){
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
            pdf(file.path(analysis.dir,i, "LOLA", paste0("EnrichLOLA_", j, "_allSig.pdf")), width=15, height=20)
            print(g)
            dev.off()
        }
    }
}








#run enrichment with custom background
library(rtracklayer)
BG <- import.bed("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/EpigenotypeHM/random_BG/random_bg_merged.sorted.bed")
#Run Enrichment
results_regionDB <- list() 
for (i in names(dmrs_final)){
    hypo <- dmrs_final[[i]][which( dmrs_final[[i]]$diff.Meth<0),]
    hyper <- dmrs_final[[i]][which( dmrs_final[[i]]$diff.Methy>0),]
    userSets<- list(hypo, hyper)
    names(userSets)<- c("hypo","hyper")
    dir.create(file.path(analysis.dir, "LOLA"))
    #set  Universe
    userUnisverse <-BG
    #Run analysis
    #results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
    results_regionDB[[i]]= runLOLA(userSets, userUnisverse, regionDB, cores=3)
    print(i)
}

dir.create(file.path(analysis.dir,"LOLA_resultTables"))
saveRDS(results_regionDB, file.path(analysis.dir, "LOLA_resultTables", "results_MSigDB_hg19_customBG.rds"))
results_regionDB<- readRDS(file.path(analysis.dir, "LOLA_resultTables", "results_MSigDB_hg19_customBG.rds"))

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
        pdf(file.path(analysis.dir,"LOLA", paste0("EnrichLOLA_customBG_", j, ".pdf")), width=15, height=15)
        print(g)
        dev.off()
    }
}


