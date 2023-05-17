
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
library(DSS)
library(homerkit)
library(ggpubr)
library(parallel)
#Directories
output.dir <- "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/bsseq"
analysis.dir <-  "/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged"
#RNAseq Directories
base.dir<- "/omics/groups/OE0219/internal/jmmlc_rnaseq/220805_rnaseq_knownDEG"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
#bsseq_all <- readRDS(file.path(output.dir , "bsseq_all_snpfil_sub_cov_repMerged.rds"))
dmrs_final<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff.rds"))
dmrs_red<- readRDS(file.path(analysis.dir, "dmrs_gr_sub_MethDiff_anno_reduced.rds"))

#load rna seq data
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
names(DEG_results_list)<- "EpigenotypeHM"

#functions
label_func <- function(x){
  breaks <- x
 # breaks[breaks>=200] <- ">=200"
  breaks
}
bubblePlot <- function(data){
    combined_data <- data
    combined_data$significant<- ifelse(combined_data$q_value_benjamini < (0.05), "Yes", "No" )
    combined_data$cellType<- c(rep("AM", nrow(combined_data)))
    combined_data$percent_of_target_sequences_with_motif <- as.numeric(sapply(strsplit(combined_data$percent_of_target_sequences_with_motif ,"%", fixed=TRUE),`[`, 1))
 # combined_data$log_p.adjusted_neg[is.infinite(combined_data$log_p.adjusted_neg)] <- 200    
    ggplot(data = as.data.frame(combined_data), aes(y=MotifName, x=direction))+coord_fixed()+
    geom_point(aes(size=log_p.adjusted_neg, fill=percent_of_target_sequences_with_motif, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "% of target\nsequences\nwith motif")+
    scale_colour_manual(values=c("grey", "black"), name="q-value < 0.05", drop=FALSE)+
    scale_size(name="p-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$MotifName))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, hjust=1, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab") 
}

#read in homer results
knownRes <- list.files(file.path(analysis.dir, "homer"), pattern="knownResults$", full.names = TRUE, recursive = TRUE, include.dirs= TRUE)
knownRes <- gsub("knownResults", "", knownRes)
knownRes <- knownRes[grep("all",knownRes, invert=TRUE)]
knownRes <- knownRes[grep("customBG",knownRes, invert=TRUE)]

names_knownRes <- strsplit(knownRes, "/", fixed=TRUE)
names_knownRes_sub <- sapply(names_knownRes, function(x) paste0(x[12]))
names(knownRes) <- names_knownRes_sub
homer_list <- lapply(knownRes, function(x)read_homer_output(x))

DAR_list <- list()
DAR_list$EpigenotypeHM <- homer_list

#plot all sig expressed motifs
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
#sub necessary comparisons
DAR_list_sel <- DAR_list[names(DEG_results_list)]

#plot results
for (i in names(DAR_list_sel)){
    DAR_list_sub[[i]]<- lapply(DAR_list_sel[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    result.df<- DEG_results_list[[i]]
    result.df$gene_id  <- toupper(rownames(result.df))
    x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$gene_id)), ]
    
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$q_value_benjamini < 0.05),]
    x
    })
    DAR_list_sub_plot[[i]] <- DAR_list_sub_sig[[i]]
    for (j in names(DAR_list_sub_plot[[i]])){
        DAR_list_sub_plot[[i]][[j]]$direction <- j
        #DAR_list_sub_plot[[i]][[j]] <-  as.data.frame(DAR_list_sub_plot[[i]][[j]]$known_motif_table)
        DAR_list_sub_plot[[i]][[j]] <- DAR_list_sub_plot[[i]][[j]][,c("motif_name", "q_value_benjamini", 
        "percent_of_target_sequences_with_motif","percent_of_background_sequences_with_motif",
        "MotifName","log_p.adjusted_neg" , "direction")]
    }

    DAR_list_sub_plot[[i]]<- do.call("rbind", DAR_list_sub_plot[[i]])
    sig <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$q_value_benjamini < 0.05, ]$MotifName
    DAR_list_sub_plot[[i]] <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$MotifName %in% sig, ]
    
    pdf(file.path(analysis.dir,"homer", paste0("HomerBubble_all_qval0.05_expressed.pdf")), height=40)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i))
    dev.off()
}

#top 10 hypo + hyper of expresed
#plot results
for (i in names(DAR_list_sel)){
    DAR_list_sub[[i]]<- lapply(DAR_list_sel[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    result.df<- DEG_results_list[[i]]
    result.df$gene_id  <- toupper(rownames(result.df))
    x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$gene_id)), ]
    
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$q_value_benjamini < 0.05),]
    x
    })
   
    sig <- list()
    DAR_list_sub_plot[[i]] <- DAR_list_sub_sig[[i]]
    for (j in names(DAR_list_sub_plot[[i]])){
        DAR_list_sub_plot[[i]][[j]]$direction <- j
        #DAR_list_sub_plot[[i]][[j]] <-  as.data.frame(DAR_list_sub_plot[[i]][[j]]$known_motif_table)
        DAR_list_sub_plot[[i]][[j]] <- DAR_list_sub_plot[[i]][[j]][,c("motif_name",  "q_value_benjamini",
        "percent_of_target_sequences_with_motif","percent_of_background_sequences_with_motif",
        "MotifName","log_p.adjusted_neg" , "direction")]
        sig[[j]]<-head(DAR_list_sub_plot[[i]][[j]]$MotifName,10)
    }
    sig <- unlist(sig)
    DAR_list_sub_plot[[i]]<- do.call("rbind", DAR_list_sub_plot[[i]])
    #sig <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$p.adjusted < 0.5, ]$MotifName

    DAR_list_sub_plot[[i]] <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$MotifName %in% sig, ]
    
    pdf(file.path(analysis.dir,"homer", paste0("HomerBubble_top10_qval0.05_expressed.pdf")))
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i))
    dev.off()
}

#plot tf motifs with tf expression being significant
#plot results
for (i in names(DAR_list_sel)){
    DAR_list_sub[[i]]<- lapply(DAR_list_sel[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    result.df<- DEG_results_list[[i]]
    result.df <- result.df[which(result.df$padj < 0.05),]
    result.df$gene_id  <- toupper(rownames(result.df))
    x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$gene_id)), ]
    
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$q_value_benjamini < 0.05),]
    x
    })
    DAR_list_sub_plot[[i]] <- DAR_list_sub_sig[[i]]
    for (j in names(DAR_list_sub_plot[[i]])){
        DAR_list_sub_plot[[i]][[j]]$direction <- j
        #DAR_list_sub_plot[[i]][[j]] <-  as.data.frame(DAR_list_sub_plot[[i]][[j]]$known_motif_table)
        DAR_list_sub_plot[[i]][[j]] <- DAR_list_sub_plot[[i]][[j]][,c("motif_name", "q_value_benjamini", 
        "percent_of_target_sequences_with_motif","percent_of_background_sequences_with_motif",
        "MotifName","log_p.adjusted_neg" , "direction")]
    }

    DAR_list_sub_plot[[i]]<- do.call("rbind", DAR_list_sub_plot[[i]])
    sig <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$q_value_benjamini < 0.05, ]$MotifName
    DAR_list_sub_plot[[i]] <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$MotifName %in% sig, ]
    
    pdf(file.path(analysis.dir,"homer", paste0("HomerBubble_all_qval0.05_TFexpressionPadj0.05.pdf")), height=14)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i))
    dev.off()
}




#plot tf expression of significant TFs
adjp <- 0.05
known_TF_sig<- list()
known_TF_sig_all<- list()
known_TF_sig_all_mgi_symbol<- list()
DAR_list_sel <- DAR_list[names(DEG_results_list)]
#remove $BMDM_healthy_vs_BMDM_tumor due to now sig combination
for(i in names(DAR_list_sel)){
    result.df<- DEG_results_list[[i]]
    result.sig <- result.df[which(result.df$padj < adjp), ]
    result.sig$mgi_symbol  <- toupper(as.character(rownames(result.sig)))

    for(j in names(DAR_list_sel[[i]])){
        known_TF <- DAR_list_sel[[i]][[j]]$known_motif_table
        known_TF <- known_TF[known_TF$q_value_benjamini < 0.05,]
        known_TF$MotifName <-  sapply(strsplit(as.character(known_TF$motif_name) ,"(", fixed=TRUE),`[`, 1)
        known_TF$MotifName <-  toupper(sapply(strsplit(known_TF$MotifName ,"/", fixed=TRUE),`[`, 1))
        known_TF_sig[[j]] <- known_TF[which(known_TF$MotifName %in%  result.sig$mgi_symbol) ,]
    }
    known_TF_sig_all[[i]] <- known_TF_sig
    known_TF_sig_all_mgi_symbol[[i]] <- lapply(known_TF_sig, function(x)x$MotifName)
}
known_TF_sig_all_mgi_symbol


#plot fold change of expression
library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]
for(i in names(DEG_results_list)){
    TF <- c(known_TF_sig_all_mgi_symbol[[i]]$hypo,known_TF_sig_all_mgi_symbol[[i]]$hyper)
    plot <- DEG_results_list[[i]][toupper(rownames(DEG_results_list[[i]]))%in%TF, ]
    plot$SYMBOL <- rownames(plot)
    plot$Change <- ifelse(plot$log2FoldChange>0, "blue","red")

    pdf(file.path(analysis.dir,"homer", paste0("DEG_TF_expr_qval",adjp,"_TF_barplot.pdf")), width=4)
    print(ggbarplot(plot, x="SYMBOL", y="log2FoldChange",fill="Change",palette=col,ylab="log2 fold change" , order=rownames(plot)[order(rownames(plot), decreasing=TRUE)])
    + rremove("ylab")+ rremove("legend")+rotate_x_text(angle = 45)+ coord_flip() )
    dev.off()
}












