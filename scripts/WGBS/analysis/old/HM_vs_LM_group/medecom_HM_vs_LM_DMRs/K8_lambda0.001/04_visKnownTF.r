#Directories
input.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/"
analysis.dir <-  "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200519_DMR_HM_vs_LM"

#Libraries
library(homerkit)
library(ggpubr)
library(parallel)

#load data

#read in homer results
knownRes <- list.files(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer"), pattern="knownResults$", full.names = TRUE, recursive = TRUE, include.dirs= TRUE)
knownRes <- gsub("knownResults", "", knownRes)
knownRes <- knownRes[grep("bg", knownRes, invert=TRUE)]

names_knownRes <- strsplit(knownRes, "/", fixed=TRUE)
names_knownRes_sub <- sapply(names_knownRes, function(x) paste0(x[16], "_", x[17]))
names(knownRes) <- names_knownRes_sub
homer_list <- mclapply(knownRes, function(x)read_homer_output(x))
#subdevide
DAR_list <- list()
for(i in unique(sapply(names_knownRes, function(x) x[16]))){
    DAR_list[[i]]<- homer_list[grep(paste0("^",i), names(homer_list))]
    temp <- strsplit(names(DAR_list[[i]]), "_", fixed=TRUE)
    temp <- sapply(temp, function(x)x[length(x)])
    names( DAR_list[[i]]) <- temp
}

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
          axis.text.x=element_text(size=12, angle = 45, hjust=1, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab") 
}

#plot all sig motifs
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
for (i in names(DAR_list)){
    DAR_list_sub[[i]]<- lapply(DAR_list[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    #result.df<- DEG_results_list[[i]]
    #result.df$mgi_symbol  <- toupper(rownames(result.df))
    #x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$mgi_symbol)), ]
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$p.adjusted < 0.05),]
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
    
    pdf(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer",i, paste0("HomerBubble_all_qval0.05.pdf")), height=40)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i))
    dev.off()
}


#top 10 hypo + hyper 
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
for (i in names(DAR_list)){
    DAR_list_sub[[i]]<- lapply(DAR_list[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    #result.df<- DEG_results_list[[i]]
    #result.df$mgi_symbol  <- toupper(rownames(result.df))
    #x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$mgi_symbol)), ]
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$p.adjusted < 0.05),]
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
    
    pdf(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer", i,paste0("HomerBubble_top10_qval0.05.pdf")))
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i))
    dev.off()

}



#do the same for a combined list
unique(sapply(names_knownRes, function(x) x[15]))
DAR_list <- list()
for(i in c("LMC_05", "LMC_07", "LMC_08") ){
    DAR_list[[i]]<- homer_list[grep(i, names(homer_list))]
    temp <- strsplit(names(DAR_list[[i]]), "LMC", fixed=TRUE)
    temp <- sapply(temp, function(x)paste0(x[3]))
    names( DAR_list[[i]]) <- paste0("LMC:",temp)
}

#plot all sig motifs
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
for (i in names(DAR_list)){
    DAR_list_sub[[i]]<- lapply(DAR_list[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    #result.df<- DEG_results_list[[i]]
    #result.df$mgi_symbol  <- toupper(rownames(result.df))
    #x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$mgi_symbol)), ]
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$p.adjusted < 0.05),]
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
    
    pdf(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer", paste0("combined_",i,"_HomerBubble_all_qval0.05.pdf")), height=60, width=9)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i))
    dev.off()
}

 
#top 10 hypo + hyper
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
for (i in names(DAR_list)){
    DAR_list_sub[[i]]<- lapply(DAR_list[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    #result.df<- DEG_results_list[[i]]
    #result.df$mgi_symbol  <- toupper(rownames(result.df))
    #x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$mgi_symbol)), ]
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$p.adjusted < 0.05),]
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
    
    pdf(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer", paste0("combined_",i,"_HomerBubble_top10_qval0.05.pdf")), width=12)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i)+
    theme(text = element_text(size=5),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
    ))
    dev.off()
}



#top 5 hypo + hyper 
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
for (i in names(DAR_list)){
    DAR_list_sub[[i]]<- lapply(DAR_list[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    #result.df<- DEG_results_list[[i]]
    #result.df$mgi_symbol  <- toupper(rownames(result.df))
    #x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$mgi_symbol)), ]
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$p.adjusted < 0.05),]
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
        sig[[j]]<-head(DAR_list_sub_plot[[i]][[j]]$MotifName,5)
    }
    sig <- unlist(sig)
    DAR_list_sub_plot[[i]]<- do.call("rbind", DAR_list_sub_plot[[i]])
    #sig <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$p.adjusted < 0.5, ]$MotifName

    DAR_list_sub_plot[[i]] <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$MotifName %in% sig, ]
    
    pdf(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer", paste0("combined_",i,"_HomerBubble_top5_qval0.05.pdf")), width=9)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i)+
    theme(text = element_text(size=8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
    ))
    dev.off()

}


#top 5 hypo + hyper 
DAR_list_sub <- list()
DAR_list_sub_sig <- list()
DAR_list_sub_plot<-list()
for (i in names(DAR_list)){
    DAR_list_sub[[i]]<- lapply(DAR_list[[i]], function(x){
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$motif_name ,"(", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- sapply(strsplit(x$known_motif_table$MotifName ,"/", fixed=TRUE),`[`, 1)
    x$known_motif_table$MotifName <- toupper(x$known_motif_table$MotifName)
    #result.df<- DEG_results_list[[i]]
    #result.df$mgi_symbol  <- toupper(rownames(result.df))
    #x$known_motif_table <- x$known_motif_table[which(x$known_motif_table$MotifName %in% unique(result.df$mgi_symbol)), ]
    #x$known_motif_table$p.adjusted <- p.adjust(exp(x$known_motif_table$log_p_value), "BH")
    #x$known_motif_table$log_p.adjusted_neg <- -(log10(x$known_motif_table$p.adjusted))
    x$known_motif_table$log_p.adjusted_neg <- -(x$known_motif_table$log_p_value)
    x
    })
    DAR_list_sub_sig[[i]]<- lapply(DAR_list_sub[[i]], function(x){
    #x <- x$known_motif_table[which(x$known_motif_table$p.adjusted < 0.1),]
    x <- as.data.frame(x$known_motif_table)
    #x <- x[which(x$p.adjusted < 0.05),]
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
        sig[[j]]<-head(DAR_list_sub_plot[[i]][[j]]$MotifName,3)
    }
    sig <- unlist(sig)
    DAR_list_sub_plot[[i]]<- do.call("rbind", DAR_list_sub_plot[[i]])
    #sig <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$p.adjusted < 0.5, ]$MotifName

    DAR_list_sub_plot[[i]] <- DAR_list_sub_plot[[i]][DAR_list_sub_plot[[i]]$MotifName %in% sig, ]
    
    pdf(file.path(analysis.dir, "medecom/downstream_combined_K3_lambda0.01/homer", paste0("combined_",i,"_HomerBubble_top3_qval0.05.pdf")), width=9)
    print(bubblePlot(DAR_list_sub_plot[[i]])+ggtitle(i)+
    theme(text = element_text(size=8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
    ))
    dev.off()

}
