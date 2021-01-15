#libraries
library(data.table)
library(ggplot2)
library(dplyr)
#Cellphone DB visualization
base.dir <- "/home/heyj/c010-datasets/Internal/JMMLC/CellphoneDB/"
data.dir <- file.path(base.dir, "out")
output.dir <- file.path(base.dir, "custom_plots")
dir.create(output.dir)

#load comparisons of interest
anno <- lapply(comp_oi, function(x){
        pvalues <- fread(file.path(base.dir,x, "output","pvalues.txt"))
        anno <- pvalues[,c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", 
            "gene_b","secreted", "receptor_a", "receptor_b","annotation_strategy", "is_integrin" )]
        anno
})
anno <- unique(do.call("rbind", anno))

comp_oi <- c("HM", "IM", "LM")
pvalues_melt <- lapply(comp_oi, function(x){
        pvalues <- fread(file.path(base.dir,x, "output","pvalues.txt"))
        pvalues_melt <- melt(pvalues[,c(1, 12:ncol(pvalues)), with=FALSE])
        colnames(pvalues_melt) <- c("id_cp_interaction", "interacting_celltype", "pvalue")
        pvalues_melt$comparison <- x
        pvalues_melt$sender  <- sapply(strsplit(as.character(pvalues_melt$interacting_celltype), "|", fixed=TRUE ), "[",1)
        pvalues_melt$receiver  <- sapply(strsplit(as.character(pvalues_melt$interacting_celltype), "|", fixed=TRUE ), "[",2)
        pvalues_melt
})
pvalues_melt_comb <- do.call("rbind", pvalues_melt)
pvalues_melt_comb <- merge(pvalues_melt_comb, anno, by="id_cp_interaction")

means_melt <- lapply(comp_oi, function(x){
        means <- fread(file.path(base.dir,x, "output","means.txt"))
        means_melt <- melt(means[,c(1, 12:ncol(means)), with=FALSE])
        colnames(means_melt) <- c("id_cp_interaction", "interacting_celltype", "mean")
        means_melt$comparison <- x
        means_melt$sender  <- sapply(strsplit(as.character(means_melt$interacting_celltype), "|", fixed=TRUE ), "[",1)
        means_melt$receiver  <- sapply(strsplit(as.character(means_melt$interacting_celltype), "|", fixed=TRUE ), "[",2)
        means_melt
})
means_melt_comb <- do.call("rbind", means_melt)
means_melt_comb <- merge(means_melt_comb, anno, by="id_cp_interaction")

for(i in unique(pvalues_melt_comb$sender)){
    #subset celltype
    pvalues_melt_comb_ct <- pvalues_melt_comb[pvalues_melt_comb$sender ==i, ]
    means_melt_comb_ct <- means_melt_comb[means_melt_comb$sender ==i, ]
    #get id interaction with at least one sign pvalue
    pval_index <- unique(pvalues_melt_comb_ct[pvalues_melt_comb_ct$pvalue < 0.05, ]$id_cp_interaction)
    length(pval_index)
    #get id ingteraction with at least one mean above 0.5
    mean_index <- unique(means_melt_comb_ct[means_melt_comb_ct$mean > 0.05, ]$id_cp_interaction)
    length(mean_index)
    #get overlap
    pval_mean_index <- intersect(pval_index, mean_index)
    length(pval_mean_index)
    #subset
    pvalues_melt_comb_sub <- pvalues_melt_comb_ct[pvalues_melt_comb_ct$id_cp_interaction %in%pval_mean_index,  ]
    means_melt_comb_sub <- means_melt_comb_ct[means_melt_comb_ct$id_cp_interaction %in%pval_mean_index,  ]


    #get ids and celltypes with strongest variance in pvalue
    pvalues_melt_comb_sub$comb_id <- paste0(pvalues_melt_comb_sub$id_cp_interaction, "_", pvalues_melt_comb_sub$interacting_celltype)
    #split by comb id
    pvalues_melt_comb_sub_split<- split(pvalues_melt_comb_sub,pvalues_melt_comb_sub$comb_id)
    length(pvalues_melt_comb_sub_split)

    #get variance
    variance_list <- lapply(pvalues_melt_comb_sub_split, function(x){
        variance_list <-var(x$pvalue)
        variance_list
    })
    variances <- unlist(variance_list)
    variances <- variances[order(variances, decreasing=TRUE)]
    variances_df <- data.frame(id_cp_interaction= sapply(strsplit(names(variances), "_", fixed=TRUE), "[", 1), interacting_celltype= sapply(strsplit(names(variances), "_", fixed=TRUE), "[", 2), variance=as.vector(variances))
    

    #select ids of interest
    id_oi <- unique(variances_df[is.na(variances_df$variance),]$id_cp_interaction)
    id_oi <- unique(variances_df[is.na(variances_df$variance) | variances_df$variance == 0.5,]$id_cp_interaction)
    length(id_oi)
    #get most variable ids of interest
    pvalues_melt_comb_sub_sub <- pvalues_melt_comb_sub[pvalues_melt_comb_sub$id_cp_interaction %in% id_oi ,]
    means_melt_comb_sub_sub <- means_melt_comb_sub[means_melt_comb_sub$id_cp_interaction %in% id_oi,]


    #prepare for plotting
    plot <- cbind(pvalues_melt_comb_sub_sub, means_melt_comb_sub_sub$mean)
    plot$interacting_pair_compaarison <- paste0(plot$interacting_pair, "_", plot$comparison)
    plot$id_cp_interaction <- ordered(plot[order(plot$id_cp_interaction, match(plot$comparison, c("HM","IM", "LM"))),]$id_cp_interaction)
    plot$logpval <- -log10(plot$pvalue)
    plot$log2means <- log2(plot$V2)
    length(unique(plot$interacting_celltype))

    #plot
    dir.create(file.path(output.dir, "LM_HM_IM_comparison", i), recursive=TRUE)
    pdf(file.path(output.dir, "LM_HM_IM_comparison" ,i, "test.pdf"), height=ifelse(length(unique(plot$interacting_pair))/2.5< 16,16,length(unique(plot$interacting_pair))/2.5) , width=ifelse(length(unique(plot$interacting_pair))/2< 10,10,length(unique(plot$interacting_pair))/2) )
            print(ggplot(data = plot, aes(x=interacting_celltype, y=interacting_pair_compaarison)) +
                geom_point(aes(size=logpval, fill=log2means), pch=21) +
                scale_fill_gradient2( midpoint = 0, low="darkblue", high="darkred", name = "Log2 mean (Molecule 1 vs. Molecule 2)")+
                scale_size(name="p-value\n(-log10)") +
                theme_bw() +
                theme(#panel.grid.minor = element_blank(),
                    #panel.grid.major = element_blank(),
                    axis.text=element_text(size=14, colour = "black"),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.text.y = element_text(size=12, colour = "black"),
                    axis.title=element_blank(),
                    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")))
    dev.off()
    print(i)
}

