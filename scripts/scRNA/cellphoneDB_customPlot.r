#libraries
library(data.table)
library(ggplot2)
library(dplyr)
#Cellphone DB visualization
base.dir <- "/home/heyj/c010-datasets/Internal/JMMLC/CellphoneDB/"
data.dir <- file.path(base.dir, "out")
output.dir <- file.path(base.dir, "custom_plots")
dir.create(output.dir)

#loop over following analysis
pvalue_files <-list.files(base.dir, "pvalues.txt",full.names=TRUE, recursive=TRUE)
set_oI <-sapply(strsplit(pvalue_files, "/", fixed=TRUE), "[",9 )
#forloop 
for(i in set_oI){
#load data
pvalues <- fread(file.path(base.dir,i, "output","pvalues.txt"))
means <- fread(file.path(base.dir,i, "output","means.txt"))
anno <- pvalues[,c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", 
    "gene_b","secreted", "receptor_a", "receptor_b","annotation_strategy", "is_integrin" )]

#loop over all celltypes
all_comp <- colnames(pvalues[,-c(1:11)])
sender_receiver  <- data.frame(comp= all_comp, sender=sapply(strsplit(all_comp, "|", fixed=TRUE ), "[",1), receiver=sapply(strsplit(all_comp, "|", fixed=TRUE ), "[",2))
all_celltypes<- unique(unlist(strsplit(all_comp, "|", fixed=TRUE )))
#select celltype interactions of interest
    for(celltype_oi in all_celltypes){
        #celltype_oi <- "HSC"
        cell_idx <- as.character(sender_receiver[sender_receiver$sender ==celltype_oi,  ]$comp)

        #order them
        cell_idx <- c("id_cp_interaction",cell_idx)
        #subset pvalues and means
        rownames(pvalues)<- pvalues$id_cp_interaction
        pvalues_sub_ct <- pvalues[,..cell_idx]
        means_sub_ct <- means[,..cell_idx]

        #select interactions that are significant in at least on of the above specified interaction
        pvalue_idx <- pvalues_sub_ct$id_cp_interaction[rowSums(pvalues_sub_ct < 0.05) >=1]
        if(length(pvalue_idx)>1){
        #select interactions that have a mean bigger than threshold
        mean_idx <- means_sub_ct$id_cp_interaction[rowSums(means_sub_ct > 0.5) >=1]
        #both passing threshold
        idx <- pvalue_idx[pvalue_idx %in% mean_idx]
        length(cell_idx)
        #subset 
        pvalues_sub <- pvalues_sub_ct[pvalues_sub_ct$id_cp_interaction %in% idx,]
        means_sub <- means_sub_ct[means_sub_ct$id_cp_interaction %in% idx,]

        #prepare for plotting
        pvalues_sub_melt <- melt(pvalues_sub)
        colnames(pvalues_sub_melt)<- c("id_cp_interaction", "Celltypes", "pvalue")
        means_sub_melt <- melt(means_sub)
        colnames(means_sub_melt)<- c("id_cp_interaction", "Celltypes", "mean")
        table(means_sub_melt$id_cp_interaction == pvalues_sub_melt$id_cp_interaction)
        table(means_sub_melt$variable == pvalues_sub_melt$variable)
        plot <- cbind(pvalues_sub_melt,means_sub_melt[,"mean"])
        #add interaction data
        plot <- merge(plot, anno, by="id_cp_interaction")
        #change values for plotting
        plot[plot$pvalue==0,]$pvalue = 0.0009
        plot$logpval <- -log10(plot$pvalue)
        plot$log2means <- log2(plot$mean)
        #order by clustering of mean
        temp <- means_sub
        temp <- temp[,!"id_cp_interaction"]
        temp <- log2(temp)
        minim <- min(unlist(lapply(temp, function(x)min(x[!is.infinite(x)]))))
        invisible(lapply(names(temp),function(.name) set(temp, which(is.infinite(temp[[.name]])), j = .name,value =minim)))
        temp <- as.matrix(t(temp))
        if(ncol(temp)>2 & sum(matrixStats::rowVars(t(temp))==0)==0){
        dissimilarity <- 1 - cor(temp, use="pairwise.complete.obs")
        distance <- as.dist(dissimilarity)
        hrld<- hclust(distance)
        dend<- hrld%>% as.dendrogram 
        plot$id_cp_interaction <- as.factor(plot$id_cp_interaction)
        plot$id_cp_interaction <- ordered(plot$id_cp_interaction,means_sub$"id_cp_interaction"[order.dendrogram(dend)])
        }
        dir.create(file.path(output.dir,i, celltype_oi), recursive=TRUE)
        pdf(file.path(output.dir, i,celltype_oi, "Interaction_allSign_meanClusteredRow.pdf"), height=ifelse(length(unique(plot$interacting_pair))/5< 7,7,length(unique(plot$interacting_pair))/5) , width=ifelse(length(unique(plot$interacting_pair))/5< 12,12,length(unique(plot$interacting_pair))/2))
        print(ggplot(data = plot, aes(x=Celltypes, y=interacting_pair)) +
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

        #order by clustering of pvalue
        temp <- pvalues_sub
        temp <- temp[,!"id_cp_interaction"]
        temp <- -log10(temp)
        maxim <- max(unlist(lapply(temp, function(x)max(x[!is.infinite(x)]))))
        invisible(lapply(names(temp),function(.name) set(temp, which(is.infinite(temp[[.name]])), j = .name,value =maxim)))
        temp <- as.matrix(t(temp))
          if(ncol(temp)>2 & sum(matrixStats::rowVars(t(temp))==0)==0){
        dissimilarity <- 1 - cor(temp, use="pairwise.complete.obs")
        distance <- as.dist(dissimilarity)
        hrld<- hclust(distance)
        dend<- hrld%>% as.dendrogram 
        plot$id_cp_interaction <- as.factor(plot$id_cp_interaction)
        plot$id_cp_interaction <- ordered(plot$id_cp_interaction,means_sub$"id_cp_interaction"[order.dendrogram(dend)])
          }
        pdf(file.path(output.dir, i,celltype_oi, "Interaction_allSign_pvalClusteredRow.pdf"), height=ifelse(length(unique(plot$interacting_pair))/5< 7,7,length(unique(plot$interacting_pair))/5) , width=ifelse(length(unique(plot$interacting_pair))/5< 12,12,length(unique(plot$interacting_pair))/2))
        print(ggplot(data = plot, aes(x=Celltypes, y=interacting_pair)) +
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
        print(celltype_oi)
        }}
                print(i)

}





#####To be adapted
#same also for columns
temp <- means_sub
temp <- temp[,!"id_cp_interaction"]
temp <- log2(temp)
minim <- min(unlist(lapply(temp, function(x)min(x[!is.infinite(x)]))))
invisible(lapply(names(temp),function(.name) set(temp, which(is.infinite(temp[[.name]])), j = .name,value =minim)))
temp <- as.matrix(temp)
dissimilarity <- 1 - cor(temp, use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
plot$Celltypes <- as.factor(plot$Celltypes)
plot$Celltypes <- ordered(plot$Celltypes,colnames(means_sub[,!"id_cp_interaction"])[order.dendrogram(dend)])

pdf(file.path(output.dir, i,"HSC_Interaction_allSign_meanClusteredRowCol.pdf"), height=25, width=10)
ggplot(data = plot, aes(x=Celltypes, y=interacting_pair)) +
    geom_point(aes(size=logpval, fill=log2means), pch=21) +
    scale_fill_gradient2( midpoint = 0, low="darkblue", high="darkred", name = "Log2 mean (Molecule 1 vs. Molecule 2)")+
    scale_size(name="p-value\n(-log10)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()

#same also for columns
temp <- pvalues_sub
temp <- temp[,!"id_cp_interaction"]
temp <- -log10(temp)
maxim <- max(unlist(lapply(temp, function(x)max(x[!is.infinite(x)]))))
invisible(lapply(names(temp),function(.name) set(temp, which(is.infinite(temp[[.name]])), j = .name,value =maxim)))
temp <- as.matrix(temp)
dissimilarity <- 1 - cor(temp, use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
plot$Celltypes <- as.factor(plot$Celltypes)
plot$Celltypes <- ordered(plot$Celltypes,colnames(pvalues_sub[,!"id_cp_interaction"])[order.dendrogram(dend)])
pdf(file.path(output.dir, i,"HSC_Interaction_allSign_pvalClusteredRowCol.pdf"), height=25, width=10)
ggplot(data = plot, aes(x=Celltypes, y=interacting_pair)) +
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
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()

pdf(file.path(output.dir, i,"HSC_Interaction_allSign_pvalClusteredRowCol_transposed.pdf"), height=10, width=33)
ggplot(data = plot, aes(x=interacting_pair, y=Celltypes)) +
    geom_point(aes(size=logpval, fill=log2means), pch=21) +
    scale_fill_gradient2( midpoint = 0, low="darkblue", high="darkred", name = "Log2 mean (Molecule 1 vs. Molecule 2)")+
    scale_size(name="p-value\n(-log10)") +
    theme_bw() +
    theme(#panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()

#select interactions of interest
gen_oI <- "TNF"
#select interactions that are significant in at least on of the above specified interaction
pvalue_idx <- pvalues_sub_ct$id_cp_interaction[rowSums(pvalues_sub_ct < 0.05) >=1]
#select interactions that have a mean bigger than threshold
mean_idx <- means_sub_ct$id_cp_interaction[rowSums(means_sub_ct > 0.5) >=1]
#both passing threshold
idx <- pvalue_idx[pvalue_idx %in% mean_idx]
#select interactions that are significant in at least on of the above specified interaction
select_idx <- pvalues$id_cp_interaction[grep(gen_oI, pvalues$interacting_pair)]
#select interactions that have a mean bigger than threshold
#both passing threshold
idx <- idx[idx %in% select_idx]
length(idx)
#subset 
pvalues_sub <- pvalues_sub_ct[pvalues_sub_ct$id_cp_interaction %in% idx,]
means_sub <- means_sub_ct[means_sub_ct$id_cp_interaction %in% idx,]


#prepare for plotting
pvalues_sub_melt <- melt(pvalues_sub)
colnames(pvalues_sub_melt)<- c("id_cp_interaction", "Celltypes", "pvalue")
means_sub_melt <- melt(means_sub)
colnames(means_sub_melt)<- c("id_cp_interaction", "Celltypes", "mean")
table(means_sub_melt$id_cp_interaction == pvalues_sub_melt$id_cp_interaction)
table(means_sub_melt$variable == pvalues_sub_melt$variable)
plot <- cbind(pvalues_sub_melt,means_sub_melt[,"mean"])
#add interaction data
plot <- merge(plot, anno, by="id_cp_interaction")
#change values for plotting
plot[plot$pvalue==0,]$pvalue = 0.0009
plot$logpval <- -log10(plot$pvalue)
plot$log2means <- log2(plot$mean)
#ordercelltype
cell_idx <- colnames(pvalues)[grep(celltype_oi,unique(plot$Celltypes))]
cell_idx <- c(cell_idx[grep(paste0(celltype_oi, "$"),cell_idx)], cell_idx[-grep(paste0(celltype_oi, "$"),cell_idx)])
plot$Celltypes <- ordered(plot$Celltypes,colnames(pvalues_sub[,!"id_cp_interaction"])[order.dendrogram(dend)])

#order by clustering of pvalue
pdf(file.path(output.dir, paste0("Mac_Interaction_",gen_oI,"_pval.pdf")), height=7, width=10)
ggplot(data = plot, aes(x=Celltypes, y=interacting_pair)) +
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
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()