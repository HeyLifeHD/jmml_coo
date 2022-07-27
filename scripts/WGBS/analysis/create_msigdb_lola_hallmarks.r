#create msigdb lola hallmarks
#libraries
library(qusage)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v75)
msigdb <- read.gmt("/home/heyj/omics/groups/OE0219/internal/jmmlc_pbat/data/h.all.v7.1.symbols.gmt")

#get promoter regions for all genes
ENS <- EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
GENES <- genes(ENS)
PROMOTERS <- promoters(GENES, upstream=1500, downstream=1500)

#annotate with symbol
msigdb_gene_regions <- lapply(msigdb, function(x){
    x <- GENES[GENES$symbol %in% x, ]
    x
})

#get promoter regions
msigdb_promoter_regions <- lapply(msigdb_gene_regions, function(x){
    x <- promoters(x, upstream=1500, downstream=1500)
    x
})

#prepare data for lola
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
dir.create( file.path(datasets.dir,"MSigDB_hg19", "Hallmarks_Promoter", "regions"), recursive=TRUE)
for(i in names(msigdb_promoter_regions)){
    export.bed(msigdb_promoter_regions[[i]], file.path(file.path(datasets.dir,"MSigDB_hg19", "Hallmarks_Promoter", "regions", paste0(i, ".bed"))))
}


#same for promoter and gene body
msigdb_gene_regions <- lapply(msigdb, function(x){
    x <- GENES[GENES$symbol %in% x, ]
    x
})


#prepare data for lola
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
dir.create( file.path(datasets.dir,"MSigDB_hg19", "Hallmarks_Gene", "regions"), recursive=TRUE)
for(i in names(msigdb_gene_regions)){
    export.bed(msigdb_gene_regions[[i]], file.path(file.path(datasets.dir,"MSigDB_hg19", "Hallmarks_Gene", "regions", paste0(i, ".bed"))))
}



#make combined object

msigdb_both_regions<- list()
for(i in names(msigdb_promoter_regions)){
   msigdb_both_regions[[i]] <- reduce(unlist(GRangesList(msigdb_promoter_regions[[i]], msigdb_gene_regions[[i]])))
}

#prepare data for lola
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
dir.create( file.path(datasets.dir,"MSigDB_hg19", "Hallmarks_PromoterGeneBody", "regions"), recursive=TRUE)
for(i in names(msigdb_both_regions)){
    export.bed(msigdb_both_regions[[i]], file.path(file.path(datasets.dir,"MSigDB_hg19", "Hallmarks_PromoterGeneBody", "regions", paste0(i, ".bed"))))
}