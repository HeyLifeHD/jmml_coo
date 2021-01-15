#create msigdb lola hallmarks
#libraries
library(qusage)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v75)
library(LOLA)
#Directories
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

msigdb <- read.gmt("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/h.all.v7.1.symbols.gmt")

msigdb <- list.files("/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/msigdb_collections", pattern="gmt", recursive=TRUE, full.names=TRUE)
temp <- sapply(strsplit(msigdb, "/"), "[", 11)
temp <- sapply(strsplit(temp, ".symbols"), "[", 1)
names(msigdb)<- temp
#get promoter regions for all genes
ENS <- EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
GENES <- genes(ENS)
PROMOTERS <- promoters(GENES, upstream=1500, downstream=1500)

#prepare data for lola
msigdb_data <- list()
msigdb_gene_regions <- list()
msigdb_promoter_regions <- list()
msigdb_both_regions<- list()
for(i in names(msigdb)){
    #read in data
    msigdb_data[[i]] <- read.gmt(msigdb[[i]])
    #annotate with symbol
    msigdb_gene_regions[[i]] <- lapply(msigdb_data[[i]], function(x){
        x <- GENES[GENES$symbol %in% x, ]
        x
    })
    #get promoter regions
    msigdb_promoter_regions[[i]] <- lapply(msigdb_gene_regions[[i]] , function(x){
        x <- promoters(x, upstream=1500, downstream=1500)
        x
})
    print(i)
    #make combined object
    for(j in names(msigdb_promoter_regions[[i]])){
         msigdb_both_regions[[i]][[j]] <- list()
        msigdb_both_regions[[i]][[j]] <- reduce(unlist(GRangesList(msigdb_promoter_regions[[i]][[j]], msigdb_gene_regions[[i]][[j]])))
    }
    dir.create( file.path(datasets.dir,"MSigDB_hg19", i, "regions"), recursive=TRUE)
    for(j in names(msigdb_both_regions[[i]])){
        export.bed(msigdb_both_regions[[i]][[j]], file.path(file.path(datasets.dir,"MSigDB_hg19", i, "regions", paste0(j, ".bed"))))
        print(j)
    }
    print(paste0(which(names(msigdb)==i), " of ", length(names(msigdb))))
    
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


regionDB = loadRegionDB(file.path(datasets.dir,"MSigDB_hg19"))
saveRDS(regionDB,file.path(datasets.dir,"MSigDB_hg19_complete.rds") )
