mkdir /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/EpigenotypeHM/random_BG/ 
#run script
Rscript /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/scripts/randomRegions.R /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/tracks/dmrs_sub_MethDiff_EpigenotypeHM.bed \
/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/EpigenotypeHM/random_BG \
14000

#sort bed files
cd /home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/EpigenotypeHM/random_BG \

conda activate betools
for file in `ls ./*bed`
do
    echo ${file}
    sort -k1,1 -k2,2n ${file} > ${file}.sorted.bed
done

#combine bed files 
cat *sorted.bed > random_bg_combined.bed
sort -k1,1 -k2,2n random_bg_combined.bed > random_bg_combined.sorted.bed
#and merge
bedtools merge -i random_bg_combined.sorted.bed > random_bg_merged.bed
sort -k1,1 -k2,2n random_bg_merged.bed > random_bg_merged.sorted.bed

######################################################################################
## this script allows to sample random sequences using a set of template sequences
## template sequences are provided as bed files
## fo each template sequence, the script samples N random locations in the same chromosome
## with equal size, then computes the frequency of a user given motif (e.g. CG)
## and keep the random sequences with similar motif frequency as the template one
## CDS regions are masked (N's) and regions with too many N's are rejected
## a target number of sequences is defined, so hopefully from the N sampled regions
## we should be able to extract the target number of sequences.
## if not, a warning is given
########################################################################################



args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  message("randomPeaks.R [bed.file] [out.dir] [min.seq] [genome.assembly]")
    exit()
  
}

#template.dir <- args[1] # contains the bed files to be used as templates
#bed.file <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/tracks/dmrs_sub_MethDiff_EpigenotypeHM.bed"
#out.dir <- "/home/heyj/icgc/dkfzlsdf/analysis/C010/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/EpigenotypeHM/random_BG " # main directory in which output will be written.
#min.seq <- 15000 # number of random sets

#template.dir <- args[1] # contains the bed files to be used as templates
bed.file <- args[1]
out.dir <- args[2] # main directory in which output will be written.
min.seq <- as.numeric(args[3]) # number of random sets
genome.assembly <- args[4]
message("=============================     ARGUMENTS             =========================================")
#message(paste0("template dir : ",template.dir))
message(paste0("bed file : ",bed.file))
message(paste0("output dir : ",out.dir))
message(paste0("number of randomizations : ",min.seq))
message("=================================================================================================")

randomPeaks <- function(bed.file,
                        genome=genome.assembly,
                        maskcds = FALSE,
                        min.seq=min.seq,
                        out.dir=out.dir) {
  
  if (genome =="mm10") {
    
    library(org.Mm.eg.db,quietly=TRUE,verbose=FALSE)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene,quietly=TRUE,verbose=FALSE)
    library(rtracklayer,quietly=TRUE,verbose=FALSE)
    library("BSgenome.Mmusculus.UCSC.mm10",quietly=TRUE,verbose=FALSE)
    ### define the genome of interest
    GENOME <- BSgenome.Mmusculus.UCSC.mm10
    TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  else if (genome =="hg19") {
    
    library(org.Hs.eg.db,quietly=TRUE,verbose=FALSE)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene,quietly=TRUE,verbose=FALSE)
    library(rtracklayer,quietly=TRUE,verbose=FALSE)
    library("BSgenome.Hsapiens.UCSC.hg19",quietly=TRUE,verbose=FALSE)
    ### define the genome of interest
    GENOME <- BSgenome.Hsapiens.UCSC.hg19
    TXDB <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  
  message(paste0("Genome : ",genome))
  #####====================================================================
  ########################      USER DEFINED   ############################
  #####====================================================================
  
  ### read template fasta files / bed files
  
  #main.dir <- '/ibios/raid6/users/herrmanc/Projects/LINACC/'
  #data.dir <- file.path(main.dir,'data/DMR/current')
  #out.dir <- file.path(main.dir,'output/shuffled_seq/stress')
  #   bed.file <- file.path(data.dir,'children','2014-06-10_DMRs_all_year0_smoking_8vs8.children.bed')
  #   bed.file <- file.path(data.dir,'mothers','2014-06-10_DMRs_all_year0_smoking_8vs8.mothers.bed')
  
  #my.bed <- import.bed(bed.file,format="bed")
  tmp <- read.table(bed.file)[,1:3]
  colnames(tmp) <- c("chr","start","end")
  my.bed <- makeGRangesFromDataFrame(tmp)
  l.my.bed <- length(my.bed)
  
  message(paste0("Number of template sequences : ",l.my.bed))
  
  ### define the target motif
  mot <- 'CG'
  
  message(paste0("here min.seq=",min.seq))
  N <- 100 * min.seq       ## number of random sampliwarnings
  eps <- 0.15       ## tolerance in the deviation in motif composition
  max.N <- 0.2     ## max. N's aloowed in sequence
  
  ##########################################################################
  ##########################################################################
  
  sl <- seqlengths(GENOME)
  
  ## select set of chromosomes to consider
  keep.chr <- seqnames(GENOME)[1:24]
  
  ## split BSgenome object into MaskedDNAString objects
  genome <- list()
  for (i in keep.chr) {genome[[i]] <- GENOME[[i]]}
  rm(GENOME)
  gc()

  ## decide whether to mask CDS with N's

  if (maskcds) {
      message('Masking CDS with Ns...')
      cds <- unlist(cdsBy(TXDB,by="gene"))
  
  ## create a mask with all coding sequences
    message("Building mask for cds...")
    cdsmask <- lapply(names(genome), function(x) {
      #message(x)
      X <- reduce(ranges(cds[seqnames(cds) == x]))
      Mask(mask.width=sl[x],
           start=start(X),
           width=width(X)
      )
    })
    names(cdsmask) <- names(genome)
    
    
    ## append CDS mask to already present masks
    maskedgenome <- list()
    for (x in names(cdsmask)) {
      message(x)
      names(cdsmask[[x]]) <- paste0("cds.",x)
      masks(genome[[x]]) <- append(masks(genome[[x]]),cdsmask[[x]])    
      maskedgenome[[x]] <- injectHardMask(genome[[x]],letter="N")
      maskedgenome[[x]] <- genome[[x]] ## do NOT mask CDS
    }   
  }
  maskedgenome <- genome
  rm(genome)
  gc()

  message("Genome is masked...")
  #######################################################################
  ##### NOW CREATE RANDOM VIEWS
  #######################################################################
  
  
  WARNINGS <- ""
  
  sets <- split(1:length(my.bed),1:length(my.bed) %% 10)
  
  random.seq <- lapply(sets, function(s) {
    message(paste0("new set of length ",length(s),"..."))
    mclapply(s, function(i) {
      
      #message(i)
      ## analyse the template (t) sequences
      t.chr <- as.character(seqnames(my.bed)[i])
      
      if (!t.chr %in% keep.chr) {return(NULL)}
      t.w <- width(my.bed)[i]  
      t.vi <- Views(maskedgenome[[t.chr]],start=start(my.bed)[i],width=t.w)
      
      ## determine motif composition
      t.mot <- countPattern(mot,t.vi)/width(t.vi)
      
      ## Now sample N random location in same chromosome
      
      l <- sl[t.chr]
      i.pos <- sample(1:(l-t.w-1),N)
      
      ## create random views
      r.vi <- Views(maskedgenome[[t.chr]],start=i.pos,width=t.w)
      
      ## remove sequences with too many N's
      r.N <- vcountPattern('N',r.vi)/t.w
      i.keep.N <- which(r.N < max.N)
      
      ## determine demotif composition
      r.mot <- vcountPattern(mot,r.vi)/t.w
      #r.mot <- vcountPattern(mot,r.vi)   ### to be used when absolute number of motifs is considered !!!!!!!!!!!
 
      #####################################################################################
      ### TWO DIFFERENT CRITERIA
      
      ## determine which sequences have a less that eps deviation from target composition
      i.keep.mot <- which(abs(r.mot-t.mot)/t.mot < eps)
      
      ## keep sequences with at least 3 CpGs
      #i.keep.mot <- which(r.mot >= 3)
      
      ####
      ######################################################################################
      i.keep <- intersect(i.keep.N,i.keep.mot)
      
      if (length(i.keep) < min.seq) { 
         message(paste0("Warning : template sequence ",i," : only ",length(i.keep)," available instead of ",min.seq,".."))
         WARNINGS <- paste0(WARNINGS,"\t","Warning : template sequence ",i," : only ",length(i.keep)," available instead of ",min.seq,"..")
      } else {
        i.keep <- i.keep[1:min.seq]
      }
      
      if (length(i.keep) >0) {
        return(GRanges(seqnames=t.chr,ranges=ranges(r.vi[i.keep]),motComp=r.mot[i.keep]))
      }
      gc()
      return(NULL)
    },mc.cores=8)
  })
  
  random.seq <- do.call("c",random.seq)
  
  #message("WARNINGS")
  warnings()
  ## now build min.seq sets of random sequences
  ## for those sequences for wich fewer random sequences could be determined
  ## sample over the existing ones
  
  random.seq <- random.seq[sapply(random.seq,length) > 0]
  
  X <- do.call("c",lapply(random.seq,function(x) {
    
    if (length(x) == min.seq) {
      i <- 1:length(x)
    } else {
      i <- sample(1:length(x),min.seq,replace=TRUE)
    }
    return(x[i])
  }))
  
  ## check for NULL elements
  
  i.null <- which(sapply(X,is.null))  
  if (length(i.null)>0) {X <- X[-i.null]}
  
  names(X) <- NULL
  X <- do.call("c",X)
  
  XX <- split(X,rep(1:min.seq,length(random.seq)))
  
  out.file <- sub(".bed","",basename(bed.file))
  
  dir.create(file.path(out.dir,out.file),showWarnings = FALSE)
  mclapply(1:length(XX),function(i) {  
    export(XX[[i]],file.path(out.dir,paste0(out.file,".shuffled.CpG.",i,".bed")))
  },mc.cores=6)
  saveRDS(XX,file=file.path(out.dir,out.file,paste0(out.file,'.ALL.shuffled.CpG.rds')))
  #  return(XX)
}

##########################################    MAIN    #########################################
#main.dir <- '/ibios/raid6/users/herrmanc/Projects/LINACC/'
#data.dir <- file.path(main.dir,'data/DMR/current')
#data.dir <- file.path(main.dir,'data/DMR/stress')

dir.create(out.dir, showWarnings = FALSE)

#files <- list.files(template.dir,pattern="bed",full.names=TRUE)

## NEW FILES FOR CHILDREN15/07/2014
#files <- list.files(file.path(data.dir,'children'),pattern=".*bed",recursive=TRUE,full.names=TRUE)


#randomPeaks(files[1],min.seq=min.seq)
#lapply(files,function(x) {message(paste0("+++++ ",x));randomPeaks(x,min.seq=min.seq,out.dir=out.dir)})

randomPeaks(bed.file,min.seq=min.seq,out.dir=out.dir)
