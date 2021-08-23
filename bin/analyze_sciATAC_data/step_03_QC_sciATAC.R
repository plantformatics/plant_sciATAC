## process At root data ##

# load libraries
library(Socrates)

# input files
sci_rep1 <- "Arabidopsis_root_sciATAC_rep1.unique.tn5.bed"
sci_rep2 <- "Arabidopsis_root_sciATAC_rep2.unique.tn5.bed"

# set-up arguments
ann <- "TAIR10_GFF3_genes.gff"
chr <- "TAIR10_reference.sizes"

# output prefixes
pre1 <- "Atroot_sciATAC_rep1"
pre2 <- "Atroot_sciATAC_rep2"

# load objects
obj1 <- loadBEDandGenomeData(sci_rep1, ann, chr)
obj2 <- loadBEDandGenomeData(sci_rep2, ann, chr)

# call ACRs
obj1 <- callACRs(obj1, genomesize=95e6, 
                 shift= -50, 
                 extsize=100,
                 fdr=0.05,
                 output="sci1_peaks", 
                 tempdir="./sci1_macs2_temp", 
                 verbose=T)

obj2 <- callACRs(obj2, genomesize=95e6,
                 shift= -50,
                 extsize=100,
                 fdr=0.05,
                 output="sci2_peaks",
                 tempdir="./sci2_macs2_temp",
                 verbose=T)

# build metadata
obj1 <- buildMetaData(obj1,
                      tss.window=2000, 
                      verbose=TRUE)
obj2 <- buildMetaData(obj2,
                      tss.window=2000, 
                      verbose=TRUE)

# filter
message(" - finding cells ...")
obj1 <- findCells(obj1, 
                  doplot=T,
                  max.cells=8000,
                  min.tn5=1000,
                  filt.tss=TRUE, 
                  filt.frip=TRUE,
                  frip.min.freq=0.5,
                  prefix=pre1)
obj2 <- findCells(obj2,
                  doplot=T,
                  max.cells=8000,
                  min.tn5=1000,
                  filt.tss=TRUE,
                  filt.frip=TRUE,
                  frip.min.freq=0.5,
                  prefix=pre2)

# generate sparse matrix
obj1 <- generateMatrix(obj1,
                       filtered=T, 
                       windows=500, 
                       peaks=F, 
                       verbose=T)
obj2 <- generateMatrix(obj2,
                       filtered=T, 
                       windows=500, 
                       peaks=F, 
                       verbose=T)

# save QC data
saveRDS(obj1, file=paste0(pre1, ".QC.rds"))
saveRDS(obj2, file=paste0(pre2, ".QC.rds"))

# convert to Socrates format for downstream analysis. 
soc.obj1 <- convertSparseData(obj1, verbose=T)
soc.obj2 <- convertSparseData(obj2, verbose=T)

# merge objects
soc.list <- list(sci_rep1=soc.obj1,
		 sci_rep2=soc.obj2)

# merge
merge.obj <- mergeSocratesRDS(obj.list=soc.list)

# save QC object
saveRDS(merge.obj, file="Atroot_filtered.merged.soc.rds")
