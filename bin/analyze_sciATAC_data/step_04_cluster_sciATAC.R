# find cell - types with Socrates ----

# requires R > 3.6.3: ml R/4.0.0-foss-2019b

## load libraries ----
library(Socrates)
library(harmony)

## create socrates object ----
obj <- readRDS("Atroot_filtered.merged.soc.rds")

## filter matrix ----
obj <- cleanData(obj, min.c=100, min.t=0.01, max.t=0.001, verbose=T, min.p=0.01, preclusterID="sampleID")

## normalize ----
obj <- tfidf(obj)

## reduce dimensions with SVD ----
obj <- reduceDims(obj, n.pcs=50, cor.max=0.8, verbose=T)

## reduce dimensions with UMAP ----
obj <- projectUMAP(obj, verbose=T, k.near=35, m.dist=0.05)

## estimate doublets ----
obj <- detectDoublets(obj, threads=30, n.pcs=50)

## cluster ----
obj <- callClusters(obj, 
                    res=0.6,
                    m.clst=50,
                    verbose=T, 
                    threshold=3, 
                    e.thresh=3, 
                    k.near=35)


################################################################################################### !
################################################################################################### !
################################################################################################### !

## plot results ----

## parameters ----
ptcex <- 0.25
opaque <- 1

## doublets ----
pdf("At_root_sci-ATAC_umap_doubletscore.pdf", width=8, height=8)
plotUMAP(obj, cluster_slotName = "meta", column = "doubletscore")
dev.off()

## clusters ----
pdf("At_root_sci-ATAC_umap_clusters.pdf", width=8, height=8)
plotUMAP(obj, cex=ptcex, opaque=opaque)
dev.off()

## libraries ----
pdf("At_root_sci-ATAC_umap_sampleID.pdf", width=8, height=8)
plotUMAP(obj, column="sampleID", cex=ptcex, opaque=opaque)
dev.off()

## read depth ----
pdf("At_root_sci-ATAC_umap_log10nSites.pdf", width=8, height=8)
plotUMAP(obj, column="log10nSites", cex=ptcex, opaque=opaque)
dev.off()

## output data ----
out.2 <- obj$Clusters
svd.2 <- obj$PCA
write.table(out.2, file="At_root_sciATAC.LouvainResults.txt", quote=F, row.names=T, col.names=T, sep="\t")
write.table(svd.2, file="At_root_sciATAC.SVD.txt", quote=F, row.names=T, col.names=T, sep="\t")
saveRDS(obj, file="At_root_sciATAC.socrates.results.rds")

