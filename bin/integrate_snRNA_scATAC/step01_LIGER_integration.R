###################################################################################################
##                          scATAC/snRNA-seq integration with LIGER                              ##
###################################################################################################

# load arguments ----------------------------------------------------------------------------------
arg <- commandArgs(T)
if(length(arg)!=1){stop("Rscript LIGER_integration_arabidopsis.R <prefix>")}
input <- as.character(arg[1])


# load libraries ----------------------------------------------------------------------------------
library(liger)
library(parallel)
library(Matrix)
library(presto)


# load functions ----------------------------------------------------------------------------------
loadData <- function(x){
    
    # load RNA data -------------------------------------------------------------------------------
    message(" - loading and processing snRNA-seq data ...")
    a <- read.table("Atroot_snRNAseq_merged.sparse")
    a <- sparseMatrix(i=as.numeric(a$V1),j=as.numeric(a$V2),x=as.numeric(a$V3),dimnames=list(levels(a$V1),levels(a$V2)))
    b <- read.table("Atroot_snRNAseq_merged.metadata.filtered.txt")
    genes <- read.table("TAIR10_geneAnnotation.bed")
    
    # specify root data
    dat <- x
    
    # harmonize inputs
    shared.cells <- intersect(rownames(b), colnames(a))
    b <- b[shared.cells,]
    a <- a[,shared.cells]
    a <- a[rownames(a) %in% as.character(genes$V4),]
    a <- a[Matrix::rowSums(a) > 0,]
    a <- a[,Matrix::colSums(a) > 0]
    b <- b[colnames(a),]
    message("   * number of cells = ",ncol(a))
    
    
    # load ATAC data ------------------------------------------------------------------------------
    message(" - loading and processing scATAC-seq data ...")
    atac.meta <- read.table("At_root_sciATAC.LouvainResults.txt")
    
    # load gene activity
    activity.matrix <- read.table("At_root.genes.sparse")
    activity.matrix <- sparseMatrix(i=as.numeric(activity.matrix$V1),
                                    j=as.numeric(activity.matrix$V2),
                                    x=as.numeric(activity.matrix$V3),
                                    dimnames=list(levels(activity.matrix$V1),levels(activity.matrix$V2)))
    activity.matrix <- activity.matrix[Matrix::rowSums(activity.matrix)>0,]
    activity.matrix <- activity.matrix[,Matrix::colSums(activity.matrix)>0]
    
    # select sites
    ids <- intersect(rownames(atac.meta), colnames(activity.matrix))
    atac.meta <- atac.meta[ids,]
    activity.matrix <- activity.matrix[,ids]
    
    # filter empty features
    activity.matrix <- activity.matrix[Matrix::rowSums(activity.matrix>0)>0,]
    activity.matrix <- activity.matrix[,Matrix::colSums(activity.matrix)>0]
    atac.meta <- atac.meta[colnames(activity.matrix),]
    message("   * number of cells = ",ncol(activity.matrix))
    
    # shared genes between data sets --------------------------------------------------------------
    shared.genes <- intersect(rownames(a), rownames(activity.matrix))
    activity.matrix <- activity.matrix[shared.genes,]
    a <- a[shared.genes,]
    a <- a[,Matrix::colSums(a)>0]
    activity.matrix <- activity.matrix[,Matrix::colSums(activity.matrix) > 0]
    atac.meta <- atac.meta[colnames(activity.matrix),]
    b <- b[colnames(a),]
    
    
    # report number of cells in each
    message(" - scRNAseq = ", ncol(a), " | scATAC-seq = ", ncol(activity.matrix))
    
    # create LIGER object RNA --------------------------------------------------------------------
    both <- createLiger(list(rna=a, atac=activity.matrix))
    
    # return
    return(both)
    
}
extractMeta <- function(x){
    
    # original variables
    rna <- read.table("Atroot_snRNAseq_merged.metadata.filtered.txt")
    atac <- read.table("At_root_sciATAC.LouvainResults.txt")
    
    # extract variables
    cd <- x@cell.data
    ump <- x@tsne.coords
    colnames(ump) <- c("umap1", "umap2")
    cl <- data.frame(iNMF_clusters=x@clusters, cellIDs=names(x@clusters))
    share <- intersect(rownames(cd), rownames(ump))
    share <- intersect(share, rownames(cl))
    cd <- cd[share,]
    ump <- ump[share,]
    cl <- cl[share,]
    
    # create integrated df
    df <- as.data.frame(cbind(cd,ump,cl))
    
    # get rna data
    #rna <- rna[rownames(rna) %in% rownames(df),]
    #atac <- atac[rownames(atac) %in% rownames(df),]
    #celltypes <- atac$celltype
    #names(celltypes) <- rownames(atac)
    #df$celltype <- celltypes[rownames(df)]
    
    # get cell-type predictions based on scATAC annotations
    #props <- table(df$iNMF_clusters, df$celltype)
    #props <- props[,!colnames(props) %in% c("unknown")]
    #props <- prop.table(props, 1)
    #props <- t(apply(props, 1, function(z){z/max(z)}))
    #calls <- apply(props, 1, function(z){names(z)[which.max(z)]})
    #if(class(calls) == "list"){
    #    calls <- do.call(c, calls)
    #}
    #df$celltype_extended <- calls[as.character(df$iNMF_clusters)]
    #print(str(df))
    
    # return
    return(df)
    
}
loadPeaks <- function(x){
    
    # load leaf or root peaks
    message(" - loading scATAC peak x cell binarized data ...")
    atac <- readRDS("At_root.merged.rds")
    peak.names <- rownames(atac)
    peak.names <- as.data.frame(do.call(rbind, strsplit(peak.names, "_")))
    peak.names <- paste0(peak.names$V1, ':', peak.names$V2, '-', peak.names$V3)
    rownames(atac) <- peak.names
    
    # create new liger object
    message(" - loading peak data into LIGER object ...")
    x.ds <- x
    atac <- atac[ ,intersect(colnames(atac),colnames(x@raw.data[['atac']]))]
    atac <- atac[Matrix::rowSums(atac>0)>(ncol(atac)*0),]
    atac <- atac[ ,intersect(colnames(atac),colnames(x@raw.data[['atac']]))]
    atac <- atac[,colnames(atac) %in% rownames(x@H.norm)]
    x.ds@raw.data[['atac']] <- atac
    x.ds <- normalize(x.ds)
    message(" - added ", nrow(atac), " peaks ...")
    
    # return data
    return(x.ds)
}
selectMarkers <- function(a, b){
    
    # subset by markers
    a <- a[as.character(a$feature) %in% as.character(b$geneID),]
    a <- a[order(-rank(a$group), rank(a$auc), decreasing=T),]
    
    # attach cell-type info and gene name
    gname <- as.character(b$name)
    names(gname) <- b$geneID
    type <- as.character(b$type)
    names(type) <- b$geneID
    
    # append
    a$name <- gname[as.character(a$feature)]
    a$celltype <- type[as.character(a$feature)]
    
    # return
    return(a)
    
}
runDifAccess <- function(a, threads=1){
    
    # iterate over peaks
    clusts1 <- a@clusters[names(a@clusters) %in% colnames(a@norm.data[['atac']])]
    clusts1 <- droplevels(clusts1)
    clusts <- split(names(clusts1), clusts1)
    dat <- a@norm.data[["atac"]][,names(clusts1)]
    dat <- dat[Matrix::rowSums(dat>0) > (ncol(dat)*0.01),]
    ids <- colnames(dat)
    row.num <- nrow(dat)
    set.seed(1111)
    control.ids <- lapply(clusts, function(z){
        sample(z, 42)
    })
    control.ids <- do.call(c, control.ids)
    
    # run
    out <- mclapply(names(clusts), function(x){
        out.1 <- lapply(seq(1:row.num), function(z){
            if((z %% 1000)==0){message(" - iterated over ",z," ACRs for cluster ",x, " ...")}
            return(wilcox.test(dat[z,clusts[[x]]], dat[z,control.ids])$p.value)
        })
        names(out.1) <- rownames(dat)
        df <- data.frame(peak=rownames(dat),
                         group=rep(x, nrow(dat)),
                         pval=do.call(c, out.1))
        return(df)
    }, mc.cores=threads)
    out <- as.data.frame(do.call(rbind, out))
    out$padj <- p.adjust(out$pval, method="fdr")
    return(out)
    
}
runWilcoxon <- function(object, data.use = "all", compare.method){
    
    # check parameter inputs
    if (missing(compare.method)) {
        stop("Parameter *compare.method* cannot be empty!")
    }
    if (compare.method != "datasets" & compare.method != "clusters") {
        stop("Parameter *compare.method* should be either *clusters* or *datasets*.")
    }
    if (compare.method == "datasets") {
        if (length(names(object@norm.data)) < 2) {
            stop("Should have at least TWO inputs to compare between datasets")
        }
        if (!missing(data.use) & length(data.use) < 2) {
            stop("Should have at least TWO inputs to compare between datasets")
        }
    }
    
    ### create feature x sample matrix
    if (data.use == "all" | length(data.use) > 1) { # at least two datasets
        if (data.use == "all") {
            print(paste0("Performing Wilcoxon test on ALL datasets: ", toString(names(object@norm.data))))
            sample.list <- attributes(object@norm.data)$names
        }
        else {
            print(paste0("Performing Wilcoxon test on GIVEN datasets: ", toString(data.use)))
            sample.list <- data.use
        }
        genes <- Reduce(intersect, lapply(sample.list, function(sample) {
            object@norm.data[[sample]]@Dimnames[[1]]
        })) # get all shared genes of every datasets
        
        feature_matrix <- Reduce(cbind, lapply(sample.list, function(sample) {
            object@norm.data[[sample]][genes, ]
        })) # get feature matrix, shared genes as rows and all barcodes as columns
        
        # get labels of clusters and datasets
        cell_source <- object@cell.data[["dataset"]] # from which dataset
        names(cell_source) <- names(object@clusters)
        cell_source <- cell_source[colnames(feature_matrix), drop = TRUE]
        clusters <- object@clusters[colnames(feature_matrix), drop = TRUE] # from which cluster
    } else { # for one dataset only
        print(paste0("Performing Wilcoxon test on GIVEN dataset: ", data.use))
        feature_matrix <- object@norm.data[[data.use]]
        clusters <- object@clusters[object@norm.data[[data.use]]@Dimnames[[2]], drop = TRUE] # from which cluster
    }
    
    ### perform wilcoxon test
    if (compare.method == "clusters") { # compare between clusters across datasets
        len <- nrow(feature_matrix)
        if (len > 50000) {
            print("Calculating Large-scale Input...")
            results <- Reduce(rbind, lapply(suppressWarnings(split(seq(len), seq(len / 50000))), function(index) {
                wilcoxauc(log(feature_matrix[index, ] + 1e-10), clusters)
            }))
        } else {
            results <- wilcoxauc(log(feature_matrix + 1e-10), clusters)
        }
    }
    
    if (compare.method == "datasets") { # compare between datasets within each cluster
        results <- Reduce(rbind, lapply(levels(clusters), function(cluster) {
            sub_barcodes <- names(clusters[clusters == cluster]) # every barcode within this cluster
            sub_label <- paste0(cluster, "-", cell_source[sub_barcodes]) # data source for each cell
            sub_matrix <- feature_matrix[, sub_barcodes]
            if (length(unique(cell_source[sub_barcodes])) == 1) { # if cluster has only 1 data source
                print(paste0("Note: Skip Cluster ", cluster, " since it has only ONE data source."))
                return()
            }
            return(wilcoxauc(log(sub_matrix + 1e-10), sub_label))
        }))
    }
    
    # return
    return(results)
}
plotMarkers <- function(both, markers, input){
    
    # iterate over marker genes
    message(" - saving marker gene profiles to list ...")
    outs <- lapply(unique(markers$geneID), function(x){
        plotGene(both, x, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
    })
    
    # plot RNA
    pdf(paste0(input,".RNA_markers.pdf"), width=24, height=264)
}
getFDR <- function(observed, expected, fdr=0.01, grid=1000, verbose=F){
    
    # convert to data.frames
    obs <- as.data.frame(summary(observed))
    obs$i <- rownames(observed)[obs$i]
    obs$j <- colnames(observed)[obs$j]
    exp <- as.data.frame(summary(expected))
    exp$i <- rownames(expected)[exp$i]
    exp$j <- colnames(expected)[exp$j]
    
    # split into +/-
    n.exp <- subset(exp, exp$x < 0)
    p.exp <- subset(exp, exp$x > 0)
    
    # get counts
    pos.nexp <- nrow(p.exp)
    neg.nexp <- nrow(n.exp)
    pos.nobs <- nrow(subset(obs, obs$x > 0))
    neg.nobs <- nrow(subset(obs, obs$x < 0))
    if(verbose){message(" - number expected links = (+) ",pos.nexp, " | (-) ",neg.nexp)}
    if(verbose){message(" - number observed links = (+) ",pos.nobs, " | (-) ",neg.nobs)}
    
    # generate range of thresholds
    p.vals <- seq(from=0, to=1, length.out=grid)
    n.vals <- seq(from=0, to= -1, length.out=grid)
    
    # iterate over grid
    if(verbose){message(" - scanning positive thresholds ...")}
    p.thresh <- c()
    for(i in p.vals){
        num.exp <- sum(p.exp$x > as.numeric(i))
        c.fdr <- num.exp/pos.nexp
        if(is.na(c.fdr)){
            c.fdr <- 0
        }
        p.thresh <- c(p.thresh, c.fdr)
        message(" - (+) correlation threshold = ", i, " | FDR = ", c.fdr)
    }
    if(verbose){message(" - scanning negative thresholds ...")}
    n.thresh <- c()
    for(i in n.vals){
        num.exp <- sum(n.exp$x < as.numeric(i))
        c.fdr <- num.exp/(neg.nexp)
        if(is.na(c.fdr)){
            c.fdr <- 0
        }
        n.thresh <- c(n.thresh, c.fdr)
        message(" - (-) correlation threshold = ", i, " | FDR = ", c.fdr)
    }
    
    # select cut-offs
    p.threshold <- min(p.vals[which(p.thresh <= fdr)])
    n.threshold <- max(n.vals[which(n.thresh <= fdr)])
    
    # filter
    obs <- subset(obs, obs$x > p.threshold | obs$x < n.threshold)
    
    # verbose number of +/- linkages
    pos.links <- nrow(subset(obs, obs$x > 0))
    neg.links <- nrow(subset(obs, obs$x < 0))
    message(" - found ",pos.links, " + and ", neg.links," - gene-peak links ...")
    
    # return
    return(obs)
}
linkGenesAndPeaks2 <- function (gene_counts, peak_counts, genes.list = NULL, dist = "spearman",
                                alpha = 0.05, max_dist=50000, path_to_coords){
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package \"GenomicRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('GenomicRanges')", call. = FALSE)
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Package \"IRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('IRanges')", call. = FALSE)
    }
    peak.names <- strsplit(rownames(peak_counts), "[:-]")
    chrs <- Reduce(append, lapply(peak.names, function(peak) {
        peak[1]
    }))
    chrs.start <- Reduce(append, lapply(peak.names, function(peak) {
        peak[2]
    }))
    chrs.end <- Reduce(append, lapply(peak.names, function(peak) {
        peak[3]
    }))
    peaks.pos <- GenomicRanges::GRanges(seqnames = chrs, ranges = IRanges::IRanges(as.numeric(chrs.start),
                                                                                   end = as.numeric(chrs.end)))
    gene.names <- read.csv2(path_to_coords, sep = "\t", header = FALSE,
                            stringsAsFactors = F)
    gene.names <- gene.names[complete.cases(gene.names), ]
    genes.coords <- GenomicRanges::GRanges(seqnames = gene.names$V1,
                                           ranges = IRanges::IRanges(as.numeric(gene.names$V2),
                                                                     end = as.numeric(gene.names$V3)))
    names(genes.coords) <- gene.names$V4
    gene_counts <- t(gene_counts)
    peak_counts <- t(peak_counts)
    if (missing(genes.list)) {
        genes.list <- colnames(gene_counts)
    }
    missing_genes <- !genes.list %in% names(genes.coords)
    if (sum(missing_genes) != 0) {
        print(paste0("Removing ", sum(missing_genes), " genes not found in given gene coordinates..."))
    }
    genes.list <- genes.list[!missing_genes]
    genes.coords <- genes.coords[genes.list]
    print("Calculating correlation for gene-peak pairs...")
    each.len <- 0
    elements <- lapply(seq(length(genes.list)), function(pos) {
        gene.use <- genes.list[pos]
        gene.loci <- GenomicRanges::trim(suppressWarnings(GenomicRanges::promoters(GenomicRanges::resize(genes.coords[gene.use],
                                                                                                         width = 1, fix = "start"), upstream = max_dist, downstream = max_dist)))
        peaks.use <- S4Vectors::queryHits(GenomicRanges::findOverlaps(peaks.pos,
                                                                      gene.loci))
        if ((x <- length(peaks.use)) == 0L) {
            return(list(NULL, as.numeric(each.len), NULL))
        }
        res <- suppressWarnings(psych::corr.test(x = gene_counts[,
                                                                 gene.use], y = as.matrix(peak_counts[, peaks.use]),
                                                 method = dist, adjust = "holm", ci = FALSE, use = "complete"))
        pick <- res[["p"]] < alpha
        pick[is.na(pick)] <- FALSE
        if (sum(pick) == 0) {
            return(list(NULL, as.numeric(each.len), NULL))
        }
        else {
            res.corr <- as.numeric(res[["r"]][pick])
            peaks.use <- peaks.use[pick]
        }
        assign("each.len", each.len + length(peaks.use), envir = parent.frame(2))
        return(list(as.numeric(peaks.use), as.numeric(each.len),
                    res.corr))
    })
    i_index <- Reduce(append, lapply(elements, function(ele) {
        ele[[1]]
    }))
    p_index <- c(0, Reduce(append, lapply(elements, function(ele) {
        ele[[2]]
    })))
    value_list <- Reduce(append, lapply(elements, function(ele) {
        ele[[3]]
    }))
    regnet <- sparseMatrix(i = i_index, p = p_index, x = value_list,
                           dims = c(ncol(peak_counts), length(genes.list)), dimnames = list(colnames(peak_counts),
                                                                                            genes.list))
    return(regnet)
}

# load data ---------------------------------------------------------------------------------------
both <- loadData(input)
markers <- read.table("markers.bed", header=T)

# presets
res <- 0.3
docenter <- F
k.num <- 25
fc.thresh <- 0.5
lambda <- 5


# analysis ----------------------------------------------------------------------------------------

# process
message(" - normalize, select, and scale ...")
both <- normalize(both)
both <- selectGenes(both, datasets.use=1)
both <- scaleNotCenter(both)

# joint matrix factorization
message(" - run iNMF ...")
both <- optimizeALS(both, k = k.num, lambda = lambda)

# quantile normalize
both <- quantile_norm(both, do.center=docenter, ref_dataset=1)

# cluster
both <- louvainCluster(both, resolution=res, k=35, eps=0, prune=1/10)
both <- runUMAP(both, distance = 'euclidean', n_neighbors = 35, min_dist = 0.01) # euclidean

# plot
pdf(paste0(input,".LIGER_integration.pdf"), width=8, height=7)
embed <- plotByDatasetAndCluster(both,
                                 axis.labels = c('UMAP_1', 'UMAP_2'),
                                 return.plots=T,
                                 pt.size=0.1)
embed[[1]]
dev.off()

pdf(paste0(input,".LIGER_clusters.pdf"), width=8, height=7)
embed[[2]]
dev.off()

# output results
df <- extractMeta(both)
inmf <- both@H.norm
colnames(inmf) <- paste0("NMF_",seq(1:ncol(inmf)))
write.table(df, file=paste0(input,"_iNMF_integration_meta.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(inmf, file=paste0(input,"_iNMF.norm.txt"), quote=F, row.names=T, col.names=T, sep='\t')
saveRDS(both, paste0(input,".iNMF_liger.rds"))

# save normalized data
atac.norm <- as.data.frame(summary(both@norm.data[['atac']]))
rna.norm <- as.data.frame(summary(both@norm.data[['rna']]))
atac.norm$i <- rownames(both@norm.data[['atac']])[atac.norm$i]
atac.norm$j <- colnames(both@norm.data[['atac']])[atac.norm$j]
rna.norm$i <- rownames(both@norm.data[['rna']])[rna.norm$i]
rna.norm$j <- colnames(both@norm.data[['rna']])[rna.norm$j]
write.table(atac.norm, file=paste0(input, "_ATAC_norm.sparse"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(rna.norm, file=paste0(input, "_RNA_norm.sparse"), quote=F, row.names=F, col.names=F, sep="\t")

# differentially accessible and expressed markers
dg <- runWilcoxon(both, data.use = 'all', compare.method = 'clusters')
write.table(dg, file=paste0(input,"_ALL_wilcoxon.txt"), quote=F, row.names=F, col.names=T, sep="\t")
dg.r <- runWilcoxon(both, data.use = 'rna', compare.method = 'clusters')
write.table(dg.r, file=paste0(input,"_RNA_wilcoxon.txt"), quote=F, row.names=F, col.names=T, sep="\t")
dg.a <- runWilcoxon(both, data.use = 'atac', compare.method = 'clusters')
write.table(dg.a, file=paste0(input,"_ATAC_wilcoxon.txt"), quote=F, row.names=F, col.names=T, sep="\t")

# get marker results
marker.dg <- selectMarkers(dg, markers)
write.table(marker.dg, file=paste0(input,"_MARKERS_ALL_wilcoxon.txt"), quote=F, row.names=F, col.names=T, sep="\t")
marker.dg.r <- selectMarkers(dg.r, markers)
write.table(marker.dg.r, file=paste0(input,"_MARKERS_RNA_wilcoxon.txt"), quote=F, row.names=F, col.names=T, sep="\t")
marker.dg.a <- selectMarkers(dg.a, markers)
write.table(marker.dg.a, file=paste0(input,"_MARKERS_ATAC_wilcoxon.txt"), quote=F, row.names=F, col.names=T, sep="\t")

# load peak data
both.x <- loadPeaks(both)

# run diff peaks
peak.wilcoxon <- runWilcoxon(both.x, data.use = 2, compare.method = 'clusters')
peak.wilcoxon$padj <- p.adjust(peak.wilcoxon$pval, method="BH")
peak.wilcoxon <- peak.wilcoxon[peak.wilcoxon$padj < 0.05,]
peak.wilcoxon <- peak.wilcoxon[peak.wilcoxon$logFC > fc.thresh,]
peaks.sel <- unique(peak.wilcoxon$feature)
both.x@raw.data[['atac']] <- both.x@raw.data[['atac']][peaks.sel, ]
message(" - number of significant peaks = ", length(peaks.sel))

# impute peak accessibility on snRNA-seq data
both.x1 <- imputeKNN(both.x, reference = 'atac')
both.x2 <- imputeKNN(both.x, reference = 'rna')

# extract peak/gene matrices for gene-peak linkages
gmat <- cbind(both.x2@norm.data[['rna']], both.x2@norm.data[['atac']])
pmat <- cbind(both.x1@norm.data[['rna']], both.x1@norm.data[['atac']])

# shuffle
sgmat <- as.data.frame(summary(gmat))
spmat <- as.data.frame(summary(pmat))

# shuffled
sgmat$i <- sgmat$i[sample(length(sgmat$i))]
sgmat$j <- sgmat$j[sample(length(sgmat$j))]
spmat$i <- spmat$i[sample(length(spmat$i))]
spmat$j <- spmat$j[sample(length(spmat$j))]
sgmat <- sparseMatrix(i=sgmat$i, j=sgmat$j, x=sgmat$x, dimnames=list(rownames(gmat), colnames(gmat)))
spmat <- sparseMatrix(i=spmat$i, j=spmat$j, x=spmat$x, dimnames=list(rownames(pmat), colnames(pmat)))

# scan for gene-peak links
gene.bed <- "TAIR10_geneAnnotation.bed"
message(" - getting gene-peak linkages ...")
regnet <- linkGenesAndPeaks2(gene_counts = gmat,
                            peak_counts = pmat,
                            dist = 'spearman',
                            alpha = 0.05,
                            path_to_coords = gene.bed)

# scan for gene-peak links shuffled
message(" - getting shuffled gene-peak linkages ...")
shufregnet <- linkGenesAndPeaks2(gene_counts=sgmat,
                                peak_counts=spmat,
                                dist="spearman",
                                alpha=0.05,
                                path_to_coords=gene.bed)

# save obs and exp gene-peak links
obs <- as.data.frame(summary(regnet))
obs$i <- rownames(regnet)[obs$i]
obs$j <- colnames(regnet)[obs$j]
exp <- as.data.frame(summary(shufregnet))
exp$i <- rownames(shufregnet)[exp$i]
exp$j <- colnames(shufregnet)[exp$j]
write.table(obs, file=paste0(input,".observed.gene-peak-links.txt"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(exp, file=paste0(input,".expected.gene-peak-links.txt"), quote=F, row.names=F, col.names=F, sep="\t")

# find FDR < 0.01
gpl_fdr01 <- getFDR(regnet, shufregnet, fdr=0.05, verbose=T)

# filter
write.table(gpl_fdr01, file=paste0(input,"_peak_gene_links.txt"), quote=F, row.names=F, col.names=F, sep="\t")

# tidy memory and save
rm(sgmat, spmat, regnet, shufregnet, gpl_fdr01)

# save imputed data
impute.gene <- as.data.frame(summary(gmat))
impute.peak <- as.data.frame(summary(pmat))
impute.gene$i <- rownames(gmat)[impute.gene$i]
impute.gene$j <- colnames(gmat)[impute.gene$j]
impute.peak$i <- rownames(pmat)[impute.peak$i]
impute.peak$j <- colnames(pmat)[impute.peak$j]
write.table(impute.gene, file=paste0(input,".coembed_RNA.sparse"), quote=F, row.names=F, col.names=F, sep="\t")
write.table(impute.peak, file=paste0(input,".coembed_ATAC.sparse"), quote=F, row.names=F, col.names=F, sep="\t")

# save interaction track
makeInteractTrack(regnet, path_to_coords = 'TAIR10_geneAnnotation.bed')

# exit
message(" - finished LIGER analysis ...")



