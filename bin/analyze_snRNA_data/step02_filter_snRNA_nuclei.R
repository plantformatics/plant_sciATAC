## clean cells ##

# load libraries
library(Matrix)
library(MASS)
library(viridis)
library(RColorBrewer)

# load functions
convert2Matrix <- function(x){
    x <- sparseMatrix(i=as.numeric(x$V1),
                      j=as.numeric(x$V2),
                      x=as.numeric(x$V3),
                      dimnames=list(levels(x$V1),levels(x$V2)))
    return(x)
}
plotDensity <- function(x, column="pMt", main=""){
    
    # get densities
    den <- kde2d(log10(x$total), x[,column], n=300)
    
    # plot
    image(den, col=colorRampPalette(c("white","grey75",rev(magma(18))))(100), useRaster=T, 
          bty="n", main=main)
    box()
}
filterCells <- function(x, column="pMt", threshold=0, direction="greater"){
    
    # estimate z-score
    x$zscore <- (x[,column]-mean(x[,column], na.rm=T))/sd(x[,column],na.rm=T)
    
    # greater than
    if(direction=="greater"){
        xx <- subset(x, x$zscore >= threshold)
        xx$zscore <- NULL
    }else{
        xx <- subset(x, x$zscore <= threshold)
        xx$zscore <- NULL
    }
    return(xx)
    
}

# load data
message(" - loading data ...")
meta <- read.table("Atroot_snRNAseq_merged.metadata.raw.txt")
all <- read.table("Atroot_snRNAseq_merged.sparse")

# convert to sparseMatrix
message(" - filtering barcodes ...")
all <- convert2Matrix(all)
ids <- intersect(rownames(meta), colnames(all))
all <- all[,ids]
meta <- meta[ids,]

# count number of genes per cell
meta$n.genes <- Matrix::colSums(all > 0)

# count number of transcripts per cell
meta$n.trx <- Matrix::colSums(all)

# filter cells
meta <- subset(meta, meta$total > 1000 & meta$n.genes >= 100 & meta$n.trx >= 500)

# prop mito, chloro, genic
meta$pMt <- meta$mt/meta$total
meta$pPt <- meta$pt/meta$total
meta$pNuc <- meta$nuclear/meta$total
meta$pTrx <- meta$n.trx/meta$total

# estimate densities
r1 <- subset(meta, meta$library=="Atroot_sn_rep1")
r2 <- subset(meta, meta$library=="Atroot_sn_rep2")
r3 <- subset(meta, meta$library=="Atroot_sn_rep3")
r4 <- subset(meta, meta$library=="Atroot_sn_rep4")
r5 <- subset(meta, meta$library=="Atroot_sn_rep5")

# verbose
message(" - # nuclei: rep 1 = ", nrow(r1))
message(" - # nuclei: rep 2 = ", nrow(r2))
message(" - # nuclei: rep 3 = ", nrow(r3))
message(" - # nuclei: rep 4 = ", nrow(r4))
message(" - # nuclei: rep 5 = ", nrow(r5))

# plot individually - mitocondrial proportion
pdf("Mitocondrial_distributions.pdf", width=8, height=6)
layout(matrix(c(1:6), ncol=3, nrow=2))
plotDensity(r1, main="rep1")
plotDensity(r2, main="rep2")
plotDensity(r3, main="rep3")
plotDensity(r4, main="rep4")
plotDensity(r5, main="rep5")
plot.new()
dev.off()

# chloroplast proportion
pdf("Chloroplast_distributions.pdf", width=8, height=6)
layout(matrix(c(1:6), ncol=3, nrow=2))
plotDensity(r1, column="pPt", main="rep1")
plotDensity(r2, column="pPt", main="rep2")
plotDensity(r3, column="pPt", main="rep3")
plotDensity(r4, column="pPt", main="rep4")
plotDensity(r5, column="pPt", main="rep5")
plot.new()
dev.off()

# nuclear proportion
pdf("Nuclear_distributions.pdf", width=8, height=6)
layout(matrix(c(1:6), ncol=3, nrow=2))
plotDensity(r1, column="pNuc", main="rep1")
plotDensity(r2, column="pNuc", main="rep2")
plotDensity(r3, column="pNuc", main="rep3")
plotDensity(r4, column="pNuc", main="rep4")
plotDensity(r5, column="pNuc", main="rep5")
dev.off()

# filtering individual
r1.f <- filterCells(r1, column="pMt", threshold=1, direction="less")
r2.f <- filterCells(r2, column="pMt", threshold=1, direction="less")
r3.f <- filterCells(r3, column="pMt", threshold=1, direction="less")
r4.f <- filterCells(r4, column="pMt", threshold=1, direction="less")
r5.f <- filterCells(r5, column="pMt", threshold=1, direction="less")

r1.f <- filterCells(r1.f, column="pPt", threshold=1, direction="less")
r2.f <- filterCells(r2.f, column="pPt", threshold=1, direction="less")
r3.f <- filterCells(r3.f, column="pPt", threshold=1, direction="less")
r4.f <- filterCells(r4.f, column="pPt", threshold=1, direction="less")
r5.f <- filterCells(r5.f, column="pPt", threshold=1, direction="less")

r1.f <- filterCells(r1.f, column="pNuc", threshold= -1, direction="greater")
r2.f <- filterCells(r2.f, column="pNuc", threshold= -1, direction="greater")
r3.f <- filterCells(r3.f, column="pNuc", threshold= -1, direction="greater")
r4.f <- filterCells(r4.f, column="pNuc", threshold= -1, direction="greater")
r5.f <- filterCells(r5.f, column="pNuc", threshold= -1, direction="greater")

# verbose
message(" - filtered: rep 1 = ", nrow(r1.f))
message(" - filtered: rep 2 = ", nrow(r2.f))
message(" - filtered: rep 3 = ", nrow(r3.f))
message(" - filtered: rep 4 = ", nrow(r4.f))
message(" - filtered: rep 5 = ", nrow(r5.f))

# specify keepers
meta$pass <- ifelse(rownames(meta) %in% rownames(r1.f), 1,
                    ifelse(rownames(meta) %in% rownames(r2.f), 1,
                           ifelse(rownames(meta) %in% rownames(r3.f), 1,
                                  ifelse(rownames(meta) %in% rownames(r4.f), 1, 
					 ifelse(rownames(meta) %in% rownames(r5.f), 1, 0)))))

# filter
message(" - final filter, outputting cleaned data ...")
f.meta <- meta[meta$pass > 0,]
all <- all[,rownames(f.meta)]
write.table(meta, file="Atroot_snRNAseq_merged.metadata.filtered.v2.txt", quote=F, row.names=T, col.names=T, sep="\t")
out <- as.data.frame(summary(all))
out$i <- rownames(all)[out$i]
out$j <- colnames(all)[out$j]
write.table(out, file="Atroot_snRNAseq_merged.v2.sparse", quote=F, row.names=F, col.names=F, sep="\t")
