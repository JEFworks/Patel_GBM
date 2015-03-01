########
## Recreating Figure 1c 
## Large copy number variation inference from single cell RNA-seq
## Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma
## Patel et al. Science, 2014
## @author: Jean Fan (jeanfan@fas.harvard.edu)
## @date: Oct 06, 2014
########


########
## Dependencies
########

library(biomaRt) ## for gene coordinates
mart.obj <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
library(caTools) ## for sliding window means


########
## Load Data
########

## load single cell data
## start with the FPM matrix derived from SCDE model, original paper used TPM
load('../RData/fpm.RData') 
## original paper used log2(TPM+1) transformation
fpm <- log2(fpm+1)
## get cell sample groups
cells1 <- colnames(fpm)[grepl('MGH26_', colnames(fpm))]
cells2 <- colnames(fpm)[grepl('MGH28_', colnames(fpm))]
cells3 <- colnames(fpm)[grepl('MGH29_', colnames(fpm))]
cells4 <- colnames(fpm)[grepl('MGH30_', colnames(fpm))]
cells5 <- colnames(fpm)[grepl('MGH31_', colnames(fpm))]

## load normal bulk data derived from normal brain samples in GTEx
## start with the FPM matrix derived from SCDE model, original paper used TPM
load('../RData/fpm_bulk.RData')
## original paper used log2(TPM+1) transformation
fpm.bulk <- log2(fpm.bulk+1)
 
 
########
## Filter Genes
########

## just make sure we have same genes in single cells and bulk 
genes.have <- intersect(rownames(fpm), rownames(fpm.bulk))
length(genes.have)
fpm <- fpm[genes.have,]
fpm.bulk <- fpm.bulk[genes.have,]

## original filtering proposed in paper 
## restrict to highly/confidently detected genes
## with average expression > 4.5 for all cells, or > 6 for one sample
vi <- rowMeans(fpm)>4.5; table(vi)
vi.c1 <- rowMeans(fpm[,cells1])>6; table(vi.c1)
vi.c2 <- rowMeans(fpm[,cells2])>6; table(vi.c2)
vi.c3 <- rowMeans(fpm[,cells3])>6; table(vi.c3)
vi.c4 <- rowMeans(fpm[,cells4])>6; table(vi.c4)
vi.c5 <- rowMeans(fpm[,cells5])>6; table(vi.c5)
vi.tot <- Reduce(union, list(rownames(fpm)[vi], rownames(fpm)[vi.c1], rownames(fpm)[vi.c2], rownames(fpm)[vi.c3], rownames(fpm)[vi.c4], rownames(fpm)[vi.c5]))
length(vi.tot)
## results in only ~2500 genes
## just try to get more comparable number of genes ~6000
vi <- rowMeans(fpm)>3; table(vi)
fpm <- fpm[vi,]
 
## also filter out poor cells
## original paper suggests filter out cells with less than 2000 detected genes
vi <- colSums(fpm>0)>2000; table(vi)
fpm <- fpm[,vi]
 
## apply same gene filtering to bulk
fpm.bulk <- fpm.bulk[rownames(fpm),]

dim(fpm)
dim(fpm.bulk)
 

########
## Normalize
########

## column mean normalize to control for cell quality
fpm <- t(t(fpm) - colMeans(fpm)) 
## row mean normalize also; ie. centering
ss.means <- rowMeans(fpm)
fpm <- fpm-ss.means
 
## column mean normalize to control for tumor quality
fpm.bulk <- t(t(fpm.bulk) - colMeans(fpm.bulk))
## paper says to subtract average normalized expression from single cell data
## similar to centering of the single cell data
fpm.bulk <- fpm.bulk-ss.means
 
## now we normalize to bulk
bulk.means <- rowMeans(fpm.bulk)
fpm.bulk <- fpm.bulk-bulk.means
fpm <- fpm-bulk.means
 

########
## Map genes to positions
########

## in paper, used start positions for genes
gos <- getBM(values=rownames(fpm),attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name","start_position","end_position"),filters=c("hgnc_symbol"),mart=mart.obj)
gos$pos <- gos$start_position
 
## organize into chromosomes
tl <- tapply(1:nrow(gos),as.factor(gos$chromosome_name),function(ii) {
    na.omit(fpm[gos$hgnc_symbol[ii[order(gos$pos[ii],decreasing=F)]],])
})
## only care about these chromosomes
tl <- tl[c(as.factor(1:22), "X")] 


########
## Plot
########

## color
pcol <- colorRampPalette(c("blue","white","red"),space="Lab")(300)
## label cells by sample
group <- unlist(lapply(colnames(fpm), function(x) strsplit(x, '[.]|_')[[1]][1]))
names(group) <- colnames(fpm)
sc <- rainbow(length(levels(as.factor(group))))
colcols <- sc[as.factor(group)]
names(colcols) <- colnames(fpm)

#png(file="../results/cnv.png",width=1000,height=600);
  ## size of box correspond to square root of size of chromosome 
  l <- layout(matrix(seq(1,length(tl)+1),1,length(tl)+1,byrow=T),heights=c(1),widths=sqrt(c(100, unlist(lapply(tl,nrow)))),FALSE)
  ## plot cell label colors first
  par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(1,1:length(group),t(as.matrix(as.numeric(as.factor(group)))),col=sc,xlab="",ylab="",axes=F); box()
  zlim <- c(-2,2);
  ## use sliding window average!!
  ## sliding window size of 100 genes
  window.size <- 101
  ## plot chromosomes
  lapply(names(tl),function(nam) {
    d <- tl[[nam]]
    d <- apply(d,2,runmean,k=window.size)
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
    image(1:nrow(d),1:ncol(d),d,col=pcol,zlim=zlim,xlab="",ylab="",axes=F,main=nam); box()
  })
#dev.off()

