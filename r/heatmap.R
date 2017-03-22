pkgs <- c('RColorBrewer', 'pvclust', 'gplots')
lapply(pkgs, require, character.only = TRUE)

palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)


# Load data
genre <- read.table("./fpkms.logs", header = FALSE, row.names = 1)
slogs <- read.table("./samples.logs")
colnames(genre) <- as.matrix(slogs)
genre <- as.matrix(genre)


## HIERARCHICAL AND BOOTSTRAP ANALYSIS
## by sample
rawdata <- t(genre)
scaledata=t(scale(genre))

## by genes
rawdata <- genre
scaledata=scale(genre)

## Clustering using dissimilarity analysis
# use "pairwise.complete.obs" when generating NAs
#hra <- hclust(as.dist(1-cor(t(scaledata), method="pearson",use = "pairwise.complete.obs")), method="average")
hra <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="ward.D2")
hca <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="ward.D2")


## CUT THE TREE
mycl <- cutree(hra, h=max(hra$height)/2)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
pdf("heatmap.1.pdf")
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=palette.green, scale="row", RowSideColors=mycolhc)
dev.off()

## BOOTSTRAPING to create pvalues
n=1000
a=0.95
pvData <- pvclust(scale(t(rawdata)), method.dist="correlation", method.hclust="ward.D2", nboot=n)

pdf("bootstrap.pdf")
plot(pvData, hang=-1)
pvrect(pvData, alpha=a)
dev.off()

## RETRIEVE MEMBERS OF SIGNIFICANT CLUSTERS.
clsig <- unlist(pvpick(pvData, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)

## Function to Color Dendrograms
dendroCol <- function(dend=dend, keys=keys, xPar="edgePar", bgr="red", fgr="blue", pch=20, lwd=1, ...) {
        if(is.leaf(dend)) {
                myattr <- attributes(dend)
                if(length(which(keys==myattr$label))==1){
                	attr(dend, xPar) <- c(myattr$edgePar, list(lab.col=fgr, col=fgr, pch=pch, lwd=lwd))
                	# attr(dend, xPar) <- c(attr$edgePar, list(lab.col=fgr, col=fgr, pch=pch))
                } else {
                	attr(dend, xPar) <- c(myattr$edgePar, list(lab.col=bgr, col=bgr, pch=pch, lwd=lwd))
                }
        }
  return(dend)
}

dend_colored <- dendrapply(as.dendrogram(pvData$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# use xPar="nodePar" to color tree labels


## PLOT HEATMAP
pdf("heatmap2.pdf")
heatmap.2(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=palette.green, scale="row", trace="none", RowSideColors = mycolhc, margins=c(8,20))
dev.off()

pdf("heatmap3.pdf")
heatmap.2(rawdata, Rowv=dend_colored, col=palette.red, symm=F,symkey=F,symbreaks=T, scale="row", trace="none", RowSideColors=mycolhc,margins=c(8,20))
dev.off()
