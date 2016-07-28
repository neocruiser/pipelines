pkgs <- c('RColorBrewer', 'pvclust', 'gplots')
lapply(pkgs, require, character.only = TRUE)

display.brewer.pal(9, "YlGnBu")
palett <- brewer.pal(9,"Greens")

# Load data
genre <- read.table("./fpkms.logs", header = FALSE, row.names = 1)
slogs <- read.table("./samples.logs")
colnames(genre) <- as.matrix(slogs)
genre <- as.matrix(genre)
#heatmap.g <- heatmap(genre, Rowv=NA, Colv=NA, col=palett, margins=c(1,40), scale="row", labCol=c(seq(1:16)))


## HIERARCHICAL AND BOOTSTRAP ANALYSIS
# Color function to generate green-red heat maps
my.colorFct <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
	if (n < 2) stop("n must be greater than 2")
	n1 <- n%/%2
	n2 <- n - n1
	c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}

## by sample
rawdata <- t(genre)
scaledata=t(scale(genre))

## by genes
rawdata <- genre
scaledata=scale(genre)


hra <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete")
hca <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row")

## CUT THE TREE
mycl <- cutree(hra, h=max(hra$height)/2)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## BOOTSTRAPING
n=5000
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
heatmap(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## PLOT HEATMAP
#x11(height=5,width =8)
pdf("heatmap.pdf")
heatmap.2(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc,margins=c(8,20))
dev.off()
