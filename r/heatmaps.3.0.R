pkgs <- c('RColorBrewer', 'pvclust', 'gplots', 'vegan')
lapply(pkgs, require, character.only = TRUE)

# choose color palettes
display.brewer.all()
palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)


# Load data in matrix form
genre <- as.matrix(read.table("expressions", header = TRUE, row.names = 1))


## debugging
## resampling
tenpercent <- c(dim(genre)[1] * .1)
selected <- sample(dim(genre)[1], tenpercent)
rawdata <- genre[selected, ]
rawdata <- decostand(x = rawdata, method = s)
scaledata=scale(rawdata)



standardize_df <- c("standardize", "range", "log")
normalize_df <- c("complete", "ward.D2", "average")
correlate_rows <- c("pearson", "spearman")
correlate_columns <- c("pearson", "spearman")


for ( s in standardize_df ) {
    for ( n in normalize_df ) {
        for ( cr in correlate_rows ) {
            for ( cc in correlate_columns ) {

                # standardization
                genre <- decostand(x = genre, method = s)
                #genre <- wisconsin(genre)

                ## HIERARCHICAL AND BOOTSTRAP ANALYSIS
                ## clustering by sample (columns)
                rawdata <- t(genre)
                scaledata=t(scale(genre))

                ## clustering by genes/species (rows)
                rawdata <- genre
                scaledata=scale(genre)

                ## Clustering using dissimilarity analysis
                # use "pairwise.complete.obs" when generating NAs
                #hra <- hclust(as.dist(1-cor(t(scaledata), method="pearson",use = "pairwise.complete.obs")), method="average")
                hra <- hclust(as.dist(1-cor(t(scaledata), method= cr)), method= n)
                hca <- hclust(as.dist(1-cor(scaledata, method= cc)), method= n)


                ## CUT THE TREE
                mycl <- cutree(hra, h=max(hra$height)/2)
                maxClusters <- length(unique(mycl))
                if ( maxClusters <= 12) {
                    mycolhc <- brewer.pal(maxClusters, name = 'Set3')
                    
                } else {
                    mycolhc <- colorRampPalette(brewer.pal(11, name="Spectral"))(maxClusters)
                }
                mycolhc <- mycolhc[as.vector(mycl)]
                
                pdf(paste("heatmap1.STD",s,".CLU",n,".VAR-CORR",cr,".FEA-CORR",cc,".pdf", sep = ""))
                heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=palette.green, scale="row", RowSideColors=mycolhc, cexRow=.1, cexCol=.1)
                dev.off()



                
                ## BOOTSTRAPING to create pvalues
                # multiscale bootstrap resampling
                bst=2000
                a=0.95
                pvData <- pvclust(scale(t(rawdata)), method.dist="correlation", method.hclust= n, nboot= bst, parallel=TRUE)

                pdf(paste("bootstrap.STD",s,".CLU",n,".VAR-CORR",cr,".FEA-CORR",cc,".BST",bst,".pdf", sep = ""))
                plot(pvData, hang=-1, cex.pv=.2, cex=.2, float=0)
                pvrect(pvData, alpha=a, pv="au", type="geq", cex.lab=.3)
                dev.off()


                pdf(paste("standardErrors.STD",s,".CLU",n,".VAR-CORR",cr,".FEA-CORR",cc,".BST",bst,".pdf", sep = ""))
                seplot(pvData, type=c("au", "bp"))
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
                pdf(paste("heatmap2.STD",s,".CLU",n,".VAR-CORR",cr,".FEA-CORR",cc,".BST",bst,".pdf", sep = ""))
                heatmap.2(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=palette.green, scale="row", trace="none", RowSideColors = mycolhc, margins=c(5,5), cexRow=.1, cexCol=.1)
                dev.off()

            }
        }
    }
}
