## INCREASING THE POWER AND THRESHOLD REDUCES THE NUMBER OF NETWORKS
pow <- seq(2, 7, 1)
th <- seq(.4, .5, .1)

#################
### RUN CODE ####
#################
pkgs <- c('limma','reshape2','gplots','WGCNA','dplyr','igraph',"RColorBrewer","vegan")
lapply(pkgs, require, character.only = TRUE)

palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)

source("./convertMatrix2graph.R")

#LOAD DATA
counts <- read.table("./logs", header = T, row.names = 1)
counts <- as.matrix(counts[, -1])
tbl_df(counts)



# standardization
standardize_df <- c("standardize", "range", "log") # hellinger
for ( s in standardize_df ) {
    counts <- decostand(x = counts, method = s)
    nco <- dim(counts)[1]

    allowWGCNAThreads()
    #create similarity matrix
    cordist <- function(dat) {
        cor_matrix  <- cor(t(dat))

        dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
        dist_matrix <- log1p(dist_matrix)
        dist_matrix <- 1 - (dist_matrix / max(dist_matrix))

        sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
    }
    sim_matrix <- cordist(counts)

    pdf(paste("similarity.matrix.SSIZE",nco,".STD",s,".heatmap.pdf",sep = ""))
    heatmap_indices <- sample(nrow(sim_matrix), c(nco*.1))
    heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
              col=palette.green,
              labRow=NA, labCol=NA,
              trace='none', dendrogram='row',
              xlab='Gene', ylab='Gene',
              main='Similarity matrix',
              density.info='none', revC=TRUE)
    dev.off()


## construct different networks based on their power analysis
    for(p in pow) {

        #Convert similarity matrix to adjacency matrix.
        adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=p, type='signed')
        gc()
        ## gene ids are serial numbers
        gene_ids <- rownames(adj_matrix)
        adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
        rownames(adj_matrix) <- gene_ids
        colnames(adj_matrix) <- gene_ids

        pdf(paste("adjacency.matrix.SSIZE",nco,".STD",s,".heatmap.pdf", sep = ""))
        heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
                  col=palette.green,
                  labRow = NA, labCol=NA,
                  trace='none', dendrogram='row',
                  xlab='Gene', ylab='Gene',
                  main='Adjacency matrix',
                  density.info='none', revC=TRUE)
        dev.off()

        ## Detect co-expression modules
        ## Hierarchical clustering first
        correlate_rows <- c("pearson", "spearman")
        normalize_df <- c("complete", "ward.D2", "average")


## redistribute genes into modules based on different correlation strategies
        for ( n in normalize_df ) {
            for ( cr in correlate_rows ) {


                gene_tree <- hclust(as.dist(1-cor(t(adj_matrix),
                                                  method= cr)), method = n)
                                        #gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")


                ## create only a dendrogram from cluster visualization
                dend <- as.dendrogram(hclust(as.dist(1-cor(t(adj_matrix),
                                                           method= cr)), method= n))

                ## Get the number of clusters (modules) and the number of genes per cluster
                d <- NULL

                # max number of genes per module            
                imax = c(nco * .1)
                ival = c(nco * .001)

                for ( i in seq(5, imax, ival) ) {
                    module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=i,
                                                       deepSplit=TRUE)
                    d <- rbind(d, data.frame(genes = i, modules = summary(module_labels)[[6]]))
                }

                ## The mean of the number of modules will be used to cut the dendrogram
                min.mods <- floor(apply(d, 2, function(x) mean(x)))
                min.mods

                ## number of genes per module
                mods_a = min.mods[[1]] - 50
                mods_b = min.mods[[1]]
                mods_c = min.mods[[1]] + 50
                
                # Iterate clustering based on the number of genes per module
                for ( fm in c(mods_a, mods_b, mods_c) ) {
                #    for ( f in c(1, 2) ) {
                #fm <- floor(min.mods[[f]])
                #fm <- floor(((imax-fm)/2.5) + fm)

                    module_labels <- cutreeDynamicTree(dendro=gene_tree,
                                                       minModuleSize=fm,
                                                       deepSplit=TRUE)
                    pdf(paste("minimum.module.SSIZE",nco,".STD",s,".var-CORR",cr,".CLU",n,".pdf", sep = ""))
                    plot(d, main = paste("Module (cluster) size selected = ", fm, sep=""))
                    abline(lm(d$modules ~ d$genes), col="red")
                    lines(lowess(d$genes,d$modules), col="blue")
                    dev.off()

                    module_colors <- labels2colors(module_labels)
                    gene_info <- data.frame(id = gene_ids, modules=module_colors)
                    gene_info$color_rgb<- col2hex(gene_info$modules)

                    ### Merge annotated contigs with coexpressed modules
                    tbl_df(gene_info)
                    dim(adj_matrix)
                    df=gene_info

                    #EXTRACT NETWORK
                    for (t in th){
                        g <- export_network_to_graphml(adj_matrix,
                                                       filename = paste("network.POW",p,
                                                                        ".Th",t,
                                                                        ".GEN",fm,
                                                                        ".STD",s,
                                                                        ".SSIZE",nco,
                                                                        ".CLU",n,
                                                                        ".var-CORR",cr,
                                                                        ".graphml",
                                                                        sep = "" ),
                                                       threshold=t,
                                                       nodeAttrDataFrame=df)
                    }
                }
            }
            
        }

        
    }
}


#save(file = "log.Rdata")
disableWGCNAThreads()
gc()

sink("r_session.info")
sessionInfo()
sink()
