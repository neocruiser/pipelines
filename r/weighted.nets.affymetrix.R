#################
### RUN CODE ####
#################
pkgs <- c('limma','reshape2','gplots','WGCNA','dplyr','igraph',"RColorBrewer","vegan")
lapply(pkgs, require, character.only = TRUE)

source("./convertMatrix2graph.R")

#LOAD DATA
counts <- as.matrix(read.table("./expressions", header = T, row.names = 1))
tbl_df(counts)

dim(counts)

nco <- dim(counts)[1]

## INCREASING THE POWER AND THRESHOLD REDUCES THE NUMBER OF NETWORKS
if ( nco >= 500 ) {
    pow <- seq(4, 12, 2)
    th <- seq(.1, .6, .1)
} else if ( nco < 500 ) {
    pow <- seq(2, 6, 1)
    th <- seq(.1, .6, .1)
}
th=.5


# sampling for heatmaps
if ( nco >= 100 ){
    heatmap.indices.sampling <- sample(nco, c(nco*.1))
    
} else {
    heatmap.indices.sampling <- seq(1:nco)
}


palette.gr <- brewer.pal(11, name = "PiYG")
palette.rd <- brewer.pal(9, name = "YlOrRd")
palette.green <- colorRampPalette(palette.gr)(n = nco)
palette.red <- colorRampPalette(palette.rd)(n = nco)


### DEBUGGING ###

#require(svMisc)
#for (i in 0:101)
#progress(i, progress.bar = TRUE)

#debug()
#browser()
#traceback()

## output analyses done for debugging purposes
### END ###

    
# standardization
standardize_df <- c("hellinger", "standardize", "range")

## Initialize variable that will contain different networks iterations
networks.summary=NULL
## Initialize variable that will contain differnt module/gene iterations
dm <- NULL
## create data frame to export module associations per gene across iterations
gene2module <- data.frame(ids = rownames(counts))
## create an iterative counter and initialize it
icc <- function(){ i=0; function(){ i <<- i + 1;  i }}
ite <- icc()
steps_done="tmp"
                    

## Create networks
## iterated across different normalization and correlation methods
## each technique generates different modules
for ( s in standardize_df ) {
    ## This normalization step will change the scores
    ## based on the sum by Margins
    ## uses ranges, logs, square roots etc.,
    ## some scores are converted to NaN, this means the method was not compatible
    ## the scores should remain non-negative to successfuly get a graph
    counts <- decostand(x = counts, method = s)

#    counts <- scale(counts, center=T, scale = T)

    allowWGCNAThreads()

    #create similarity matrix
    cordist <- function(dat) {
        cor_matrix  <- cor(t(dat))

        dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))

#        dist_matrix <- as.matrix(dist(scale(dat, center=T, scale=T), diag=TRUE, upper=TRUE))

        # Calculates the log-likelihood
#        dist_matrix <- log1p(dist_matrix)
        # Scale negative values
        dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
#        sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
    }
    sim_matrix <- cordist(counts)

    pdf(paste("similarity.matrix.standardize.FeatureSIZE",nco,".STD",s,".heatmap.pdf",sep = ""))
    heatmap_indices = heatmap.indices.sampling
    heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
              col=palette.green,
              labRow=NA, labCol=NA,
              trace='none', dendrogram='row',
              xlab='Gene', ylab='Gene',
              main='Similarity matrix',
              density.info='none', revC=TRUE)
    dev.off()


    
## Start a progress bar
#pb <- txtProgressBar(min = 0, max = nco, style = 3)
#k=0


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

        pdf(paste("adjacency.matrix.FeatureSIZE",nco,".STD",s,".heatmap.pdf", sep = ""))
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
#        correlate_rows <- c("pearson", "spearman")
        correlate_rows <- c("pearson")        
#        normalize_df <- c("complete", "ward.D2", "average")
        normalize_df <- c("complete")


## redistribute genes into modules based on different correlation strategies
        for ( n in normalize_df ) {
            for ( cr in correlate_rows ) {


                gene_tree <- hclust(as.dist(1-cor(t(adj_matrix),
                                                  method= cr,
                                                  use = "pairwise.complete.obs")), method = n)
                                        #gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")


                ## create only a dendrogram from cluster visualization
                dend <- as.dendrogram(hclust(as.dist(1-cor(t(adj_matrix),
                                                           method= cr,
                                                           use = "pairwise.complete.obs")), method= n))

                ## Get the number of clusters (modules) and the number of genes per cluster
                # max number of genes per module            
                if ( nco >= 100 ) {
                    imax = floor(nco * .1)
                    ival = floor(nco * .01)
                } else if ( nco < 100 ) {
                    imax = nco
                    if ( nco >= 50 ) {
                        ival = floor(nco * .05)           
                    } else {
                        ival = floor(nco * .1)
                    }
                }


                for ( i in seq(5, imax, ival) ) {
                    ## Create a tabulated summary of networks topology
                    module_labels <- cutreeDynamicTree(dendro=gene_tree,
                                                       minModuleSize=i,
                                                       deepSplit=TRUE)

                    # during some iterations on higher power thresholds
                    # the number of modules is 0
                    # so R throws an error, the pipe is then aborted
                    # below is to amend the broken pipe and keep the analysis going
                    nbmods = NULL
                    nbmods = try(summary(module_labels)[[6]], silent = TRUE)
                    if ( is.numeric(nbmods) == TRUE ) {
                        NbModules = nbmods
                    } else {
                        NbModules = NULL
                    }
                    
                    dm <- rbind(dm, data.frame(MaxGenesPerModule = i,
                                               NbModules = NbModules,
                                               Normalization=n, Correlation=cr,
                                               Standardization=s,
                                               SimilaritySize=nco,
                                               EdgeThreshold=th,
                                               CorrelationPower=p))
                }



                ## The mean of the number of modules will be used to cut the dendrogram
                min.mods <- floor(apply(dm[, 1:2], 2, function(x) mean(x)))
                min.mods

                ## number of genes per module
                if ( nco >= 500 ) {
                    mods_a = c((min.mods[[1]] - c( round(min.mods[[1]]/min.mods[[2]]) * 10 ) ) + 1)
                    mods_b = min.mods[[1]]
                    mods_c = c((min.mods[[1]] + c( round(min.mods[[1]]/min.mods[[2]]) * 12 ) ))
                } else if ( nco < 500 ) {
                    mods_a = min(dm[1])
                    mods_b = min.mods[[1]]
                    mods_c = c ( imax - ival )
                }

                # Iterate Clustering Based on the number of genes per module
                for ( fm in c(mods_a, mods_b, mods_c) ) {

                    module_labels <- cutreeDynamicTree(dendro=gene_tree,
                                                       minModuleSize=fm,
                                                       deepSplit=TRUE)

                    pdf(paste("minimum.module.FeatureSIZE",nco,".STD",s,".var-CORR",cr,".CLU",n,".pdf", sep = ""))
                    plot(dm, main = paste("Module (cluster) size selected = ", fm, sep=""))
                    abline(lm(dm$NbModules ~ dm$MaxGenesPerModule), col="red")
                    lines(lowess(dm$MaxGenesPerModule,dm$NbModules), col="blue")
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
                                                                        ".th",t,
                                                                        ".GENperMOD",fm,
                                                                        ".STD",s,
                                                                        ".FeatureSIZE",nco,
                                                                        ".CLU",n,
                                                                        ".varCORR",cr,
                                                                        ".graphml",
                                                                        sep = "" ),
                                                       threshold=t,
                                                       nodeAttrDataFrame=df,
                                                       verbose = TRUE)

                    }


                    ## Count all genes and their connections per method used, threshold, power.
                    ## Group all nodes and edges
                    adj_matrix[abs(adj_matrix) < th] <- 0
                    orphaned <- (colSums(adj_matrix) == 0)
                    fit_matrix <- adj_matrix[!orphaned, !orphaned]

                    fx <- graph.adjacency(fit_matrix, mode='undirected', weighted=TRUE, diag=FALSE)

                    degree(fx)[1:10]

                    x=degree(fx)
                    for ( l in seq(1, nco, c(floor(nco * .01) + 1) )) {
                        a <- length(x[ x > l ])
                        networks.summary <- rbind(networks.summary,
                                                  data.frame(MaxEdgesPerGene= l, NbNodes= a,
                                                             Normalization=n, Correlation=cr,
                                                             Standardization=s,
                                                             MaxGenesPerModule=fm,
                                                             SimilaritySize=nco,
                                                             EdgeThreshold=th,
                                                             CorrelationPower=p))
                    }

                    ## output analyses done for debugging purposes
                    if ( file.exists(steps_done) )
                        file.remove(steps_done)
                    
                    ## increase counter by 1 and amend new methods succesfully executed
                    ii <- ite()
                    steps_done=paste0("network.iteration_",ii,".po",p,
                                      ".th",t,".GENperMOD",fm,".STD",s,".FeatureSIZE",nco,
                                      ".CLU",n,".varCORR",cr,".tmp")

                    if ( !file.exists(steps_done) )
                        file.create(steps_done)

                    ## add modules, one iteration at a time
                    gene2module <- cbind(gene2module, module_labels)
                    names(gene2module)[c( 1 + ii)] <- paste0("GMOD",fm,"PO",p,"TH",t,"_STD",s,"CLU",n,"COR",cr)
                
                }
            }
        }
    }
}


## Create a summarized tabulated file about clustering, network organization, and correlation
write.table(gene2module, "ids2modules.summary.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(dm, "modules.summary.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(networks.summary, "networks.summary.txt", quote=FALSE, sep="\t", row.names=FALSE)



#save(file = "log.Rdata")
disableWGCNAThreads()
gc()

sink("r_session.info")
sessionInfo()
sink()
