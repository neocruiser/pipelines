pkgs <- c('limma','reshape2','gplots','WGCNA','dplyr','igraph')
lapply(pkgs, require, character.only = TRUE)

pow <- 10
th <- 0.5
#load data
#counts <- read.table("../../R/ganglia/data/diffExpr.P1e-4_C2.matrix.log2.dat", header = T)
counts <- read.table("./diffExpr.P1e-4_C2.matrix.log2.dat", header = T)
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
#Convert similarity matrix to adjacency matrix.
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=pow, type='signed')
rm(sim_matrix)
gc()
## gene ids are Trinity IDs
gene_ids <- rownames(adj_matrix)
adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids


## Detect co-expression modules
## Hierarchical clustering first
gene_tree <- hclust(as.dist(1-cor(t(adj_matrix),
                                  method="pearson")),
                    method="average")
#gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")


## create only a dendrogram from cluster visualization
dend <- as.dendrogram(hclust(as.dist(1-cor(t(adj_matrix),
                                           method="pearson")),
                             method="average"))
#x11(); plot(dend);gc()

## Get the number of clusters (modules) and the number of genes per cluster
d <- NULL
imax=200
for ( i in seq(15,imax,10) ) {
    module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=i,
                                       deepSplit=TRUE)
    d <- rbind(d, data.frame(genes = i, modules = summary(module_labels)[[6]]))
}
## The mean of the number of clusters will be used to cut the dendrogram
min.mods <- apply(d, 2, function(x) mean(x))
fm <- floor(min.mods[[1]])
fm <- floor(((imax-fm)/1.5) + fm)
fm
module_labels <- cutreeDynamicTree(dendro=gene_tree,
                                   minModuleSize=fm,
                                   deepSplit=TRUE)
pdf("minimum.module.sizes.pdf")
plot(d, main = paste("Module (cluster) size selected = ", fm, sep=""))
abline(lm(d$modules ~ d$genes), col="red")
lines(lowess(d$genes,d$modules), col="blue")
dev.off()

module_colors <- labels2colors(module_labels)
gene_info <- data.frame(id = gene_ids, modules=module_colors)
gene_info$color_rgb<- col2hex(gene_info$modules)


#cat contigs.deseq2.p4.c2.prot.fa.tsv | sed 's/ /./g' | cut -f1,9,13 | awk '{if($2<=0.00001)print$1,$3}' | sed 's/_. / /g' | sort - | uniq -u | wc -l


### Merge annotated contigs with coexpressed modules
dim(gene_info)
annotations <- read.table("./contigs.deseq2.p4.c2.tsv_id2description.txt", fill = TRUE, na.strings = c("", "NA")) %>% na.omit()
tbl_df(annotations)
annotations <- annotations[!duplicated(annotations[,1]), ]

df <- merge(gene_info, annotations, by.x = "id", by.y = "V1", all.x = T)
df <- df[!duplicated(df$id),]
dim(df)
colnames(df) <- c('gene_id','modules','colors_rgb','description')
df$description <- as.character(df$description)
head(df)


# Export network into file
export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {

    # Determine filename to use
    if (is.null(filename)) {
        filename='network.graphml'
    }

    # TODO 2015/04/09
    # Add option to rescale correlations for each module before applying
    # threshold (this is simpler than the previous approach of trying to
    # determine a different threshold for each module)
    #
    # Still, modules with very low correlations should be given somewhat
    # less priority than those with very high correlations.

    #module_colors <- unique(nodeAttrDataFrame$color)
    #module_genes <- which(nodeAttrDataFrame$color == color)
    #module_adjmat <- adj_mat[module_genes,]
    #num_genes <- length(module_genes)

    # Adjust threshold if needed to limit remaining edges
    max_edges <- max_edge_ratio * nrow(adj_mat)

    edge_to_total_ratio <- max_edges / length(adj_mat)
    edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))

    # Also choose a minimum threshold to make sure that at least some edges
    # are left
    min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))

    threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))

    # Remove edges with weights lower than the cutoff
    adj_mat[abs(adj_mat) < threshold] <- 0

    # Drop any genes with no edges (TODO: Make optional)
    orphaned <- (colSums(adj_mat) == 0)
    adj_mat <- adj_mat[!orphaned, !orphaned]

    # Also remove annotation entries
    if (!is.null(nodeAttr)) {
        nodeAttr <- nodeAttr[!orphaned]
    }
    if (!is.null(nodeAttrDataFrame)) {
        nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
    }

    # Keep track of non-positive edges and rescale to range 0,1
    is_zero     <- adj_mat == 0
    is_negative <- adj_mat < 0

    adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
    adj_mat[is_zero] <- 0
    adj_mat[is_negative] <- -adj_mat[is_negative]

    if (verbose) {
        message(sprintf("Outputting matrix with %d nodes and %d edges",
                        nrow(adj_mat), sum(adj_mat > 0)))
    }

    # Create a new graph and add vertices
    # Weighted graph
    if (weighted) {
        g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
    } else {
        adj_mat[adj_mat != 0] <- 1
        g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
    }

    # Add single node annotation from vector
    if (!is.null(nodeAttr)) {
        g <- set.vertex.attribute(g, "attr", value=nodeAttr)
    }

    # Add node one or more node annotations from a data frame
    if (!is.null(nodeAttrDataFrame)) {
        for (colname in colnames(nodeAttrDataFrame)) {
            g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
        }
    }

    edge_correlation_negative <- c()

    # neg_correlations[edge_list]
    edge_list <- get.edgelist(g)

    for (i in 1:nrow(edge_list)) {
        from <- edge_list[i, 1]
        to   <- edge_list[i, 2]
    }

    # Save graph to a file
    write.graph(g, filename, format='graphml')

    # return igraph
    return(g)
}
#extract network
g <- export_network_to_graphml(adj_matrix, filename='./network.graphml',
                               threshold=th, nodeAttrDataFrame=df)

disableWGCNAThreads()
gc()
