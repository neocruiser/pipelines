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
