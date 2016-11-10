## INCREASING THE POWER AND THRESHOLD REDUCES THE NUMBER OF NETWORKS
pow <- seq(6, 8, 1)
th <- seq(.4, .5, .1)

#################
### RUN CODE ####
#################
pkgs <- c('limma','reshape2','gplots','WGCNA','dplyr','igraph')
lapply(pkgs, require, character.only = TRUE)

source("./convertMatrix2graph.R")

#load data
#counts <- read.table("../../R/ganglia/data/diffExpr.P1e-4_C2.matrix.log2.dat", header = T)
counts <- read.table("./diffExpr.P1e-4_C2.matrix.log2.dat", header = T)
counts <- data.frame(contigs = rownames(counts), counts) %>%
    arrange(desc(contigs))
rownames(counts) <- counts$contigs
counts <- counts[, -1]
tbl_df(counts)

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

pdf("similarity.matrix.sample.heatmap.pdf")
heatmap_indices <- sample(nrow(sim_matrix), 500)
heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
            col=redgreen(75),
            labRow=NA, labCol=NA, 
            trace='none', dendrogram='row',
            xlab='Gene', ylab='Gene',
            main='Similarity matrix',
            density.info='none', revC=TRUE)
dev.off()

for(p in pow) {

#Convert similarity matrix to adjacency matrix.
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=p, type='signed')
gc()
## gene ids are Trinity IDs
gene_ids <- rownames(adj_matrix)
adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids

    pdf("adjacency.matrix.heatmap.pdf")
    heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
              col=redgreen(75),
              labRow = NA, labCol=NA, 
              trace='none', dendrogram='row',
              xlab='Gene', ylab='Gene',
              main='Adjacency matrix',
              density.info='none', revC=TRUE)
    dev.off()

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
# change the number of genes per cluster
    for ( fm in c(25, 50, 100) ) {
#    for ( f in c(1, 2) ) {        
#fm <- floor(min.mods[[f]])
#fm <- floor(((imax-fm)/2.5) + fm)
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


### Merge annotated contigs with coexpressed modules
tbl_df(gene_info)
dim(adj_matrix)
# this simply removes annotated genes without a description content being found in any gene database.
# however there might be another of the same annotated gene with a description. this gene is a duplicate and will remain in the data frame
### To make the annotation file, merge IPS output and Panther output
annotations <- read.table("./contigs.deseq2.p4.c2.prot.fa.tsv.id2description.NR-PTHR-IPS.diamond5.LEN20.EVAL5.txt", fill = TRUE, na.strings = c("", "NA"))
tbl_df(annotations)
        
df <- merge(gene_info, annotations, by.x = "id", by.y = "V1", all.x = T)
df <- df[!duplicated(df$id),]
colnames(df) <- c('gene_id','modules','colors_rgb','description')
df$description <- as.character(df$description)

df <- arrange(df, desc(gene_id))

        tbl_df(df)
        
#extract network
    for (t in th){
        g <- export_network_to_graphml(adj_matrix,
                                       filename = paste("network.PVAL4.FOLD2.POW",p,".Th",t,".GEN",fm,".graphml",sep = "" ),
                                       threshold=t, nodeAttrDataFrame=df)
    }
    }
}

save(file = "log.Rdata")
disableWGCNAThreads()
gc()
