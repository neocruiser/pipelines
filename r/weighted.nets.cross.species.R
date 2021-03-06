## INCREASING THE POWER AND THRESHOLD REDUCES THE NUMBER OF NETWORKS
pow <- seq(1, 4, 1)
th <- seq(.2, .5, .1)

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

#load data
# a = species 1
# b = species 2
counts.a <- read.table("./logs.a", header = T)
counts.b <- read.table("./logs.b", header = T)
counts.a <- counts.a[, c(1:3,6,5,4,8,9,7)]
counts.b <- counts.b[, c(3,2,1,4,5,6,8,7,9)]
counts <- rbind(counts.a, counts.b)

standardize_df <- c("standardize", "range", "log", "hellinger")
for(std in standardize_df){
counts <- decostand(x = counts, method = std)

counts <- data.frame(contigs = rownames(counts), counts) %>%
    arrange(desc(contigs))
rownames(counts) <- counts$contigs
counts <- counts[, -1]
tbl_df(counts)

number.genes <- dim(counts)[1]

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

pdf(paste("./pdf/similarity.matrix.SSIZE",number.genes,".STD",std,".heatmap.pdf",sep = ""))
heatmap_indices <- sample(nrow(sim_matrix), number.genes/2)
heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
            col=palette.red,
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

    pdf(paste("./pdf/adjacency.matrix.SSIZE",number.genes,".STD",std,".heatmap.pdf", sep = ""))
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
    correlate_rows <- c("pearson","spearman")
    normalize_df <- c("complete","ward.D2","average")
    for(n in normalize_df){
        for(cr in correlate_rows){

    gene_tree <- hclust(as.dist(1-cor(t(adj_matrix),
                                  method=cr)),
                    method=n)


## create only a dendrogram from cluster visualization
dend <- as.dendrogram(hclust(as.dist(1-cor(t(adj_matrix),
                                           method=cr)),
                             method=n))
#x11(); plot(dend);gc()

## Get the number of clusters (modules) and the number of genes per cluster
d <- NULL
imax=20
for ( i in seq(5,imax,5) ) {
    module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=i,
                                       deepSplit=TRUE)
    d <- rbind(d, data.frame(genes = i, modules = summary(module_labels)[[6]]))
}
## The mean of the number of clusters will be used to cut the dendrogram
min.mods <- apply(d, 2, function(x) mean(x))
# change the number of genes per cluster
    for ( fm in c(5, 10, 20) ) {
#    for ( f in c(1, 2) ) {
#fm <- floor(min.mods[[f]])
#fm <- floor(((imax-fm)/2.5) + fm)
fm
module_labels <- cutreeDynamicTree(dendro=gene_tree,
                                   minModuleSize=fm,
                                   deepSplit=TRUE)
pdf(paste("./pdf/minimum.module.SSIZE",number.genes,".STD",std,".var-CORR",cr,".CLU",n,".pdf", sep = ""))
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
        annotations <- read.table("./qpx.clam.ids", fill = TRUE, na.strings = c("", "NA"))
#        annotations <- read.table("./contigs.deseq2.p4.c2.prot.fa.tsv.id2description.NR-PTHR-IPS.diamond5.LEN20.EVAL5.txt", fill = TRUE, na.strings = c("", "NA"))
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
                                       filename = paste("network.POW",p,
                                                        ".Th",t,
                                                        ".GEN",fm,
                                                        ".STD",std,
                                                        ".SSIZE",number.genes,
                                                        ".CLU",n,
                                                        ".varCORR",cr,
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
