# Load packages
pkgs <- c('RColorBrewer', 'genefilter', 'gcrma', 'limma', 'simpleaffy',
          'affyQCReport', 'plier', "affycoretools", 'pdInfoBuilder',
          'pvclust', 'foreach', 'affy', 'affyPLM', 'oligo', 'pd.hta.2.0',
          'dplyr', 'plyr', 'reshape', 'tidyr')
lapply(pkgs, require, character.only = TRUE)

# Multicore processing
#mp = 8
#Sys.setenv(R_THREADS = mp)

# Parallel processing
ldPath()
ocSamples(100)
ocProbesets(200)
require(ff)
registerDoMC(2)

# Choose charts colors
palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)

# Microarray files loaded into array
# Robust Multi-Chip average for background correction, normalization
# expressions are log2 transformed
cel.normalized <- list.celfiles() %>%
    read.celfiles() %>%
    oligo::rma()
write.exprs(cel.normalized,file="data.txt")

# moderated t-statistics and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3)))
colnames(design) <- c("group1", "group2", "group3")
contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
fit2 <- lmFit(cel.normalized, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()
top.genes <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=dim(cel.normalized)[[1]]) 
write.table(top.genes, file="moderated.tstat.bayes.limma.txt", row.names=F, sep="\t") 

# export significant genes and charts
pdf("venn.tstat.bayes.limma.pdf")
decideTests(fit2, p.value=0.05) %>%
    vennDiagram()
dev.off()

print("Number of genes significant (adjP < 0.05) in this list:")
topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=dim(cel.normalized)[[1]]) %>%
    filter(adj.P.Val < 0.05) %>%
    dim

print("Number of genes significant (adjP < 0.05, folds over 2, avg expression over 10) in this list:")
topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=dim(cel.normalized)[[1]]) %>%
    filter(adj.P.Val < 0.01) %>%
    filter(logFC > 1 | logFC < -1) %>%
    filter(AveExpr > 10) %>%
    dim

pdf("heatmap.tstat.bayes.limma.pdf")
decideTests(fit2, p.value=0.000005) %>%
    heatDiagram(fit2$coef, primary=1)
dev.off()
