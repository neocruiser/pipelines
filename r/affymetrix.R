# Load packages
pkgs <- c('RColorBrewer', 'genefilter', 'limma',
          'pvclust', 'foreach', 'oligo', 'pd.hta.2.0',
          'dplyr', 'plyr', 'reshape', 'tidyr', 'doMC')
#pkgs2 <- c('gcrma', 'simpleaffy', 'affyQCReport', 'plier', "affycoretools", 'affy', 'affyPLM')
lapply(pkgs, require, character.only = TRUE)

# Multicore processing
#mp = 8
#Sys.setenv(R_THREADS = mp)

# Parallel processing
#registerDoMC(12)
#ldPath()
#ocSamples(50)
#ocProbesets(200)
#dir.create('ffObjs')
#ldPath('ffObjs')

# Choose charts colors
mypar <- function(row=1, col=1)
    par(mar=c(2.5, 2.5, 1.6, 1.1),
        mgp=c(1.5, 0.5, 0),
        mfrow=c(row, col))
palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)

# Microarray files loaded into array
# Robust Multi-Chip average for background correction, normalization
# expressions are log2 transformed
cel.raw <- list.celfiles(full=TRUE, listGzipped=FALSE) %>%
    read.celfiles()

## Samples classification and experimental designs
sampleNames(cel.raw)
ids <- read.table("sampleIDs.txt")
metadata <- read.table("phenodata.txt", sep = "\t", header = T) %>%
  dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction, IPI_GROUP) %>%
  mutate(GROUP_NEW = ifelse(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                         "CNS_RELAPSE_CHOPorEQUIVALENT",
                                         "CNS_DIAGNOSIS"), "CNS",
                     ifelse(GROUP == "TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE", GROUP))) %>%
  filter(!SAMPLE_ID %in% c("CNR6038T1", "CNR7010T1")) %>%
  as.data.frame() %>%
  filter(!GROUP %in% c("NORMAL_GCB_CONTROL", "NORMAL_ABC_CONTROL")) %>%
  filter(Timepoint != "T2") %>%
  mutate(SITE_NEW = ifelse(SITE == "LN", "LN", "EN"))
pd <- new('AnnotatedDataFrame', data=metadata)
sampleNames(pd) <- metadata$SAMPLE_ID
phenoData(cel.raw) <- pd

# RMA normalization
trx.normalized <- oligo::rma(cel.raw, target='core')
probe.normalized <- oligo::rma(cel.raw, target='probeset')
write.exprs(trx.normalized, file="normalized.trx.expression.txt")
write.exprs(probe.normalized, file="normalized.probe.expression.txt")

# get probeset info
pInfo <- getProbeInfo(cel.raw, target="probeset",field=c("fid","type"),sortBy="none")
dim(cel.raw); dim(pInfo)
table(pInfo$type[pInfo$fid %in% rownames(cel.raw)])

# Log intensities
pdf("boxplot.raw.pdf")
boxplot(cel.raw, target="core")
dev.off()
pdf("histogram.raw.pdf")
hist(cel.raw, target="core")
dev.off()
pdf("ma.raw.pdf")
MAplot(cel.raw[, 1:4], pairs=TRUE)
dev.off()

smallData <- rawData[, 1:12]
grps <- as.character(smallData$tissue)
grps <- as.factor(grps)
MAplot(smallData, pairs=TRUE, groups=grps)


plmFit <- fitProbeLevelModel(cel.raw, target='core')
#cols <- rep(darkColors(nlevels(rawData$tissue)), each=3)
mypar(2, 1)
pdf("nuse.plm.raw.pdf")
NUSE(plmFit, col=palette.green)
RLE(plmFit, col=palette.red)
dev.off()






# moderated t-statistics and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3)))
colnames(design) <- c("group1", "group2", "group3")
contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)


design <- model.matrix(~ -1+factor(c(1,1,2,2)))
colnames(design) <- c("group1", "group2")
contrast.matrix <- makeContrasts(group2-group1, levels=design)
fit2 <- lmFit(trx.normalized, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()
top.genes <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=dim(trx.normalized)[[1]]) 
write.table(top.genes, file="moderated.tstat.bayes.limma.txt", row.names=TRUE, sep="\t", quote=FALSE) 

# export significant genes and charts
# F-test p-values rather than t-test p-values (NESTEDF)
# Benjamini and Hochberg’s method to control the false discovery rate (GLOBAL)
index <- c("up", "down")
pval <- c(0.001, 0.0001, "global", "nestedF")
for (i in index) {
    for (p in pval) {
        pdf(paste(i, ".", p, ".venn.tstat.bayes.limma.pdf", sep=""))
        if (p>0) {
            decideTests(fit2, p.value=p) %>%
                vennDiagram(include=i)
        } else {
            decideTests(fit2, method=p) %>%
                vennDiagram(include=i)
        }
        dev.off()
    }
}


print("Number of genes significant (adjP < 0.05) in this list:")
topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=dim(trx.normalized)[[1]]) %>%
    filter(adj.P.Val < 0.05) %>%
    dim

print("Number of genes significant (adjP < 0.05, folds over 2, avg expression over 10) in this list:")
topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=dim(trx.normalized)[[1]]) %>%
    filter(adj.P.Val < 0.01) %>%
    filter(logFC > 1 | logFC < -1) %>%
    filter(AveExpr > 10) %>%
    dim

pdf("heatmap.tstat.bayes.limma.pdf")
decideTests(fit2, p.value=0.000005) %>%
    heatDiagram(fit2$coef, primary=1)
dev.off()
