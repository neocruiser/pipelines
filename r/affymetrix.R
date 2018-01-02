# Load packages
pkgs <- c('RColorBrewer', 'genefilter', 'limma',
          'pvclust', 'foreach', 'oligo', 'pd.hta.2.0',
          'dplyr', 'plyr', 'reshape', 'tidyr', 'doMC')
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

# Microarray files loaded into array
# Robust Multi-Chip average for background correction, normalization
# expressions are log2 transformed
cel.raw <- list.celfiles(full=TRUE, listGzipped=FALSE) %>%
    read.celfiles()
sampleNames(cel.raw)
#ids <- read.table("summary/sampleIDs.txt")
ids <- read.table("sampleIDs")

# Choose charts colors
mypar <- function(row=1, col=1)
    par(mar=c(2.5, 2.5, 1.6, 1.1),
        mgp=c(1.5, 0.5, 0),
        mfrow=c(row, col))
palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)

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

# get probeset info
pInfo <- getProbeInfo(cel.raw, target="probeset",field=c("fid","type"),sortBy="none")
dim(cel.raw); dim(pInfo)
sink("probe.info.affymetrix.chip.txt")
table(pInfo$type[pInfo$fid %in% rownames(cel.raw)])
sink()


# Probe level model fitted to the raw data
# normalized unscaled standard errors and relative log expression
plmFit <- fitProbeLevelModel(cel.raw, target='core')
cols <- rep(darkColors(nlevels(cel.raw$Prediction)), each=2)
mypar(2, 1)
pdf("nuse.plm.raw.pdf")
NUSE(plmFit, col=cols)
RLE(plmFit, col=cols)
dev.off();dev.off()

## Samples classification and experimental designs
#metadata <- read.table("summary/phenodata.txt", sep = "\t", header = T) %>%
metadata <- read.table("phenodata.txt", sep = "\t", header = T) %>%
  dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction, ABClikelihood) %>%
  filter(Timepoint != "T2") %>%
  mutate(Relapse = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                         "CNS_RELAPSE_CHOPorEQUIVALENT",
                                         "CNS_DIAGNOSIS") ~ 1,
                             GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ 0,
                             GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ 0,
                             TRUE ~ 2)) %>%
  mutate(ABClassify = case_when(ABClikelihood >= .9 ~ 1,
                                ABClikelihood <= .1 ~ 0,
                                TRUE ~ 2)) %>%
  mutate(Lymphnodes = case_when(SITE == "LN" ~ 1, TRUE ~ 0))


# make sure all samples preserve their ID
metadata$Relapse <- as.factor(metadata$Relapse)
metadata$ABClassify <- as.factor(metadata$ABClassify)
metadata$Lymphnodes <- as.factor(metadata$Lymphnodes)
metadata <- metadata[metadata$SAMPLE_ID %in% ids$V1, ]
row.names(metadata) = metadata$SAMPLE_ID
colnames(cel.raw) = metadata$SAMPLE_ID
pd <- AnnotatedDataFrame(data=metadata)
sampleNames(pd) <- metadata$SAMPLE_ID
phenoData(cel.raw) <- pd
sink("metadata.samples.info.txt")
summary(metadata)
sink()

# RMA normalization
trx.normalized <- oligo::rma(cel.raw, target='core')
probe.normalized <- oligo::rma(cel.raw, target='probeset')
write.exprs(trx.normalized, file="normalized.trx.expression.txt")
write.exprs(probe.normalized, file="normalized.probe.expression.txt")


# Pull affymetrix annotations for genes and exons
featureData(trx.normalized) <- getNetAffx(trx.normalized, 'transcript')
sink("annotation.gene.exon.affy.txt")
with(fData(trx.normalized), table(seqname, category))
sink()


# moderated t-statistics and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
groups = c("relapse", "cells", "nodes")
for (g in groups) {
    if (g == "relapse") {
        design <- model.matrix(~ -1 + metadata$Relapse)
        colnames(design) <- c("noRelapse", "relapse", "control")
        contrast.matrix <- makeContrasts(noRelapse-relapse,
                                         noRelapse-control,
                                         relapse-control,
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast
    } else if (g == "cells") {
        design <- model.matrix(~ -1 + metadata$ABClassify)
        colnames(design) <- c("GCB", "ABC", "control")
        contrast.matrix <- makeContrasts(GCB-ABC,
                                         GCB-control,
                                         ABC-control,
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast
    } else if (g == "nodes") {
        design <- model.matrix(~ -1 + metadata$Lymphnodes)
        colnames(design) <- c("Extranodal", "Lymphnodes")
        contrast.matrix <- makeContrasts(Extranodal-Lymphnodes,
                                         levels=design)
        coef <- c(1) # refrence to each contrast
    }

    # get the description of the two sample groups being compared
    contrast.group <- gsub(" ","",colnames(contrast.matrix)) # diana 3atetne el wa7e

    # Fit Bayesian model and extract Differential Genes (sorted by significance)
    for (f in coef) {
        cel.fit <- lmFit(trx.normalized, design) %>%
            contrasts.fit(contrast.matrix) %>%
            eBayes()
        topTable(cel.fit, coef=f, adjust="fdr", sort.by="B",
                     number=dim(trx.normalized)[[1]]) %>%
            write.table(file=paste0(contrast.group[f],".moderated.tstat.bayes.limma.txt"),
                row.names=TRUE, sep="\t", quote=FALSE) 

# summarizing differential expressions
        print("Number of genes significant (adjP < 0.05) in this list:")
        topTable(cel.fit, coef=f, adjust="fdr",
                 sort.by="P",
                 number=dim(trx.normalized)[[1]]) %>%
            filter(adj.P.Val < 0.05) %>%
            dim

        print("Number of genes significant (adjP < 0.05, folds over 2, avg expression over 10) in this list:")
        topTable(cel.fit, coef=f, adjust="fdr", sort.by="P", number=dim(trx.normalized)[[1]]) %>%
            filter(adj.P.Val < 0.01) %>%
            filter(logFC > 1 | logFC < -1) %>%
            filter(AveExpr > 10) %>%
            dim

# heatmap and volcanoplots based on p-vals without clustering and bootrstapping
        pdf(paste0(contrast.group[f],".heatmap.tstat.bayes.limma.pdf"))
        decideTests(cel.fit, p.value=0.000005) %>%
            heatDiagram(cel.fit$coef, primary=f) # is primary the coef??????????
        dev.off()

        pdf(paste0(contrast.group[f],".volvano.tstat.bayes.limma.pdf"))
        volcanoplot(cel.fit,coef=f,highlight=10)
        dev.off()




# export significant genes and charts
# F-test p-values rather than t-test p-values (NESTEDF)
# Benjamini and Hochberg’s method to control the false discovery rate (GLOBAL)
        index <- c("up", "down")
        pval <- c(0.01, 0.001, 0.0001, "global", "nestedF")
        for (i in index) {
            for (p in pval) {
                pdf(paste0(contrast.group[f],".",i,".P",p,".venn.tstat.bayes.limma.pdf"))
                if (p>0) {
                    decideTests(cel.fit, p.value=p) %>%
                        vennDiagram(include=i)
                } else {
                    decideTests(cel.fit, method=p) %>%
                        vennDiagram(include=i)
                }
                dev.off()
            }
        }

    }

}
