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
sampleNames(cel.raw)
#ids <- read.table("summary/sampleIDs.txt")
ids <- read.table("sampleIDs")

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
metadata <- read.table("phenodata", sep = "\t", header = T) %>%
  dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction, ABClikelihood) %>%
  filter(Timepoint != "T2") %>%
#  filter(!PATIENT_ID %in% c(paste0("CNR800",1:4),paste0("CNR900",1:2))) %>% # remove controls
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


# moderated t-statistics (of standard errors) and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
groups = c("relapseCellsFactorial", "multiGrpRelapse", "twoGrpNodes", "multiGrpCells")

g = c("relapseCellsFactorial")
f=1

for (g in groups) {

    if (g == "multiGrpRelapse") {
        design <- model.matrix(~ -1 + metadata$Relapse)
        colnames(design) <- c("noRelapse", "relapse", "control")
        contrast.matrix <- makeContrasts(noRelapse-relapse,
                                         noRelapse-control,
                                         relapse-control,
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast

    } else if (g == "relapseCellsFactorial") {
        sample.factors <- paste(metadata$Prediction, metadata$Relapse, sep=".")
        sample.factors <- factor(sample.factors,
                                 levels = c("ABC.0", "GCB.0", "U.0",
                                            "ABC.1", "GCB.1", "U.1",
                                            "ABC.2", "GCB.2"))
        design <- model.matrix(~0 + sample.factors)
        colnames(design) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(Relapse2ABC=ABC.1-ABC.0,
                                         Relapse2GCB=GCB.1-GCB.0,
                                         ABC2GCB=(ABC.1-ABC.0)-(GCB.1-GCB.0),
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast

    } else if (g == "multiGrpCells") {
        design <- model.matrix(~ -1 + metadata$ABClassify)
        colnames(design) <- c("GCB", "ABC", "Redundant")
        contrast.matrix <- makeContrasts(GCB-ABC,
                                         GCB-Redundant,
                                         ABC-Redundant,
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast

    } else if (g == "twoGrpNodes") {
        design <- model.matrix(~ -1 + metadata$Lymphnodes)
        colnames(design) <- c("Extranodal", "Lymphnodes")
        contrast.matrix <- makeContrasts(Extranodal-Lymphnodes, levels=design)
        coef <- c(1) # refrence to each contrast
    }

    
    # get the description of the two sample groups being compared
    contrast.group <- gsub(" ","",colnames(contrast.matrix)) # diana 3atetne el wa7e

    # get 15% of differentially expressed genes
    selected <- round(dim(trx.normalized)[[1]] * .15)

    
    # Fit Bayesian model and extract Differential Genes (sorted by significance)
    # Benjamini & Hochberg (1995): adjusted "fdr"
    for (f in coef) {
        cel.fit <- lmFit(trx.normalized, design) %>%
            contrasts.fit(contrast.matrix) %>%
            eBayes()
        topTable(cel.fit, coef=f, adjust="fdr", sort.by="B",
                     number=selected) %>%
            write.table(file=paste0(contrast.group[f],".moderated-tstat-bayes.limma.",g,".txt"),
                row.names=FALSE, sep="\t", quote=FALSE) 

        # heatmap and volcanoplots based on p-vals without clustering and bootrstapping
        pdf(paste0(contrast.group[f],".heatmap.tstat-bayes.limma.",g,".pdf"))
        decideTests(cel.fit, p.value=0.000005) %>%
            heatDiagram(cel.fit$coef, primary=f) # is primary the coef??????????
        dev.off()

        pdf(paste0(contrast.group[f],".volvano.tstat-bayes.limma.",g,".pdf"))
        volcanoplot(cel.fit,coef=f,highlight=10)
        dev.off()

        # export significant genes and charts
        # F-test p-values rather than t-test p-values (NESTEDF)
        # Benjamini and Hochberg’s method to control the false discovery rate (GLOBAL)
        index <- c("up", "down")
        pval <- c(0.01, 0.001, 0.0001, "global", "nestedF")
        for (i in index) {
            for (p in pval) {
                pdf(paste0(contrast.group[f],".",i,".P",p,".venn.tstat-bayes.limma.",g,".pdf"))
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


# Colomn names of the Annotated Limma TOptable
#   [1] "transcriptclusterid" "probesetid"          "seqname"            
#   [4] "strand"              "start"               "stop"               
#   [7] "totalprobes"         "geneassignment"      "mrnaassignment"     
#   [10] "swissprot"           "unigene"             "gobiologicalprocess"
#   [13] "gocellularcomponent" "gomolecularfunction" "pathway"            
#   [16] "proteindomains"      "category"            "locustype"          
#   [19] "notes"               "logFC"               "AveExpr"            
#   [22] "t"                   "P.Value"             "adj.P.Val"          
#   [25] "B"



