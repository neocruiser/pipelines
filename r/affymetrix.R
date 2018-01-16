# Load packages
pkgs <- c('RColorBrewer', 'genefilter', 'limma',
          'pvclust', 'foreach', 'oligo', 'pd.hta.2.0',
          'plyr', 'dplyr', 'reshape', 'tidyr', 'doMC')
lapply(pkgs, require, character.only = TRUE)

# Multicore processing
#mp = 8
#Sys.setenv(R_THREADS = mp)

# Parallel processing
#registerDoMC(12)
#require(ff)
#ldPath()
#ocSamples(200)
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


#########################
####  Read samples   ####
#########################
# Microarray files loaded into array
cel.raw <- list.celfiles("../raw", full=TRUE, listGzipped=FALSE) %>%
    read.celfiles()
length(sampleNames(cel.raw))
ids <- read.table("summary/sampleIDs")
gc()


###########################
####  Quality control  ####
###########################
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
gc()

# get probeset info
pInfo <- getProbeInfo(cel.raw, target="probeset",field=c("fid","type"),sortBy="none")
dim(cel.raw); dim(pInfo)
sink("probe.info.affymetrix.chip.txt")
table(pInfo$type[pInfo$fid %in% rownames(cel.raw)])
sink()
gc()

# Probe level model fitted to the raw data
# scaling the standard errors and relative log expression
#plmFit <- fitProbeLevelModel(cel.raw, target='core')
#cols <- rep(darkColors(nlevels(cel.raw$Prediction)), each=2)
#mypar(2, 1)
#pdf("nuse.plm.raw.pdf")
#NUSE(plmFit, col=cols)
#RLE(plmFit, col=cols)
#dev.off();dev.off()
#gc()


########################################
####  Differential gene expression  ####
########################################
## Wrap of essential functions for fitting gene expression
## standard error standardization
## output plots: volcano, heatmaps, venn
## output summary: log transformed expressions
moderatedFit <- function(data=trx.normalized, contrasts=contrast.matrix, labels="unsetGroup", pval=0.001, coef=coef, percent=.15){
# Fit Bayesian model and extract Differential Genes (sorted by significance)
# Benjamini & Hochberg (1995): adjusted "fdr"
# get the description of the two sample groups being compared
    contrast.group <- gsub(" ","",colnames(contrasts)) # diana 3atetne el wa7e
    
# get a proportion of differentially expressed genes
    selected <- round(dim(data)[[1]] * percent)

    
    for (f in coef) {

        cel.fit <- lmFit(data, design) %>%
            contrasts.fit(contrasts) %>%
            eBayes()
        
        topTable(cel.fit, coef=f, adjust="fdr", sort.by="B",
                 number=selected) %>%
            write.table(file=paste0(contrast.group[f],".",labels,
                                    ".moderated-tstat-bayes.limma.txt"),
                        row.names=FALSE, sep="\t", quote=FALSE) 

# heatmap and volcanoplots based on p-vals without clustering and bootrstapping
        pdf(paste0(contrast.group[f],".",labels,".heatmap.tstat-bayes.limma.pdf"))
        decideTests(cel.fit, p.value=pval) %>%
            heatDiagram(cel.fit$coef, primary=f) # is primary the coef?
        dev.off()

        pdf(paste0(contrast.group[f],".",labels,".volvano.tstat-bayes.limma.pdf"))
        volcanoplot(cel.fit,coef=f,highlight=10)
        dev.off()

# export significant genes and charts
# F-test p-values rather than t-test p-values (NESTEDF)
# Benjamini and Hochberg’s method to control the false discovery rate (GLOBAL)
        index <- c("up", "down")
        pvalue <- c(0.01, 0.001, 0.0001, "global", "nestedF")

        for (i in index) {
            for (p in pvalue) {
                pdf(paste0(contrast.group[f],".",labels,".",i,".P",
                           p,".venn.tstat-bayes.limma.pdf"))
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



##################################
####  First samples grouping  ####
##################################
## Samples classification and experimental designs
#  filter(!PATIENT_ID %in% c(paste0("CNR800",1:4),paste0("CNR900",1:2))) %>% # remove controls
metadata <- read.table("summary/phenodata", sep = "\t", header = T) %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction, ABClikelihood) %>%
    filter(Timepoint != "T2") %>%
    mutate(Relapse = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "R",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NR",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "NR",
                               TRUE ~ "C")) %>%
    mutate(ABClassify = case_when(ABClikelihood >= .9 ~ 1,
                                  ABClikelihood <= .1 ~ 0,
                                  TRUE ~ 2)) %>%
    mutate(Lymphnodes = case_when(SITE == "LN" ~ 1, TRUE ~ 0)) %>%
    mutate(Nodes = case_when(SITE == "LN" ~ "LN", TRUE ~ "EN"))

sink("metadata.3groups.samples.info.txt")
summary(metadata)
sink()

# make sure all samples have correct IDs accross all the array
metadata$ABClassify <- as.factor(metadata$ABClassify)
metadata$Lymphnodes <- as.factor(metadata$Lymphnodes)
metadata <- metadata[metadata$SAMPLE_ID %in% ids$V1, ]
row.names(metadata) = metadata$SAMPLE_ID
colnames(cel.raw) = metadata$SAMPLE_ID
pd <- AnnotatedDataFrame(data=metadata)
sampleNames(pd) <- metadata$SAMPLE_ID
phenoData(cel.raw) <- pd

# Robust Multi-Chip average for background correction, normalization
# expressions are log2 transformed
# RMA normalization
trx.normalized <- oligo::rma(cel.raw, target='core')
write.exprs(trx.normalized, file="normalized.trx.expression.txt")
#probe.normalized <- oligo::rma(cel.raw, target='probeset')
#write.exprs(probe.normalized, file="normalized.probe.expression.txt")
gc()

# Pull affymetrix annotations for genes and exons
featureData(trx.normalized) <- getNetAffx(trx.normalized, 'transcript')
sink("annotation.gene.exon.affy.txt")
with(fData(trx.normalized), table(seqname, category))
sink()
gc()

# moderated t-statistics (of standard errors) and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
groups = c("relapseCOOFactorial", "multiGrpRelapse", "twoGrpNodes", "multiGrpCOO", "relapseNodesFactorial")

for (g in groups) {

    if (g == "multiGrpRelapse") {
        design <- model.matrix(~ -1 + metadata$Relapse)
        colnames(design) <- c("noRelapse", "relapse", "control")
        contrast.matrix <- makeContrasts(noRelapse-relapse,
                                         noRelapse-control,
                                         relapse-control,
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)
            
    } else if (g == "relapseCOOFactorial") {
        sample.factors <- paste(metadata$Prediction, metadata$Relapse, sep=".")
        sample.factors <- factor(sample.factors,
                                 levels = c("ABC.NR", "GCB.NR", "U.NR",
                                            "ABC.R", "GCB.R", "U.R",
                                            "ABC.C", "GCB.C"))
        design <- model.matrix(~0 + sample.factors)
        colnames(design) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(Relapse2ABC=ABC.R-ABC.NR,
                                         Relapse2GCB=GCB.R-GCB.NR,
                                         ABC2GCB=(ABC.R-ABC.NR)-(GCB.R-GCB.NR),
                                         levels=design)
        coef <- rep(1:3) # refrence to each contrast
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)
        
    } else if (g == "relapseNodesFactorial") {
        sample.factors <- paste(metadata$Nodes, metadata$Relapse, sep=".")
        sample.factors <- factor(sample.factors,
                                 levels = c("LN.NR", "EN.NR",
                                            "LN.R", "EN.R",
                                            "LN.C", "EN.C"))
        design <- model.matrix(~0 + sample.factors)
        colnames(design) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(Relapse2LN=LN.R-LN.NR,
                                         Relapse2EN=EN.R-EN.NR,
                                         LN2EN=(LN.R-LN.NR)-(EN.R-EN.NR),
                                         levels=design)
        coef <- rep(1:3) # refrence to each contrast
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)

    } else if (g == "multiGrpCOO") {
        design <- model.matrix(~ -1 + metadata$ABClassify)
        colnames(design) <- c("GCB", "ABC", "Redundant")
        contrast.matrix <- makeContrasts(GCB-ABC,
                                         GCB-Redundant,
                                         ABC-Redundant,
                                         levels=design)
        coef <- rep(1:ncol(design)) # refrence to each contrast
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)
        
    } else if (g == "twoGrpNodes") {
        design <- model.matrix(~ -1 + metadata$Lymphnodes)
        colnames(design) <- c("Extranodal", "Lymphnodes")
        contrast.matrix <- makeContrasts(Extranodal-Lymphnodes, levels=design)
        coef <- c(1) # refrence to each contrast
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)
        
    }
}


###################################
####  Second samples grouping  ####
###################################
## Samples classification and experimental designs
metadata <- read.table("summary/phenodata", sep = "\t", header = T) %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction, ABClikelihood) %>%
    filter(Timepoint != "T2") %>%
    mutate(Relapse = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "R",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NR",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "S",
                               TRUE ~ "C")) %>%
    mutate(ABClassify = case_when(ABClikelihood >= .9 ~ 1,
                                  ABClikelihood <= .1 ~ 0,
                                  TRUE ~ 2)) %>%
    mutate(Lymphnodes = case_when(SITE == "LN" ~ 1, TRUE ~ 0)) %>%
    mutate(Nodes = case_when(SITE == "LN" ~ "LN", TRUE ~ "EN"))

sink("metadata.4groups.samples.info.txt")
summary(metadata)
sink()

# make sure all samples preserve their ID
metadata$ABClassify <- as.factor(metadata$ABClassify)
metadata$Lymphnodes <- as.factor(metadata$Lymphnodes)
metadata <- metadata[metadata$SAMPLE_ID %in% ids$V1, ]
row.names(metadata) = metadata$SAMPLE_ID
colnames(cel.raw) = metadata$SAMPLE_ID
pd <- AnnotatedDataFrame(data=metadata)
sampleNames(pd) <- metadata$SAMPLE_ID
phenoData(cel.raw) <- pd
gc()

# RMA normalization
trx.normalized <- oligo::rma(cel.raw, target='core')
write.exprs(trx.normalized, file="normalized.systemic.trx.expression.txt")
#probe.normalized <- oligo::rma(cel.raw, target='probeset')
#write.exprs(probe.normalized, file="normalized.systemic.probe.expression.txt")
gc()

# Pull affymetrix annotations for genes and exons
featureData(trx.normalized) <- getNetAffx(trx.normalized, 'transcript')

# moderated t-statistics (of standard errors) and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
groups = c("systemicRelapse", "systemicCOOFactorial")

for (g in groups) {

    if (g == "systemicRelapse") {
        design <- model.matrix(~ -1 + metadata$Relapse)
        colnames(design) <- c("noRelapse", "relapse", "systemic", "control")
        contrast.matrix <- makeContrasts(noRelapse-relapse,
                                         relapse-systemic,
                                         noRelapse-systemic,
                                         relapse-control,
                                         noRelapse-control,
                                         systemic-control,
                                         levels=design)
        coef <- rep(1:6)
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)
            
    } #else if (g == "systemicCOOFactorial") {

}

sample.factors <- paste(metadata$Prediction, metadata$Relapse, sep=".")
        sample.factors <- factor(sample.factors,
                                 levels = c("ABC.NR", "GCB.NR", "U.NR",
                                            "ABC.R", "GCB.R", "U.R",
                                            "ABC.C", "GCB.C",
                                            "ABC.S", "GCB.S"))
        design <- model.matrix(~0 + sample.factors)
        colnames(design) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(Relapse2ABC2S=ABC.R-ABC.S,
                                         Relapse2GCB2S=GCB.R-GCB.S,
                                         ABC2GCB2S=(ABC.R-ABC.S)-(GCB.R-GCB.S),
                                         levels=design)
        coef <- rep(1:3) # refrence to each contrast
#        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, pval=.001, coef=coef, percent=.15)
        
#    }
#}

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



        cel.fit <- lmFit(trx.normalized, design) %>%
            contrasts.fit(contrast.matrix) %>%
            eBayes()
        
        topTable(cel.fit, coef=f, adjust="fdr", sort.by="B",
                 number=selected) %>%
            write.table(file=paste0(contrast.group[f],
                                    ".moderated-tstat-bayes.limma.","systemicCOOFactorial",".txt"),
                        row.names=FALSE, sep="\t", quote=FALSE) 
