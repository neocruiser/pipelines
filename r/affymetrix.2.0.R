# Load packages
pkgs <- c('RColorBrewer', 'genefilter', 'limma',
          'pvclust', 'foreach', 'oligo', 'pd.hta.2.0',
          'plyr', 'dplyr', 'reshape', 'tidyr', 'doMC',
          'vegan')
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
cel.raw <- list.celfiles("../../raw", full=TRUE, listGzipped=FALSE) %>%
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
moderatedFit <- function(data=trx.normalized, contrasts=contrast.matrix, labels="unsetGroup", coef=coef, percent=.15, group=g){
# Fit Bayesian model and extract Differential Genes (sorted by significance)
# Benjamini & Hochberg (1995): adjusted "fdr"
# get the description of the two sample groups being compared
    contrast.group <- gsub(" ","",colnames(contrasts)) # diana 3atetne el wa7e
    
# get a proportion of differentially expressed genes
    selected <- round(dim(data)[[1]] * percent)

    for (f in coef) {

        cel.fit <- lmFit(data, strategy) %>%
            contrasts.fit(contrasts) %>%
            eBayes()
        
        topTable(cel.fit, coef=f, adjust="fdr", sort.by="B",
                 number=selected) %>%
            write.table(file=paste0(contrast.group[f],".",labels,
                                    ".moderated-tstat-bayes.limma.txt"),
                        row.names=FALSE, sep="\t", quote=FALSE) 

    }

    sink(paste0(group,".regression.OK")); sink()
    
}



###################################
####  Sample grouping  ####
###################################
## Samples classification and experimental designs
metadata <- read.table("summary/phenodata", sep = "\t", header = T) %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Score, Prediction, ABClikelihood) %>%
    filter(Timepoint != "T2") %>%
    mutate(Groups = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "CNS",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NOREL",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "SYST",
                               TRUE ~ "CTRL")) %>%
    mutate(ABClassify = case_when(ABClikelihood >= .9 ~ "ABC",
                                  ABClikelihood <= .1 ~ "GCB",
                                  TRUE ~ "U")) %>%
    mutate(ABCScore = case_when(Score > 2412 ~ "ABC",
                                Score <= 1900 ~ "GCB",
#                                Score == NA ~ "NA",
                                TRUE ~ "U")) %>%
    mutate(Nodes = case_when(SITE == "LN" ~ "LN",
                             SITE == "TO" ~ "LN",
                             SITE == "SP" ~ "LN",
                             TRUE ~ "EN")) %>%
    mutate(Lymphnodes = case_when(Nodes == "LN" ~ 1, TRUE ~ 0))

# make sure all samples preserve their ID
metadata$Groups <- as.factor(metadata$Groups)
metadata$ABClassify <- as.factor(metadata$ABClassify)
metadata$ABCScore <- as.factor(metadata$ABCScore)
metadata$Nodes <- as.factor(metadata$Nodes)
metadata$Lymphnodes <- as.factor(metadata$Lymphnodes)
metadata <- metadata[metadata$SAMPLE_ID %in% ids$V1, ]
row.names(metadata) = metadata$SAMPLE_ID
colnames(cel.raw) = metadata$SAMPLE_ID
pd <- AnnotatedDataFrame(data=metadata)
sampleNames(pd) <- metadata$SAMPLE_ID
phenoData(cel.raw) <- pd
gc()


summary(metadata)
sink("metadata.samples.info.txt")
summary(metadata)
sink()


# RMA normalization
trx.normalized <- oligo::rma(cel.raw, target='core')
write.exprs(trx.normalized, file="normalized.systemic.trx.expression.txt")

probe.normalized <- oligo::rma(cel.raw, target='probeset')
write.exprs(probe.normalized, file="normalized.systemic.probe.expression.txt")
gc()



##### REMOVE HIGH VARIANCE GENES
if ( file.exists("./ids.wo.ncrna") ) {

    ## file already constructed based on RNA pattern occurence in the annotated array
    ids.wo.ncrna <- read.table("ids.wo.ncrna", header = FALSE)
    
    # REMOVE NCRNAS
    x <- as.matrix(trx.normalized[as.matrix(ids.wo.ncrna), ])
    colnames(x) <- colnames(trx.normalized)
    dim(x)
    sink("create.subset.wo.ncRNA.OK"); sink()
    
} else {
    
    x <- as.matrix(trx.normalized)
    dim(x)
    sink("create.subset.w.ncRNA.OK"); sink()

}

######################
## FUNCTION CALLING ##
######################
get.var <- function(dat, n, from = 1, to = (dim(dat)[2])*0.1, silent = FALSE, remove.hi = 1 ){
    ## GET THE RANGE OF VARIANCE ACROSS ALL THE DATASET
    locus.var <- apply(t(dat), n, var)

    if ( remove.hi == 1 ) {
        ## discard high variance genes
        hi.var <- order(abs(locus.var), decreasing = T)[from:to]
    } else if ( remove.hi == 0 ) {
        ## discard low variance genes
        hi.var <- order(abs(locus.var), decreasing = F)[from:to]
    }
    
    if ( silent == FALSE ) {
        cat("Number of selected high-variance genes:",length(hi.var),"\n")    
    }
    
    return(hi.x <- dat[,hi.var])
}


######################
## SUBSET SELECTION ##
######################
set.seed(15879284)

# data transformation. Difference between gene expressions
# has better variance interpretation
xs <- decostand(x, "standardize")

## get number of discarded high variance genes
## iterate multiple thresholds, maximum 20 iterations
gv <- NULL
start_th=floor( (nrow(xs) * 0.1) / 2)
end_th=nrow(xs)
increment_th=floor( (end_th - start_th) / 20 )

for (nset in seq(start_th, end_th, increment_th)) {
    ## UNSUPERVISED GENE SELECTION BASED ON HIGH VARIANCE
    ## TO DISCARD EXTREME VARIANCE
    hi.x <- get.var(t(xs), 1, from = 1, to = nset)

    # get mean (dm) and maximum value of variance (dmv) of the whole dataset
    dm_old<- summary(apply(hi.x, 2, var))[[4]]
    dmv_old<- summary(apply(hi.x, 2, var))[[6]]
    offset = c( (dm_old* (log(dmv_old)) + (dmv_old* 0.1)) )

    # get mean standard deviation (dms) and max SD (dsv)
    dms_old <- summary(apply(hi.x, 2, sd))[[4]]
    dsv_old <- summary(apply(hi.x, 2, sd))[[6]]

    # get the mean for each gene
    selected <- data.frame(locus = colnames(hi.x),
                           var = apply(hi.x, 2, var)) %>%
        filter(var > offset) %>%
        nrow

    # recalculate variance based on adujusted new thresholds
    hi.x <- get.var(t(xs), 1, from = c(selected+1), to = nset, silent = TRUE)
    dm_new <- summary(apply(hi.x, 2, var))[[4]]
    dmv_new <- summary(apply(hi.x, 2, var))[[6]]
    dms_new <- summary(apply(hi.x, 2, sd))[[4]]
    dsv_new <- summary(apply(hi.x, 2, sd))[[6]]
    gv <- rbind(gv, data.frame(dimension=nset,
                               meanVariance=dm_old,
                               maxVariance=dmv_old,
                               meanSD=dms_old,
                               maxSD=dsv_old,
                               discarded=selected,
                               adj.meanVariance=dm_new,
                               adj.maxVariance=dmv_new,
                               adj.meanSD=dms_new,
                               adj.maxSD=dsv_new))

}

write.table(gv, "summary.adjusted.means.subsetting.txt", quote=FALSE, sep="\t", row.names=F)


## subset the dataset based on a selected mean and SD
means2subset <- gv %>%
    filter(adj.meanVariance > 0.03 & adj.meanVariance <= 0.032) %>%
    select(dimension, discarded)

from.m=c(means2subset$discarded[[1]] + 1)
to.m=means2subset$dimension[[1]]

## dimension minus the discarded high variance genes
adj.x <- get.var(t(xs), 1, from = from.m, to = to.m, silent = FALSE, remove.hi = 0)


##   subset and selcet normal variance genes
trx.normalized <- trx.normalized[colnames(adj.x), ]
write.exprs(trx.normalized, file=paste0("normalized.subset.",to.m,".systemic.trx.expression.txt"))

# Pull affymetrix annotations for genes and exons
featureData(trx.normalized) <- getNetAffx(trx.normalized, 'transcript')

sink("subset.annotated.transcripts.OK"); sink()

# Test the matrix construction process
# get structure of the matrices
# example design
sample.factors <- paste(metadata$Groups, metadata$Nodes, sep=".")
sample.factors <- factor(sample.factors, levels=c(unique(sample.factors)))
strategy <- model.matrix(~0 + sample.factors)
colnames(strategy) <- levels(sample.factors)
contrast.matrix <- makeContrasts(CNSvsNOREL_LN = CNS.LN-NOREL.LN, levels=strategy)

cat("\n\n\n\nSummary Metadata:\n")
summary(metadata)
cat("\n\n\n\nAll possible permutations between sample cases:\n")
sample.factors           
cat("\n\n\n\nAll possible associations between pemutations and samples:\n")
strategy
cat("\n\n\n\nNumber of samples and number of permutations:\n")
dim(strategy)            
cat("\n\n\n\nDesign of the analysis and which cases to be used:\n")
contrast.matrix          
cat("\n\n\n\nNumber columns must be identical to the number of rows in trx.normalized:\n")
dim(contrast.matrix)     
cat("\n\n\n\nDimensions of the gene expressions x sample cases (trx.normalized):\n")
dim(trx.normalized)      



# moderated t-statistics (of standard errors) and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors
groups = c("systemicRelapse", "systemicRelapseNodes", "systemicRelapseCOOclasses", "systemicRelapseCOOscores", "systemicRelapseCOOprediction")

for (g in groups) {

    if (g == "systemicRelapse") {
        strategy <- model.matrix(~ -1 + metadata$Groups)
        colnames(strategy) <- c("CNS", "NOREL", "SYST", "CTRL")
        contrast.matrix <- makeContrasts(CNSvsNOREL = CNS-NOREL,
                                         SYSTvsNOREL = SYST-NOREL,
                                         CNSvsSYST = CNS-SYST,
                                         levels = strategy)
        coef <- rep(1:3)
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, coef=coef, percent=1)

        
    } else if (g == "systemicRelapseNodes") {
        sample.factors <- paste(metadata$Groups, metadata$Nodes, sep=".")
        sample.factors <- factor(sample.factors, levels=c(unique(sample.factors)))
        strategy <- model.matrix(~0 + sample.factors)
        colnames(strategy) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(CNSvsNOREL_LN = CNS.LN-NOREL.LN,
                                         CNSvsNOREL_EN = CNS.EN-NOREL.EN,
                                         SYSTvsNOREL_LN = SYST.LN-NOREL.LN,
                                         SYSTvsNOREL_EN = SYST.EN-NOREL.EN,
                                         CNSvsSYST_LN = CNS.LN-SYST.LN,
                                         CNSvsSYST_EN = CNS.EN-SYST.EN,
                                         diffCNSvsNOREL_LNvsEN = (CNS.LN-NOREL.LN)-(CNS.EN-NOREL.EN),
                                         diffSYSTvsNOREL_LNvsEN = (SYST.LN-NOREL.LN)-(SYST.EN-NOREL.EN),
                                         diffCNSvsSYST_LNvsEN = (CNS.LN-SYST.LN)-(CNS.EN-SYST.EN),
                                         levels = strategy)
        coef <- rep(1:9)
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, coef=coef, percent=1)
        
    } else if (g == "systemicRelapseCOOclasses") {
        ## experimental
        sample.factors <- paste(metadata$Groups, metadata$ABClassify, sep=".")
        sample.factors <- factor(sample.factors, levels=c(unique(sample.factors)))        
        strategy <- model.matrix(~0 + sample.factors)
        colnames(strategy) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(CNSvsNOREL_ABC = CNS.ABC-NOREL.ABC,
                                         CNSvsNOREL_GCB = CNS.GCB-NOREL.GCB,
                                         SYSTvsNOREL_ABC = SYST.ABC-NOREL.ABC,
                                         SYSTvsNOREL_GCB = SYST.GCB-NOREL.GCB,
                                         CNSvsSYST_ABC = CNS.ABC-SYST.ABC,
                                         CNSvsSYST_GCB = CNS.GCB-SYST.GCB,
                                         diffCNSvsNOREL_ABCvsGCB = (CNS.ABC-NOREL.ABC)-(CNS.GCB-NOREL.GCB),
                                         diffSYSTvsNOREL_ABCvsGCB = (SYST.ABC-NOREL.ABC)-(SYST.GCB-NOREL.GCB),
                                         diffCNSvsSYST_ABCvsGCB = (CNS.ABC-SYST.ABC)-(CNS.GCB-SYST.GCB),
                                         levels = strategy)
        coef <- rep(1:9)
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, coef=coef, percent=1)
        
    } else if (g == "systemicRelapseCOOprediction") {
        ## selected
        sample.factors <- paste(metadata$Groups, metadata$Prediction, sep=".")
        sample.factors <- factor(sample.factors, levels=c(unique(sample.factors)))        
        strategy <- model.matrix(~0 + sample.factors)
        colnames(strategy) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(CNSvsNOREL_ABC = CNS.ABC-NOREL.ABC,
                                         CNSvsNOREL_GCB = CNS.GCB-NOREL.GCB,
                                         SYSTvsNOREL_ABC = SYST.ABC-NOREL.ABC,
                                         SYSTvsNOREL_GCB = SYST.GCB-NOREL.GCB,
                                         CNSvsSYST_ABC = CNS.ABC-SYST.ABC,
                                         CNSvsSYST_GCB = CNS.GCB-SYST.GCB,
                                         diffCNSvsNOREL_ABCvsGCB = (CNS.ABC-NOREL.ABC)-(CNS.GCB-NOREL.GCB),
                                         diffSYSTvsNOREL_ABCvsGCB = (SYST.ABC-NOREL.ABC)-(SYST.GCB-NOREL.GCB),
                                         diffCNSvsSYST_ABCvsGCB = (CNS.ABC-SYST.ABC)-(CNS.GCB-SYST.GCB),
                                         levels = strategy)
        coef <- rep(1:9)
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, coef=coef, percent=1)
        
    } else if (g == "systemicRelapseCOOscores") {
        ## experimental
        sample.factors <- paste(metadata$Groups, metadata$ABCScore, sep=".")
        sample.factors <- factor(sample.factors, levels=c(unique(sample.factors)))        
        strategy <- model.matrix(~0 + sample.factors)
        colnames(strategy) <- levels(sample.factors)
        contrast.matrix <- makeContrasts(CNSvsNOREL_ABC = CNS.ABC-NOREL.ABC,
                                         CNSvsNOREL_GCB = CNS.GCB-NOREL.GCB,
                                         SYSTvsNOREL_ABC = SYST.ABC-NOREL.ABC,
                                         SYSTvsNOREL_GCB = SYST.GCB-NOREL.GCB,
                                         CNSvsSYST_ABC = CNS.ABC-SYST.ABC,
                                         CNSvsSYST_GCB = CNS.GCB-SYST.GCB,
                                         diffCNSvsNOREL_ABCvsGCB = (CNS.ABC-NOREL.ABC)-(CNS.GCB-NOREL.GCB),
                                         diffSYSTvsNOREL_ABCvsGCB = (SYST.ABC-NOREL.ABC)-(SYST.GCB-NOREL.GCB),
                                         diffCNSvsSYST_ABCvsGCB = (CNS.ABC-SYST.ABC)-(CNS.GCB-SYST.GCB),
                                         levels = strategy)
        coef <- rep(1:9)
        moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels=g, coef=coef, percent=1)
        
    }
}
warnings()



sink("r_session.info")
sessionInfo()
sink()

#####################
##### DEBUGGING #####
#####################
# save and load data. Best to set it after RMA
# save.image("../debug.RData")
# load("../debug.RData")

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

