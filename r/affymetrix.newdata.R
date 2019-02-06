# Load packages
pkgs <- c('RColorBrewer', 'genefilter', 'limma',
          'pvclust', 'foreach', 'oligo', 'pd.hta.2.0',
          'plyr', 'dplyr', 'reshape', 'tidyr', 'doMC',
          'vegan', 'WGCNA', 'readr', 'hgu219.db')
lapply(pkgs, require, character.only = TRUE)

## http://web.mit.edu/~r/current/lib/R/library/AnnotationDbi/doc/AnnotationDbi.pdf
## number of probes 49363
mapped.probes <- mappedkeys(hgu219ACCNUM)
## number of annotated genes 123135
aliases <-  as.list(hgu219ALIAS2PROBE)
symbols <-  as.list(hgu219SYMBOL)

summary(hgu219SYMBOL)
head(toTable(hgu219ALIAS2PROBE))
toTable(hgu219ALIAS2PROBE)[1:5, ]

cel.raw <- list.celfiles(full=TRUE, listGzipped=FALSE) %>%
    read.celfiles()
length(sampleNames(cel.raw))
colnames(cel.raw)


## load new metadata
ids <- read.table("metadata.spain", header = T) %>%
    mutate(Groups = case_when(groups %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "CNS",
                               groups %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NOREL",
                               groups == "SYTEMIC_RELAPSE_NO_CNS" ~ "SYST",
                              TRUE ~ "CTRL"))


## gene expression normalization and background correction
trx.normalized <- oligo::rma(cel.raw, background=TRUE, normalize=TRUE)
## backup all array expressions
write.table(trx.normalized, "expressions.normalized.spain", sep = "\t", quote = FALSE)
## searching normalized data
as.matrix(trx.normalized)["11728189_a_at", ]
geneid <- c("CXCR4", "NFAM1", "TRDJ1", "ALDH1A3")
toTable(subset(hgu219ALIAS2PROBE, Rkeys = geneid))


## commented out because no linear fitting was used
## moderatedFit <- function(data=trx.normalized, contrasts=contrast.matrix, labels="unsetGroup", coef=coef, percent=.15, group=g){
## # Fit Bayesian model and extract Differential Genes (sorted by significance)
## # Benjamini & Hochberg (1995): adjusted "fdr"
## # get the description of the two sample groups being compared
##     contrast.group <- gsub(" ","",colnames(contrasts)) # diana 3atetne el wa7e
    
## # get a proportion of differentially expressed genes
##     selected <- round(dim(data)[[1]] * percent)

##     for (f in coef) {

##         ## set proportion of differentially expression genes to 1% of the genome
##         cel.fit <- lmFit(data, strategy) %>%
##             contrasts.fit(contrasts) %>%
##             eBayes(proportion=0.01)
##         ## Benjamini & Hochberg (1995) control the false discovery rate,
##         ## the expected proportion of false discoveries amongst the rejected hypotheses
##         ## p-values are adjusted for multiple testing (fdr-adjusted p-value)
##         ## logFC = log2-fold-change
##         topTable(cel.fit, coef=f, adjust="fdr", sort.by="B", number=selected) %>%
##             write.table(file=paste0(contrast.group[f],".",labels,
##                                     ".moderated-tstat-bayes.limma.txt"),
##                         row.names=TRUE, sep="\t", quote=FALSE) 

##     }

##     sink(paste0(group,".regression.OK")); sink()
    
## }

## strategy <- model.matrix(~ -1 + ids$Groups)
## colnames(strategy) <- c("CNS", "NOREL", "SYST")        
## contrast.matrix <- makeContrasts(CNSvsNOREL = CNS-NOREL,
##                                  SYSTvsNOREL = SYST-NOREL,
##                                  CNSvsSYST = CNS-SYST,
##                                  levels = strategy)
## coef <- rep(1:3)
## moderatedFit(data=trx.normalized, contrasts=contrast.matrix, labels="systemicRelapse", coef=coef, percent=1)







####################
## classification ##
####################
## To assess missing genes in new data
## get selected gene panel after running LASSO
colnames(dat) = gsub(".TC.*$","",colnames(dat))
write.table(colnames(dat[, -1]), "panel.gene", quote=F, sep="\t")

################
## affymetrix ##
################
### remove gene aliases with special characters
gene.names <- read.table("panel.genes", header = T)
#cp classification/437017/reports/log.gene.names_bestLambda.txt
#gene.names <- read.table("log.gene.names_bestLambda.txt", header = T)

## get panel ids from array
all <- NULL
for (gen in gene.names$x){
    fd <- toTable(subset(hgu219ALIAS2PROBE, Rkeys = gen))
    all <- rbind(all, fd)
    return
}

## get GENE expressions from array
## combine probe id and gene expression
## collapse same gene names
expression <- as.matrix(trx.normalized)[all[, 1], ]
expression <- data.frame(probe_id = rownames(expression), expression)
expression <- full_join(expression, all, by = "probe_id")
expression <- expression[, -1]
dat <- aggregate(expression[, 1:26], list(expression$alias_symbol), mean)
rownames(dat)  <- dat[, 1]
dat <-  dat[ , -1]
means <- t(dat)
means <- data.frame(y = ids$Groups, means)
dim(means)
means[1:5,1:5]

## export expressions
write.table(means, "expression.spain", sep = "\t", quote = FALSE)



####################
## classification ##
####################
## import training data
new.data <- read.table("expression.spain", header=T)
## find missing genes
panel <- list(test = colnames(dat[,-1]), new = colnames(new.data[, -1]))
panel.intersect <- attr(venn(panel), "intersections")

## get missing genes
## run imputation
panel.missing <- dat[, colnames(dat) %in% panel.intersect$test]
panel.missing <- data.frame(y = dat[, 1], panel.missing)
panel.missing <- aggregate(panel.missing[, 2:dim(panel.missing)[2]], list(panel.missing$y), mean) %>%
    rename_at("Group.1", ~"y")



## full data frame with new imputated missing genes
## new.df will be used in the classification function modelTune.clas
new.df <- full_join(new.data, panel.missing, by = "y")
row.names(new.df)  <- row.names(new.data)








