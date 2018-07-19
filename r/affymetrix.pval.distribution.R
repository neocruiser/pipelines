pkgs <- c('RColorBrewer',
          'pvclust', 'gplots',
          'plyr', 'dplyr', 'reshape', 'tidyr',
          'vegan', 'ggplot2')
lapply(pkgs, require, character.only = TRUE)


## load summaries
ids <- read.table("summary/sampleIDs")
metadata.raw <- read.table("summary/phenodata", sep = "\t", header = TRUE)
abc.gcb <- read.table("summary/abc_gcb.RMA.txt", header = TRUE, fill = TRUE)
pvals <- read.table("summary/summary.lmfit.all.txt", header = TRUE, fill = TRUE)
abc.gcb.meta <- read.table("summary/abc_gcb.genes.counts", header = TRUE, fill = TRUE)
abc.gcb.folds <- read.table("summary/abc_gcb.lmFolds.txt", header = TRUE, fill = TRUE)


grouping = c("systemicRelapse", "systemicRelapseNodes", "systemicRelapseCOOprediction")


## PLOT 1
## plot pre-selected genes known to be preferentially expressed in either
## ABC or GCB cell of origin subtype
## combine gene serial number with gene symbol
dm <- left_join(abc.gcb, abc.gcb.meta, by = "ID")
rownames(abc.gcb) <- paste0(dm$ID,"-", dm$symbol)

## restructure the dataset
metadata <- metadata.raw %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction) %>%
    filter(Timepoint != "T2") %>%
    mutate(Groups = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "CNS",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NOREL",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "SYST",
                               TRUE ~ "CTRL"))


## factorize
metadata <- arrange(metadata, factor(SAMPLE_ID, levels = ids$V1))
metadata$Groups <- as.factor(metadata$Groups)
y <- metadata$Prediction

df <- data.frame(y, t(abc.gcb[, -1])) %>%
    gather("id", "expression", 2:c(dim(abc.gcb)[1]+1))


pdf("boxplots.abc_gcb.RMA.pdf")
df %>%
    ggplot(aes(x = reorder(paste0(id), expression),
               y = expression,
               fill = y)) +
    geom_boxplot(outlier.colour = NA, lwd = .1) +
    coord_flip() +
    scale_color_brewer(palette="Dark2") +
    theme_minimal() +
    theme(legend.position = "top",
          text = element_text(size = 7),
          axis.text.y = element_text(size = rel(.5))) +
    ggtitle(paste0("Intersection between genes")) +
    xlab("") +
    ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")
dev.off()
try(dev.off(), silent = TRUE)



## Plot 2
## plot ABC GCB preselected genes based on lmfit limma eBayes analysis
abc.gcb.folds$LogFC <- as.numeric(as.character(abc.gcb.folds$LogFC))

pdf("barplots.abc_gcb.lmfit.folds.pdf", onefile = TRUE)
for ( g in grouping ) {
    abc.bars <- abc.gcb.folds %>%
	filter(Contrast == g) %>%
        ggplot(aes(x = reorder(paste0(ID,"-",Symbol), LogFC),
                   y = LogFC,
                   fill = factor(Comparison))) +
        geom_bar(stat = "identity",
                 position = "dodge") +
        coord_flip() +
        scale_fill_brewer(palette="Set1") +
        theme_minimal() +
        theme(legend.position = "top",
              text = element_text(size = 7),
              axis.text.y = element_text(size = rel(1.5)),
              strip.background = element_rect(linetype = "blank",
                                              fill = "white"),
              panel.border = element_rect(linetype = "blank",
                                          fill = NA),
              panel.grid.major = element_line(linetype = "blank")) +
        ggtitle(paste0("Expression of ABC vs GCB genes (QC) for\n",g)) +
        xlab("") +
        ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")

    print(abc.bars)
}
dev.off()
try(dev.off(), silent = TRUE)



## PLOT 3
## density and histogram plots using p-values from lmfit adjusted with FDR comaprisons
pdf(paste0("pval.distribution.postRMA.postClean.affymetrix.pdf"), onefile = TRUE)

for (g in grouping) {
    for (p in c(8,9)) {

        dc <- pvals %>%
            filter(Contrast == g)
        density.plot <- dc %>%
            ggplot(aes(x = dc[, p],
                       fill = Comparison)) +
            geom_density() +
            theme_minimal() +
            theme(legend.position = "top") +
            scale_fill_brewer(palette = "Paired") +
            labs(x = names(pvals)[p],
                 y = "Density")


        dc <- pvals %>%
            filter(Contrast == g)
        histo.plot <- dc %>%
            ggplot(aes(x = dc[, p],
                       fill = Comparison)) +
            geom_histogram(binwidth = 0.02,
                           col=I("white")) +
            theme_minimal() +
            facet_wrap( ~ Comparison,
                       ncol = 2, scales = "free") +
            theme(legend.position = "none") +
            scale_fill_brewer(palette = "Paired") +
            labs(x = names(pvals)[p],
                 y = "Density")

        print(density.plot)
        print(histo.plot)
        
    }    
}

dev.off()
try(dev.off(), silent = TRUE)

## PLOT 4
## Venn diagrams
pvalues <- c(10^seq(-10,-2, length = 9), 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
pdf("venn.contrasts.lmfit.pdf", onefile = TRUE)

par(mfrow=c(3,2), cex = .4)
du <- NULL
for (g in grouping) {
    for (p in pvalues) {

        df <- pvals %>%
            filter(Contrast == g) %>%
            filter(FDRadjPval <= p)

        selgenes = NULL

        ## venn w more than 5 intersections are not plotted
        if ( dim(df)[1] > 0 ) {

                    print(dim(df))
                    du <- unique(df$Comparison)
                    print(as.character(du))

            for (i in 1:length(du)) {

                if ( i == 1 ){
                    selgenes <- list(df[ df$Comparison %in% du[[i]] , 3])
                    names(selgenes)[i] <- as.character(du[[i]])
                } else {
                    selgenes <- c(selgenes, list(df[ df$Comparison %in% du[[i]] , 3]))
                    names(selgenes)[i] <- as.character(du[[i]])
                }
            }
        } else { cat(">> Warning, no intersections found!\n") }
        
        try(venn(selgenes), silent = TRUE)
        try(title(paste0(g, '\n', "@maximum FDR-adjusted p-value of ", p)), silent = TRUE)
    }
}
dev.off()
try(dev.off(), silent = TRUE)



