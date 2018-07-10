pkgs <- c('RColorBrewer',
          'pvclust', 'gplots',
          'plyr', 'dplyr', 'reshape', 'tidyr',
          'vegan', 'ggplot2')
lapply(pkgs, require, character.only = TRUE)



pvals <- read.table("./summary/summary.lmfit.all.txt", header = TRUE, fill = TRUE)


grouping = c("systemicRelapse", "systemicRelapseNodes", "systemicRelapseCOOprediction")


## PLOT 1
## density and histogram plots using p-values from lmfit adjusted with FDR comaprisons
pdf(paste0("pval.distribution.postRMA.postClean.affymetrix.pdf"), onefile = TRUE)

for (g in grouping) {

    density.plot <- pvals %>%
	filter(Contrast == g) %>%
	ggplot(aes(x = FDRadjPval,
		   fill = Comparison)) +
	geom_density() +
	theme_minimal() +
	theme(legend.position = "top") +
	scale_fill_brewer(palette = "Paired") +
	labs(x = "False discovery rate adjusted P-values",
             y = "Density")


    histo.plot <- pvals %>%
	filter(Contrast == g) %>% 
	ggplot(aes(x = FDRadjPval,
		   fill = Comparison)) +
	geom_histogram(binwidth = 0.02,
                       col=I("white")) +
	theme_minimal() +
	facet_wrap( ~ Comparison,
                   ncol = 2, scales = "free") +
	theme(legend.position = "none") +
	scale_fill_brewer(palette = "Paired") +
	labs(x = "False discovery rate adjusted P-values",
             y = "Density")

    print(density.plot)
    print(histo.plot)
    
    
}

dev.off()


## PLOT 2
## Venn diagrams
minPval <- min(pvals$FDRadjPval)
pvalues <- c(0.01, 0.025, 0.05, 0.1, 0.25)
pdf("venn.contrasts.lmfit.pdf", onefile = TRUE)

par(mfrow=c(3,2), cex = .4)
for (g in grouping) {
    for (p in pvalues) {

        df <- pvals %>%
            filter(Contrast == g) %>%
            filter(FDRadjPval <= p)

        print(dim(df))
        du <- unique(df$Comparison)
        print(du)

        selgenes = NULL
        for (i in 1:length(du)) {

            if ( i == 1 ){
                selgenes <- list(df[ df$Comparison %in% du[[i]] , 3])
                names(selgenes)[i] <- as.character(du[[i]])
            } else {
                selgenes <- c(selgenes, list(df[ df$Comparison %in% du[[i]] , 3]))
                names(selgenes)[i] <- as.character(du[[i]])
            }
        }
        try(venn(selgenes), silent = TRUE)
        title(paste0(g, '\n', "@maximum FDR-adj p-value of ", round(p, 5)))
    }
}
dev.off()




## PLOT 3
## plot pre-selected genes known to be preferentially expressed in either
## ABC or GCB cell of origin subtype
abc.gcb <- read.table("summary/abc_gcb.RMA.txt", header = TRUE, fill = TRUE, row.names = 1)
abc.gcb <- read.table("abc_gcb.RMA.txt", header = TRUE, fill = TRUE, row.names = 1)

## restructure the dataset
metadata <- read.table("phenodata", sep = "\t", header = T) %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction) %>%
    filter(Timepoint != "T2") %>%
    mutate(Groups = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "CNS",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NOREL",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "SYST",
                               TRUE ~ "CTRL"))


# factorize
metadata$Groups <- as.factor(metadata$Groups)
y <- metadata$Prediction


df <- data.frame(y, t(abc.gcb)) %>%
    gather("id", "expression", 2:c(dim(abc.gcb)[1]+1))


pdf("boxplots.abc_gcb.RMA.pdf")
df %>%
    ggplot(aes(x = reorder(paste0(id), expression),
               y = expression)) +
    geom_boxplot(outlier.colour = NA, lwd = .1) +
    coord_flip() +
    scale_color_brewer(palette="Dark2") +
    facet_wrap( ~ y ) + 
    theme_minimal() +
    theme(legend.position = "top",
          text = element_text(size = 5),
          axis.text.y = element_text(size = rel(.5))) +
    ggtitle(paste0("Intersection between genes")) +
    xlab("") +
    ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")
dev.off()
try(dev.off(), silent = TRUE)
