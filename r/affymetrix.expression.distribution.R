pkgs <- c('RColorBrewer',
          'pvclust', 'gplots',
          'plyr', 'dplyr', 'reshape', 'tidyr',
          'vegan', 'ggplot2')
lapply(pkgs, require, character.only = TRUE)


## load summaries
ids <- read.table("summary/sampleIDs")
metadata.raw <- read.table("summary/phenodata", sep = "\t", header = TRUE)
pvals <- read.table("summary/summary.lmfit.all.txt", header = TRUE, fill = TRUE)
rma.all <- read.table("summary/normalized.subsetCleaned_GEN19477.systemic.trx.expression.txt", header = TRUE, fill = TRUE)


grouping = c("systemicRelapse", "systemicRelapseNodes", "systemicRelapseCOOprediction")


## PLOT 1
p=10^-1
genes.significant <- pvals %>%
    filter(Contrast == "systemicRelapseCOOprediction") %>%
    filter(FDRadjPval <= p)
dim(genes.significant)

gene.ids <- genes.significant[1:21, 3 ]
gene.pvals <- pvals[ pvals$ID %in% gene.ids,  ]
dim(gene.pvals)

pdf("bars.significant.genes.lmfit.pdf", onefile = TRUE)
for ( g in grouping ) {
    gene.bars <- gene.pvals %>%
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
              text = element_text(size = 6),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              strip.background = element_rect(linetype = "blank",
                                              fill = "white"),
              panel.border = element_rect(linetype = "blank",
                                          fill = NA),
              panel.grid.major = element_line(linetype = "blank")) +
        ggtitle(paste0("Differential genes for\n",g, " @FDR-adjusted p-val ", p)) +
        xlab("") +
        ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")

    print(gene.bars)
}
dev.off()



## PLOT 2
## plotting single gene RMA expression across sites and cell of origin
## restructure the dataset
metadata <- metadata.raw %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction) %>%
    filter(Timepoint != "T2") %>%
    mutate(Groups = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "CNS",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NOREL",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "SYST",
                              TRUE ~ "CTRL")) %>%
    filter(Groups != "CTRL")

metadata <- arrange(metadata, factor(SAMPLE_ID, levels = ids$V1))
metadata$Groups <- as.factor(metadata$Groups)
y <- metadata$Prediction


pdf("boxplots.per.gene.RMA.pdf", onefile = TRUE)
for (gene.name in gene.ids) {
    rma.selected <- rma.all %>%
        filter(row.names(rma.all) == gene.name)

    x <- as.data.frame(t(rma.selected))
    x$SAMPLE_ID <- rownames(x)
    df <-  right_join(x, metadata[, c(1,2:6)], by = "SAMPLE_ID")

    gene.boxes <- df %>%
        ggplot(aes(x = reorder(paste0(SITE), V1),
                   y = V1,
                   fill = Prediction)) +
        geom_boxplot(outlier.colour = NA, lwd = .1) +
        coord_flip() +
        scale_color_brewer(palette="Dark2") +
        theme_minimal() +
        facet_wrap( ~ GROUP + Groups,
                   ncol = 4, scales = "free") +
        theme(legend.position = "top",
              text = element_text(size = 7),
              axis.text.y = element_text(size = rel(.5))) +
        ggtitle(paste0(gene.name, " differentially expressed @FDR-adjusted p-val", p)) +
        xlab("") +
        ylab("RMA quantile normalization (2 is 4-fold up)")

    print(gene.boxes)
}
dev.off()
