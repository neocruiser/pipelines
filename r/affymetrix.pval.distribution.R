pkgs <- c('RColorBrewer',
          'pvclust', 'foreach',
          'plyr', 'dplyr', 'reshape', 'tidyr',
          'vegan', 'ggplot2')
lapply(pkgs, require, character.only = TRUE)



pvals <- read.table("./summary.lmfit.all.txt", header = TRUE, fill = TRUE)


groups = c("systemicRelapse", "systemicRelapseNodes", "systemicRelapseCOOprediction")

pdf(paste0("pval.distribution.postRMA.postClean.affymetrix.pdf"), onefile = TRUE)

for (g in groups) {

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
