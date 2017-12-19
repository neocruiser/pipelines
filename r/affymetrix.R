# Load packages
pkgs <- c('RColorBrewer', 'oligo', 'genefilter', 'pd.hta.2.0', 'gcrma', 'limma', 'affy', 'affyPLM', 'simpleaffy', 'affyQCReport', 'plier', "affycoretools", 'pdInfoBuilder', 'pvclust', 'xps')
lapply(pkgs, require, character.only = TRUE)

# Choose charts colors
palette.gr <- brewer.pal(11, name = "PRGn")
palette.rd <- brewer.pal(11, name = "RdYlBu")
palette.green <- colorRampPalette(palette.gr)(n = 200)
palette.red <- colorRampPalette(palette.rd)(n = 200)

# Microarray files loaded into array
# Robust Multi-Chip average for background correction, normalization
# expressions are log2 transformed
cel.files <- list.celfiles()
cel.raw <- read.celfiles(cel.files)
cel.normalized <- rma(cel.raw)
write.exprs(cel.normalized,file="data.txt")


#pdf("bootstrap.pdf")
#plot(pvData, hang=-1)
#pvrect(pvData, alpha=a)
#dev.off()

data:bgSequence(data:pmSequence)
mmSequence(data:pmSequence)

