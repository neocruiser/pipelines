pkgs <- c('RColorBrewer', 'gplots',
          'dplyr', 'cn.mops',
          'doSNOW', 'lattice', 'vegan',
          'reshape2', 'ggplot2', 'tidyr', 'plyr')

lapply(pkgs, require, character.only = TRUE)


## load baited capture exon and dna segments
chromosome.segments <- read.table("padded.clean.bed",sep="\t",as.is=TRUE)
genomic.ranges <- GRanges(chromosome.segments[,1],
                          IRanges(chromosome.segments[,2],
                                  chromosome.segments[,3]))

## load aligned read to refernce bam files
raw.bam <- list.files(pattern=".bam$")
bam.ranges <- getSegmentReadCountsFromBAM(raw.bam,
                                          GR = genomic.ranges,
                                          sampleNames = paste0("S_",1:2),
                                          parallel = 32)


## cnv detection and segmentation
## https://www.rdocumentation.org/packages/cn.mops/versions/1.18.0/topics/exomecn.mops
## read counts are scaled within samples
## segmentation algorithms: fast or DNAcopy
cnv.estimation <- calcIntegerCopyNumbers(exomecn.mops(bam.ranges,
                                                      normType = "poisson",
                                                      segAlgorithm = "fast",
                                                      parallel = 32))


## plotting
## Segment means in Red
## Zeroline in grey
## plots configured based on DNAcopy R package
pdf("cnvs.by.chromosome.pdf", onefile = TRUE)
for(chr in c(seq(1:22), "X", "Y")) {
    sample.plots <- segplot(cnv.estimation, sampleIdx = 1, seqnames = chr)
    print(sample.plots)
}
dev.off()
try(dev.off(), silent = TRUE)


pdf("cnvs.by.samples.pdf", onefile = TRUE)
for(samp in 1:2){
    sample.plots <- segplot(cnv.estimation, sampleIdx = samp, plot.type = "w",
                            pt.cols = c("#202556", "#A52A2A"),
                            segcol = "#007BA7")
    print(sample.plots)
}
dev.off()
try(dev.off(), silent = TRUE)


## data export
as.data.frame(segmentation(cnv.estimation)) %>%
    write.table("cnv.segmentation.cnmops.txt", sep = '\t', row.names = FALSE, quote = FALSE)

as.data.frame(cnvs(cnv.estimation)) %>%
    write.table("cnv.samples.cnmops.txt", sep = '\t', row.names = FALSE, quote = FALSE)

as.data.frame(cnvr(cnv.estimation)) %>%
    write.table("cnv.regions.cnmops.txt", sep = '\t', row.names = FALSE, quote = FALSE)
