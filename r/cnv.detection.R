pkgs <- c('RColorBrewer', 'gplots',
          'dplyr', 'cn.mops',
          'doSNOW', 'lattice', 'vegan',
          'reshape2', 'ggplot2', 'tidyr', 'plyr')

lapply(pkgs, require, character.only = TRUE)


## number of parallel threads
threads=32

## load baited capture exon and dna segments
chromosome.segments <- read.table("/cluster/projects/kridelgroup/relapse/mutations/targets/padded.clean.bed",sep="\t",as.is=TRUE)
genomic.ranges <- GRanges(chromosome.segments[,1],
                          IRanges(chromosome.segments[,2],
                                  chromosome.segments[,3]))

## load aligned read to refernce bam files
raw.bam <- list.files("/cluster/projects/kridelgroup/relapse/mutations/raw/new",
                      pattern=".realigned.bam$", full.names = TRUE)
raw.bam
##name.bam <- gsub("_.*$","",raw.bam)
##name.bam <- gsub("^.*A61","A61",name.bam)
name.bam <- gsub(".bam$","",raw.bam)
name.bam <- gsub("^.*A95","A95",name.bam)
name.bam

## read counts and stat aggregation
bam.ranges <- getSegmentReadCountsFromBAM(raw.bam,
                                          GR = genomic.ranges,
                                          sampleNames = paste0(name.bam),
                                          parallel = threads)


## cnv detection and segmentation
## https://www.rdocumentation.org/packages/cn.mops/versions/1.18.0/topics/exomecn.mops
## read counts are scaled within samples
## segmentation algorithms: fast or DNAcopy
cnv.estimation <- calcIntegerCopyNumbers(exomecn.mops(bam.ranges,
                                                      normType = "poisson",
                                                      segAlgorithm = "fast",
                                                      parallel = threads))


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
    write.table("output.2_r.cnmops.multi_segmentation.cnv.txt", sep = '\t', row.names = FALSE, quote = FALSE)

as.data.frame(cnvs(cnv.estimation)) %>%
    write.table("output.2_r.cnmops.multi_samples.cnv.txt", sep = '\t', row.names = FALSE, quote = FALSE)

as.data.frame(cnvr(cnv.estimation)) %>%
    write.table("output.2_r.cnmops.multi_regions.cnv.txt", sep = '\t', row.names = FALSE, quote = FALSE)
