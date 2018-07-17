

## https://stackoverflow.com/questions/25781284/simplest-way-to-plot-changes-in-ranking-between-two-ordered-lists-in-r 
plotRanks <- function(a, b, labels.offset=0.15, arrow.len=0.05)
  {
  old.par <- par(mar=c(1,1,1,1), cex = .25)

  # Find the length of the vectors
  len.1 <- length(a)
  len.2 <- length(b)

  # Plot two columns of equidistant points
  plot(rep(1, len.1), 1:len.1, pch=20, cex=0.8, 
       xlim=c(0, 3), ylim=c(0, max(len.1, len.2)),
       axes=F, xlab="", ylab="") # Remove axes and labels
  points(rep(2, len.2), 1:len.2, pch=20, cex=0.8)

  # Put labels next to each observation
  text(rep(1-labels.offset, len.1), 1:len.1, a)
  text(rep(2+labels.offset, len.2), 1:len.2, b)

  # Now we need to map where the elements of a are in b
  # We use the match function for this job
  a.to.b <- match(a, b)

  # Now we can draw arrows from the first column to the second
  arrows(rep(1.02, len.1), 1:len.1, rep(1.98, len.2), a.to.b, 
         length=arrow.len, angle=20)
  par(old.par)
  }

ids2description <- read.table("./ids2description.summary.txt", header = FALSE, fill = TRUE)
colnames(ids2description) <- c("genes", "chromosome", "ensembl", "symbol", "function", "site", "symbol2")

lmfit.ids <- read.table("summary/ranking.lmfit", header = TRUE)
ranking.ids <- lmfit.ids[ lmfit.ids$transcriptclusterid %in% selgenes[[1]], 1, drop = FALSE] %>%
    rename(genes = transcriptclusterid)

ranking.adjpval <- ids2description[ids2description$genes %in% selgenes[[1]], c(1,4)] %>%
    mutate(ids = paste0(genes,"-",symbol)) %>%
    right_join(ranking.ids, by = "genes")


ranking.modules <- read.table("summry/ranking.degree", header = FALSE)

ranking.ids2 <- ranking.modules[ ranking.modules$V2 %in% selgenes[[1]], , drop = FALSE] %>%
    rename(modules = V1) %>%
    rename(genes = V2)

ranking.degree <- ids2description[ids2description$genes %in% selgenes[[1]], c(1,4)] %>%
    mutate(ids = paste0(genes,"-",symbol)) %>%
    right_join(ranking.ids2, by = "genes")

pdf("line.ranking.fdradjpval_degree.pdf")
plotRanks(rev(ranking.adjpval[,3]), rev(ranking.degree[,3]))
dev.off()

