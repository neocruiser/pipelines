# Generate expression matrix from Affymetrix HTA 2.0 arrays.
# n=240 arrays
# This script should be run from the ~/dlbcl_cns_relapse folder.
# Author: Robert Kridel <robert.kridel@gmail.com>

library("oligo")
library("genefilter")
library("ggplot2")
library("plyr")
library("dplyr")

cel.files <- list.celfiles("Affy/cel", full.names = TRUE)

# Read in cel files
eset <- read.celfiles(cel.files)

# Read in phenodata
phenodata <- read.table(file = "meta_data/phenodata.txt", sep = "\t", header = T)

phenodata <- read_tsv(file = "meta_data/phenodata.txt") %>%
  dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Prediction, IPI_GROUP) %>%
  mutate(GROUP_NEW = ifelse(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                         "CNS_RELAPSE_CHOPorEQUIVALENT",
                                         "CNS_DIAGNOSIS"), "CNS",
                     ifelse(GROUP == "TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE", GROUP))) %>%
  filter(!SAMPLE_ID %in% c("CNR6038T1", "CNR7010T1")) %>% as.data.frame() %>%
  filter(!GROUP %in% c("NORMAL_GCB_CONTROL", "NORMAL_ABC_CONTROL")) %>%
  filter(Timepoint != "T2") %>%
  mutate(SITE_NEW = ifelse(SITE == "LN", "LN", "EN"))

# QC of non-normalized data
# Need to filter first
ffun <- filterfun(pOverA(p = 0.05, A = 200), cv(a = 0.5, b = Inf))
filt <- genefilter(exprs(eset), ffun)
eset_filt <- eset[filt, ]

# PCA of non-normalized data
exp_eset_filt <- log2(exprs(eset_filt))
PCA_raw <- prcomp(t(exp_eset_filt), scale = FALSE)

PCA_raw_for_plot <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                               Group = phenodata$GROUP,
                               COO = phenodata$Prediction,
                               Sample = phenodata$SAMPLE_ID)

PCA_raw_for_plot %>%
  ggplot(aes(PC1, PC2, col = Group),
      main = "PCA plot of the raw data (log-transformed)") +
  geom_point()

# Boxplot of intensitied of non-normalized data
exp_eset_filt_for_boxplot <- stack(as.data.frame(exp_eset_filt))

exp_eset_filt_for_boxplot %>%
  ggplot() +
  geom_boxplot(aes(x = ind, y = values))

#---
# Normalize eset
#---
eset_norm <- rma(eset, target = "core")

saveRDS(eset_norm, "Affy/eset_norm.rds")

# QC of normalized data
# Need to filter first
ffun_2 <- filterfun(pOverA(p = 0.05, A = 8), cv(a = 0.05, b = Inf))
filt_norm <- genefilter(exprs(eset_norm), ffun_2)
eset_norm_filt <- eset_norm[filt_norm, ]

# PCA of non-normalized data
exp_eset_norm_filt <- log2(exprs(eset_norm_filt))
PCA_raw <- prcomp(t(exp_eset_norm_filt), scale = FALSE)

PCA_raw_for_plot <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                               Group = phenodata$GROUP,
                               COO = phenodata$Prediction,
                               Sample = phenodata$SAMPLE_ID)

PCA_raw_for_plot %>%
  ggplot(aes(PC1, PC2, col = Group),
         main = "PCA plot of the raw data (log-transformed)") +
  geom_point()

# Boxplot of intensitied of non-normalized data
exp_eset_norm_filt_for_boxplot <- stack(as.data.frame(exp_eset_norm_filt))

exp_eset_norm_filt_for_boxplot %>%
  ggplot() +
  geom_boxplot(aes(x = ind, y = values))
