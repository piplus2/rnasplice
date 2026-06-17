#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

# Load required packages
library(BiocParallel)
library(parallel)
library(DEXSeq)

BPPARAM <- MulticoreParam(workers = parallel::detectCores() - 1)

# Read command arguments
argv <- commandArgs(trailingOnly = TRUE)
argc <- length(argv)
if (argc < 5) {
  stop("Usage: run_dexseq_exon.R <counts_dir> <annotation_file> <samples_table> <contrasts_table> <ntop>")
}
# Parse command arguments
counts <- argv[1]
annotation <- argv[2]
samples <- argv[3]
contrasts <- argv[4]
ntop <- as.integer(argv[5])

#############################
## Define helper functions ##
#############################

splitByContrast <- function(object, contrast) {
  keep <- colData(object)$condition %in% contrast
  object <- object[, keep, drop = FALSE]
  colData(object)$condition <- factor(
    x = colData(object)$condition,
    levels = contrast
  )
  object
}

DEXSeq_pipeline <- function(object, BPPARAM = BiocParallel::bpparam()) {
  object <- estimateSizeFactors(object)
  object <- estimateDispersions(object, BPPARAM = BPPARAM)
  object <- testForDEU(object, BPPARAM = BPPARAM)
  object <- estimateExonFoldChanges(object, BPPARAM = BPPARAM)
  object
}

keepValidColumns <- function(x) {
  keep <- vapply(x, function(col) is.numeric(col) || is.character(col), FUN.VALUE = logical(1))
  x <- x[, keep, drop = FALSE]
  x
}

vectorToDataFrame <- function(x) {
  x <- data.frame(groupID = names(x), padj = unname(x))
  x
}

savePlotDEXSeq <- function(x, file, ntop = 10) {
  x_valid <- x[is.finite(x$pvalue) & is.finite(x$padj), , drop = FALSE]
   if (nrow(x_valid) < ntop) {
    ntop <- nrow(x_valid)
  }
  if (nrow(x_valid) == 0) {
    warning("No valid results to plot.")
    return()
  }

  ntop <- min(ntop, nrow(x_valid))
  ind <- order(x_valid$pvalue)[1:ntop]
  ids <- x_valid$groupID[ind]

  pdf(file = file, width = 11.7, height = 8.3)
  for (i in ids) {
    tryCatch(
      {
        plotDEXSeq(x, geneID = i, displayTranscripts = TRUE)
      },
      error = function(e) {
        warning(paste("Error plotting gene", i, ":", e$message))
      }
    )
  }
  dev.off()
}

# Read samples table
samples <- read.csv(samples, stringsAsFactors = TRUE, check.names = FALSE)
samples <- samples[, c("sample", "condition"), drop = FALSE]
samples <- unique(samples)

# Read contrasts table
contrasts <- read.csv(contrasts, check.names = FALSE)
contrasts <- contrasts[, c("contrast", "treatment", "control"), drop = FALSE]
contrasts <- unique(contrasts)

# List counts files
files <- paste0(counts, "/", samples$sample, ".clean.count.txt")

# Create DEXSeqDataSet object
dexseq_dataset <- DEXSeqDataSetFromHTSeq(
  countfiles = files,
  sampleData = samples,
  flattenedfile = annotation
)

# Define contrasts
groups <- levels(colData(dexseq_dataset)$condition)
contrast_names <- contrasts$contrast
contrasts <- data.frame(A = contrasts$treatment, B = contrasts$control)
contrasts <- contrasts[contrasts$A != contrasts$B, , drop = FALSE]
contrasts <- asplit(contrasts, MARGIN = 1)

# Create objects list
objects <- mapply(
  FUN = splitByContrast,
  contrast = contrasts,
  MoreArgs = list(object = dexseq_dataset),
  SIMPLIFY = FALSE
)

# Run DEXSeq workflow
objects <- lapply(objects, DEXSeq_pipeline, BPPARAM = BPPARAM)
names(objects) <- contrast_names

# Run DEXSeqResults
results_exon <- lapply(objects, DEXSeqResults)
names(results_exon) <- contrast_names

# Run perGeneQValue
results_gene <- lapply(results_exon, perGeneQValue)
names(results_gene) <- contrast_names

##########################
## Save objects to disk ##
##########################

mapply(
  saveRDS,
  object = objects,
  file = paste0("DEXSeqDataSet.", contrast_names, ".rds")
)

mapply(
  saveRDS,
  object = results_exon,
  file = paste0("DEXSeqResults.", contrast_names, ".rds")
)

mapply(
  saveRDS,
  object = results_gene,
  file = paste0("perGeneQValue.", contrast_names, ".rds")
)

##########################
## Save results to disk ##
##########################

mapply(
  write.csv,
  x = lapply(results_exon, keepValidColumns),
  file = paste0("DEXSeqResults.", contrast_names, ".csv"),
  MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
  write.csv,
  x = lapply(results_gene, vectorToDataFrame),
  file = paste0("perGeneQValue.", contrast_names, ".csv"),
  MoreArgs = list(quote = FALSE, row.names = FALSE)
)

########################
## Save plots to disk ##
########################

mapply(
  savePlotDEXSeq,
  x = results_exon,
  file = paste0("plotDEXSeq.", contrast_names, ".pdf"),
  MoreArgs = list(ntop = ntop)
)

###############################
## Print session information ##
###############################

citation("DEXSeq")

sessionInfo()
