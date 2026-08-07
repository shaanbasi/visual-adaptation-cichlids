library(DESeq2)
library(WGCNA)
library(dplyr)
library(tibble)

options(stringsAsFactors = FALSE)

input_dir  <- "Data/Processed"
output_dir <- "Data/Processed"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

## load filtered count matrix and sample metadata
count_mat <- readRDS(
  file.path(input_dir, "count_mat_filt.rds")
)

coldata <- readRDS(
  file.path(input_dir, "coldata.rds")
)

count_mat <- as.matrix(count_mat)
coldata <- as.data.frame(coldata)

## confirm that samples match
stopifnot(
  all(colnames(count_mat) %in% rownames(coldata))
)

## put metadata in the same order as the count matrix
coldata <- coldata[colnames(count_mat), , drop = FALSE]

stopifnot(
  identical(colnames(count_mat), rownames(coldata))
)

## create deseq2 object for variance stabilising transformation
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = coldata,
  design = ~ species
)

## remove extremely low-information genes for WGCNA
keep_wgcna <- rowSums(counts(dds) >= 10) >= 2

dds_wgcna <- dds[keep_wgcna, ]

## variance stabilising transformation
vsd_wgcna <- vst(
  dds_wgcna,
  blind = FALSE
)

## WGCNA requires rows = samples, columns = genes
datExpr <- t(assay(vsd_wgcna))
datExpr <- as.data.frame(datExpr)

## remove genes and samples failing wgcna quality checks
quality_check <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

if (!quality_check$allOK) {
  
  if (sum(!quality_check$goodGenes) > 0) {
    message(
      "Removing ",
      sum(!quality_check$goodGenes),
      " problematic genes."
    )
  }
  
  if (sum(!quality_check$goodSamples) > 0) {
    message(
      "Removing ",
      sum(!quality_check$goodSamples),
      " problematic samples."
    )
  }
  
  datExpr <- datExpr[
    quality_check$goodSamples,
    quality_check$goodGenes,
    drop = FALSE
  ]
}

##match metadata to retained samples
speciesTraits <- coldata[
  rownames(datExpr),
  ,
  drop = FALSE
]

##keep the original categorical species variable
speciesTraits$species <- factor(speciesTraits$species)

## save inputs required by run_wgcna
saveRDS(
  datExpr,
  file.path(output_dir, "datExpr.rds")
)

saveRDS(
  speciesTraits,
  file.path(output_dir, "speciesTraits.rds")
)
