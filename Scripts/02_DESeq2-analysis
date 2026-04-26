## STEP 3 - DESEQ2
library(DESeq2)
library(ggplot2)
library(pheatmap)

## set species factor and reference
coldata$species <- factor(coldata$species)
coldata$species <- relevel(coldata$species, ref = "Ab")

##create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_mat_filt,
  colData = coldata,
  design = ~ species
)

## quick normalization check
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

## variance stabilizing transform for QC
vsd <- vst(dds, blind = TRUE)



library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(apeglm)


## ensure species is a factor
coldata$species <- factor(coldata$species)

## initial DESeq object with Ab as reference
coldata$species <- relevel(coldata$species, ref = "Ab")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat_filt,
  colData = coldata,
  design = ~ species
)

## optional normalization / QC
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

vsd <- vst(dds, blind = TRUE)

## PCA
pcaData <- plotPCA(vsd, intgroup = "species", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = species)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()

## sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
sample_annot <- as.data.frame(colData(vsd)[, "species", drop = FALSE])

pheatmap(
  sampleDistMatrix,
  annotation_col = sample_annot,
  annotation_row = sample_annot
)

## run DESeq for Ab-reference model
dds <- DESeq(dds, test = "Wald")

resultsNames(dds)


## STEP 4- MAKE PAIRWISE RESULTS TABLE and VOLCANO PLOST
## species_a = numerator
## species_b = denominator/reference
make_pairwise_res <- function(dds, species_a, species_b, alpha = 0.05) {
  
  dds_tmp <- dds
  dds_tmp$species <- relevel(dds_tmp$species, ref = species_b)
  dds_tmp <- DESeq(dds_tmp, test = "Wald")
  
  coef_name <- paste0("species_", species_a, "_vs_", species_b)
  
  if (!coef_name %in% resultsNames(dds_tmp)) {
    stop(paste("Coefficient not found:", coef_name))
  }
  
  res_shrunk <- lfcShrink(
    dds_tmp,
    coef = coef_name,
    type = "apeglm"
  )
  
  res_tbl <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene_name") %>%
    as_tibble() %>%
    mutate(
      comparison = paste0(species_a, "_vs_", species_b),
      sig = case_when(
        !is.na(padj) & padj < alpha & log2FoldChange > 0 ~ paste(species_a, "higher"),
        !is.na(padj) & padj < alpha & log2FoldChange < 0 ~ paste(species_b, "higher"),
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  deg_tbl <- res_tbl %>%
    filter(!is.na(padj), padj < alpha)
  
  list(
    dds = dds_tmp,
    coef_name = coef_name,
    results = res_tbl,
    degs = deg_tbl
  )
}
