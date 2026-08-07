library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(ggrepel)
library(apeglm)
library(grid)


## ANALYSIS SETTINGS
## wald-test significance threshold
alpha <- 0.05

## visual fold-change guide for volcano plots
## this is not used to define significant DE genes
volcano_lfc_guide <- 1

## number of genes labelled in each expression direction
labels_per_side <- 5

## species levels
species_levels <- c("Ab", "Mz", "Nb", "On", "Pn")

## folder containing prepared data
prep_output_dir <- "Data_prep_results"

## output folder
output_dir <- "DESeq2_results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)



## LOAD PREPARED DATA

count_mat_filt <- readRDS(
  file.path(
    prep_output_dir,
    "count_mat_filt.rds"
  )
)

coldata <- readRDS(
  file.path(
    prep_output_dir,
    "coldata.rds"
  )
)


## CHECK PREPARED DATA

count_mat_filt <- as.matrix(count_mat_filt)
coldata <- as.data.frame(coldata)

## check species column
if (!"species" %in% colnames(coldata)) {
  stop("coldata must contain a column named 'species'.")
}

## check sample names
if (!setequal(colnames(count_mat_filt), rownames(coldata))) {
  stop(
    paste(
      "Sample names do not match between",
      "count_mat_filt and coldata."
    )
  )
}

## reorder metadata to match count matrix
coldata <- coldata[
  colnames(count_mat_filt),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(count_mat_filt),
    rownames(coldata)
  )
)

## check counts
if (anyNA(count_mat_filt)) {
  stop("count_mat_filt contains missing values.")
}

if (any(count_mat_filt < 0)) {
  stop("count_mat_filt contains negative values.")
}

if (any(count_mat_filt != round(count_mat_filt))) {
  stop("count_mat_filt must contain raw integer counts.")
}

storage.mode(count_mat_filt) <- "integer"


## SET SPECIES FACTOR AND REFERENCE

## check expected species labels
if (anyNA(coldata$species)) {
  stop(
    paste0(
      "Unexpected species labels were found. Expected: ",
      paste(species_levels, collapse = ", ")
    )
  )
}

## initial reference species
coldata$species <- relevel(
  coldata$species,
  ref = "Ab"
)


## CREATE DESEQ2 DATASET
dds <- DESeqDataSetFromMatrix(
  countData = count_mat_filt,
  colData = coldata,
  design = ~ species
)


## NORMALISATION CHECKS
dds <- estimateSizeFactors(dds)

norm_counts <- counts(
  dds,
  normalized = TRUE
)

## inspect size factors
sizeFactors(dds)

## inspect normalized library totals
colSums(norm_counts)

## save normalized counts
write.csv(
  as.data.frame(norm_counts),
  file.path(
    output_dir,
    "DESeq2_normalized_counts.csv"
  ),
  row.names = TRUE
)


##variance stabilising transform for QC
vsd <- vst(
  dds,
  blind = TRUE
)


##PCA
pcaData <- plotPCA(
  vsd,
  intgroup = "species",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar"),
  digits = 1
)

p_pca <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    color = species
  )
) +
  geom_point(size = 3) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +
  labs(
    title = "PCA of variance-stabilized counts",
    color = "Species"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

p_pca

ggsave(
  filename = file.path(
    output_dir,
    "PCA_species.png"
  ),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)


##SAMPLE DISTANCE HEATMAP
sampleDists <- dist(
  t(assay(vsd))
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

sample_annot <- as.data.frame(
  colData(vsd)[
    ,
    "species",
    drop = FALSE
  ]
)

colnames(sample_annot) <- "Species"

png(
  filename = file.path(
    output_dir,
    "Sample_distance_heatmap.png"
  ),
  width = 2400,
  height = 2200,
  res = 300
)

pheatmap(
  sampleDistMatrix,
  annotation_col = sample_annot,
  annotation_row = sample_annot,
  main = "Sample distances"
)

dev.off()



##FIT DESEQ2 MODELS

fit_reference_models <- function(
    dds,
    species_levels
) {
  
  fitted_models <- vector(
    mode = "list",
    length = length(species_levels)
  )
  
  names(fitted_models) <- species_levels
  
  for (reference_species in species_levels) {
    
    message(
      "Fitting DESeq2 model with ",
      reference_species,
      " as reference..."
    )
    
    dds_tmp <- dds
    
    dds_tmp$species <- relevel(
      factor(
        dds_tmp$species,
        levels = species_levels
      ),
      ref = reference_species
    )
    
    ## run DESeq using the Wald test
    fitted_models[[reference_species]] <- DESeq(
      dds_tmp,
      test = "Wald",
      quiet = TRUE
    )
  }
  
  return(fitted_models)
}


dds_by_reference <- fit_reference_models(
  dds = dds,
  species_levels = species_levels
)


## inspect coefficient names
lapply(
  dds_by_reference,
  resultsNames
)


## MAKE PAIRWISE RESULTS TABLES

make_pairwise_res <- function(
    dds_by_reference,
    species_a,
    species_b,
    alpha = 0.05
) {
  
  if (species_a == species_b) {
    stop("species_a and species_b must be different.")
  }
  
  if (!species_b %in% names(dds_by_reference)) {
    stop(
      paste0(
        "No model was found with ",
        species_b,
        " as the reference."
      )
    )
  }
  
  dds_tmp <- dds_by_reference[[species_b]]
  
  coef_name <- paste0(
    "species_",
    species_a,
    "_vs_",
    species_b
  )
  
  if (!coef_name %in% resultsNames(dds_tmp)) {
    stop(
      paste0(
        "Coefficient not found: ",
        coef_name,
        "\nAvailable coefficients: ",
        paste(
          resultsNames(dds_tmp),
          collapse = ", "
        )
      )
    )
  }
  
  ## obtain Wald-test results
  res_raw <- results(
    dds_tmp,
    name = coef_name,
    alpha = alpha
  )
  
  ## shrink log2 fold changes using apeglm
  ## p-values and padj remain based on the Wald test
  res_shrunk <- lfcShrink(
    dds = dds_tmp,
    coef = coef_name,
    res = res_raw,
    type = "apeglm"
  )
  
  res_tbl <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene_name") %>%
    as_tibble() %>%
    mutate(
      comparison = paste0(
        species_a,
        "_vs_",
        species_b
      ),
      
      ## significance is defined using FDR < 0.05 only
      sig = case_when(
        !is.na(padj) &
          padj < alpha &
          log2FoldChange > 0 ~
          paste(species_a, "higher"),
        
        !is.na(padj) &
          padj < alpha &
          log2FoldChange < 0 ~
          paste(species_b, "higher"),
        
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(
      is.na(padj),
      padj,
      desc(abs(log2FoldChange))
    )
  
  ## retain significant DE genes
  deg_tbl <- res_tbl %>%
    filter(
      !is.na(padj),
      padj < alpha
    )
  
  return(
    list(
      coef_name = coef_name,
      results = res_tbl,
      degs = deg_tbl
    )
  )
}


## VOLCANO PLOTS - editted for poster
make_volcano_plot <- function(
    res_tbl,
    species_a,
    species_b,
    alpha = 0.05,
    lfc_guide = 1,
    labels_per_side = 5
) {
  
  sig_levels <- c(
    paste(species_a, "higher"),
    paste(species_b, "higher"),
    "Not significant"
  )
  
  volcano_df <- res_tbl %>%
    mutate(
      ## prevent infinite values when adjusted p-value equals zero
      neg_log10_padj = case_when(
        is.na(padj) ~ NA_real_,
        
        TRUE ~ -log10(
          pmax(
            padj,
            .Machine$double.xmin
          )
        )
      ),
      
      sig = factor(
        sig,
        levels = sig_levels
      )
    )
  
  ## choose the most significant genes in each direction
  label_genes <- volcano_df %>%
    filter(
      !is.na(padj),
      padj < alpha,
      sig != "Not significant"
    ) %>%
    group_by(sig) %>%
    arrange(
      padj,
      desc(abs(log2FoldChange)),
      .by_group = TRUE
    ) %>%
    slice_head(
      n = labels_per_side
    ) %>%
    ungroup()
  
  color_values <- setNames(
    c(
      "#9ecae1",
      "#fdae6b",
      "grey80"
    ),
    sig_levels
  )
  
  p <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj
    )
  ) +
    
    ## plot non-significant genes first
    geom_point(
      data = filter(
        volcano_df,
        sig == "Not significant"
      ),
      aes(color = sig),
      alpha = 0.65,
      size = 2
    ) +
    
    ## plot significant genes on top
    geom_point(
      data = filter(
        volcano_df,
        sig != "Not significant"
      ),
      aes(color = sig),
      alpha = 0.9,
      size = 2.5
    ) +
    
    ## label top genes in each direction
    geom_text_repel(
      data = label_genes,
      aes(label = gene_name),
      size = 5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "black",
      min.segment.length = 0,
      seed = 1
    ) +
    
    ## FDR significance threshold
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    ## visual effect-size guides only
    geom_vline(
      xintercept = c(
        -lfc_guide,
        lfc_guide
      ),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    scale_color_manual(
      values = color_values,
      breaks = sig_levels,
      drop = FALSE,
      name = "Expression change"
    ) +
    
    scale_y_continuous(
      expand = expansion(
        mult = c(0.02, 0.12)
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_minimal(
      base_size = 16
    ) +
    
    labs(
      title = paste(
        "Volcano Plot:",
        species_a,
        "vs",
        species_b
      ),
      
      subtitle = paste0(
        "Significance: FDR < ",
        alpha,
        "; vertical lines show |log2FC| = ",
        lfc_guide
      ),
      
      x = paste0(
        "Shrunken log2 fold change (",
        species_a,
        " / ",
        species_b,
        ")"
      ),
      
      y = expression(
        -log[10]("adjusted p-value")
      )
    ) +
    
    theme(
      plot.title = element_text(
        size = 20,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        size = 13
      ),
      
      axis.title = element_text(
        size = 16
      ),
      
      axis.text = element_text(
        size = 14
      ),
      
      legend.title = element_text(
        size = 14
      ),
      
      legend.text = element_text(
        size = 13
      ),
      
      legend.key.size = unit(
        0.6,
        "cm"
      ),
      
      plot.margin = margin(
        t = 10,
        r = 25,
        b = 10,
        l = 10
      )
    )
  
  return(
    list(
      plot = p,
      volcano_df = volcano_df,
      label_genes = label_genes
    )
  )
}


## LOOP THROUGH ALL PAIRWISE COMPARISONS
comparisons <- combn(
  species_levels,
  2,
  simplify = FALSE
)

all_results_list <- list()
all_deg_list <- list()
all_volcano_list <- list()
all_comparison_summary_list <- list()


for (comp in comparisons) {
  
  species_a <- comp[1]
  species_b <- comp[2]
  
  comparison_name <- paste0(
    species_a,
    "_vs_",
    species_b
  )
  
  message(
    "Processing comparison: ",
    comparison_name
  )
  
  ## make pairwise results
  pair_res <- make_pairwise_res(
    dds_by_reference = dds_by_reference,
    species_a = species_a,
    species_b = species_b,
    alpha = alpha
  )
  
  res_tbl <- pair_res$results
  deg_tbl <- pair_res$degs
  
  ## store all results
  all_results_list[[comparison_name]] <- res_tbl
  
  ## clean DEG table
  deg_tbl_clean <- deg_tbl %>%
    select(
      gene_name,
      baseMean,
      log2FoldChange,
      lfcSE,
      pvalue,
      padj,
      comparison,
      sig
    )
  
  all_deg_list[[comparison_name]] <- deg_tbl_clean
  
  ## summarise DE genes
  comparison_summary <- tibble(
    comparison = comparison_name,
    
    species_a = species_a,
    
    species_b = species_b,
    
    n_tested = sum(
      !is.na(res_tbl$pvalue)
    ),
    
    n_DEGs = nrow(deg_tbl_clean),
    
    n_species_a_higher = sum(
      deg_tbl_clean$sig ==
        paste(species_a, "higher")
    ),
    
    n_species_b_higher = sum(
      deg_tbl_clean$sig ==
        paste(species_b, "higher")
    )
  )
  
  all_comparison_summary_list[[comparison_name]] <-
    comparison_summary
  
  ## save all results
  write.csv(
    res_tbl,
    file.path(
      output_dir,
      paste0(
        "DESeq2_",
        comparison_name,
        "_all_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  ## save significant DE genes
  write.csv(
    deg_tbl_clean,
    file.path(
      output_dir,
      paste0(
        "DESeq2_",
        comparison_name,
        "_DEGs_FDR0.05.csv"
      )
    ),
    row.names = FALSE
  )
  
  ## make volcano plot
  volcano_out <- make_volcano_plot(
    res_tbl = res_tbl,
    species_a = species_a,
    species_b = species_b,
    alpha = alpha,
    lfc_guide = volcano_lfc_guide,
    labels_per_side = labels_per_side
  )
  
  p <- volcano_out$plot
  
  all_volcano_list[[comparison_name]] <- p
  
  ## save volcano plot
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Volcano_",
        comparison_name,
        "_clean.png"
      )
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


##COMBINE AND EXPORT ALL RESULTS

all_pairwise_results <- bind_rows(
  all_results_list
)

write.csv(
  all_pairwise_results,
  file.path(
    output_dir,
    "DESeq2_all_pairwise_results.csv"
  ),
  row.names = FALSE
)


## combine all significant DE genes
all_pairwise_DEGs <- bind_rows(
  all_deg_list
)

write.csv(
  all_pairwise_DEGs,
  file.path(
    output_dir,
    "DESeq2_all_pairwise_DEGs_FDR0.05.csv"
  ),
  row.names = FALSE
)


## combine pairwise summaries
pairwise_DEG_summary <- bind_rows(
  all_comparison_summary_list
)

write.csv(
  pairwise_DEG_summary,
  file.path(
    output_dir,
    "DESeq2_pairwise_DEG_summary.csv"
  ),
  row.names = FALSE
)

pairwise_DEG_summary

