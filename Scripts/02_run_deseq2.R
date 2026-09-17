library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(ggrepel)
library(apeglm)
library(grid)
library(ashr)
library(patchwork)
library(readr)


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
prep_output_dir <- "Data/Processed"

## output folder
output_dir <- "Results/DESeq2"

## figure folder
figure_dir <- "Figures/DESeq2"


dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figure_dir,
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

## CHECK PROTEIN-CODING FILTER

cat(
  "Number of genes:",
  nrow(count_mat_filt),
  "\n"
)

cat(
  "Number of samples:",
  ncol(count_mat_filt),
  "\n"
)

noncoding_genes <- c(
  "5S_rRNA",
  "Metazoa_SRP",
  "U6"
)

if (any(noncoding_genes %in% rownames(count_mat_filt))) {
  stop(
    "Non-protein-coding genes are still present."
  )
}

if (nrow(count_mat_filt) != 8522) {
  warning(
    paste(
      "Expected 8,522 protein-coding genes, but found",
      nrow(count_mat_filt)
    )
  )
}

message(
  "Protein-coding matrix check passed: ",
  nrow(count_mat_filt),
  " genes x ",
  ncol(count_mat_filt),
  " samples."
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


## FIGURE 1: PCA + SAMPLE-DISTANCE HEATMAP

## consistent species colours used across both panels
species_colours <- c(
  "Ab" = "#0072B2",
  "Mz" = "#E69F00",
  "Nb" = "#009E73",
  "On" = "#CC79A7",
  "Pn" = "#D55E00"
)


## PCA

pcaData <- plotPCA(
  vsd,
  intgroup = "species",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar"),
  digits = 1
)

pcaData$species <- factor(
  pcaData$species,
  levels = species_levels
)

p_pca <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    colour = species
  )
) +
  geom_point(
    size = 3.2,
    alpha = 0.9
  ) +
  scale_colour_manual(
    values = species_colours,
    breaks = species_levels,
    name = "Species"
  ) +
  labs(
    x = paste0(
      "PC1 (",
      percentVar[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percentVar[2],
      "%)"
    )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text = element_text(
      size = 10
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      10, 10, 10, 10
    )
  ) +
  guides(
    colour = guide_legend(
      title = "Species",
      nrow = 1
    )
  )


## Extract shared species legend

legend_plot <- p_pca +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    colour = guide_legend(
      title = "Species",
      nrow = 1
    )
  )

shared_legend <- cowplot::get_legend(legend_plot)

p_pca_nolegend <- p_pca +
  theme(
    legend.position = "none"
  )



## Sample-distance heatmap

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

sample_annot$Species <- factor(
  sample_annot$Species,
  levels = species_levels
)

annotation_colours <- list(
  Species = species_colours
)

p_heatmap <- pheatmap(
  sampleDistMatrix,
  annotation_col = sample_annot,
  annotation_row = sample_annot,
  treeheight_row = 70,
  treeheight_col = 70,
  cellheight = 24,
  annotation_colors = annotation_colours,
  annotation_legend = FALSE,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  silent = TRUE
)


## Main A + B panels with tags
main_panels <- (
  p_pca_nolegend |
    patchwork::wrap_elements(
      full = p_heatmap$gtable
    )
) +
  patchwork::plot_layout(
    widths = c(0.85, 1.15)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )


## figure with spacer + centred shared legend

figure1 <- (
  patchwork::wrap_elements(
    full = main_panels
  ) /
    patchwork::plot_spacer() /
    patchwork::wrap_elements(
      full = shared_legend
    )
) +
  patchwork::plot_layout(
    heights = c(1, 0.025, 0.07)
  )

figure1


## Save

ggsave(
  filename = file.path(
    figure_dir,
    "Figure1_PCA_sample_distance.png"
  ),
  plot = figure1,
  width = 12,
  height = 7,
  dpi = 600
)



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


## VOLCANO PLOTS - edited for poster
make_volcano_plot <- function(
    res_tbl,
    species_a,
    species_b,
    alpha = 0.05,
    lfc_guide = 1,
    labels_per_side = 5
) {
  
  sig_levels <- c(
    "Focal species higher",
    "Other species higher",
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
        "Log2 fold change (",
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
      figure_dir,
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




## SPECIES-VS-REST DIFFERENTIAL EXPRESSION

species_vs_rest_output_dir <- file.path(
  output_dir,
  "species_vs_rest"
)

species_vs_rest_figure_dir <- file.path(
  figure_dir,
  "species_vs_rest"
)

dir.create(
  species_vs_rest_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  species_vs_rest_figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


## use the Ab-reference fitted model
dds_species <- dds_by_reference[["Ab"]]

resultsNames(dds_species)


## construct one species-vs-rest numeric contrast
make_species_vs_rest_contrast <- function(
    dds,
    focal_species,
    species_levels
) {
  
  coef_names <- resultsNames(dds)
  
  contrast_vector <- rep(
    0,
    length(coef_names)
  )
  
  names(contrast_vector) <- coef_names
  
  other_species <- setdiff(
    species_levels,
    focal_species
  )
  
  ## Ab is the reference species in this model therefore its contribution is 
  ## represented by the intercept, which cancels from the contrast
  
  if (focal_species != "Ab") {
    
    focal_coef <- paste0(
      "species_",
      focal_species,
      "_vs_Ab"
    )
    
    contrast_vector[focal_coef] <- 1
  }
  
  for (other_species_name in other_species) {
    
    if (other_species_name != "Ab") {
      
      other_coef <- paste0(
        "species_",
        other_species_name,
        "_vs_Ab"
      )
      
      contrast_vector[other_coef] <-
        contrast_vector[other_coef] -
        1 / length(other_species)
    }
  }
  
  contrast_vector
}


## make species-vs-rest result table
make_species_vs_rest_res <- function(
    dds,
    focal_species,
    species_levels,
    alpha = 0.05
) {
  
  contrast_vector <- make_species_vs_rest_contrast(
    dds = dds,
    focal_species = focal_species,
    species_levels = species_levels
  )
  
  res_raw <- results(
    dds,
    contrast = contrast_vector,
    alpha = alpha
  )
  
  ## shrink composite contrast using ashr
  res_shrunk <- lfcShrink(
    dds = dds,
    contrast = contrast_vector,
    res = res_raw,
    type = "ashr"
  )
  
  focal_label <- paste(
    focal_species,
    "higher"
  )
  
  rest_label <- paste(
    "Other species higher"
  )
  
  res_tbl <- as.data.frame(
    res_shrunk
  ) %>%
    tibble::rownames_to_column(
      "gene_name"
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      comparison = paste0(
        focal_species,
        "_vs_rest"
      ),
      
      sig = dplyr::case_when(
        !is.na(padj) &
          padj < alpha &
          log2FoldChange > 0 ~
          focal_label,
        
        !is.na(padj) &
          padj < alpha &
          log2FoldChange < 0 ~
          rest_label,
        
        TRUE ~
          "Not significant"
      )
    ) %>%
    dplyr::arrange(
      is.na(padj),
      padj,
      dplyr::desc(
        abs(log2FoldChange)
      )
    )
  
  deg_tbl <- res_tbl %>%
    dplyr::filter(
      !is.na(padj),
      padj < alpha
    )
  
  list(
    contrast = contrast_vector,
    results = res_tbl,
    degs = deg_tbl
  )
}



## SPECIES VS REST VOLCANO PLOT

## load final curated visual system candidate gene list

visual_genes_plot <- readr::read_csv(
  "Results/targeted_visual_genes/curated_visual_genes_all_species_vs_rest.csv",
  show_col_types = FALSE
) %>%
  dplyr::pull(gene_name.y) %>%
  trimws() %>%
  tolower() %>%
  unique()


## check visual-system candidate representation in DESeq2 dataset
visual_genes_found <- intersect(
  visual_genes_plot,
  tolower(
    rownames(dds)
  )
)

visual_genes_missing <- setdiff(
  visual_genes_plot,
  tolower(
    rownames(dds)
  )
)

cat(
  "Visual candidate genes represented in DESeq2 dataset:",
  length(visual_genes_found),
  "\n"
)

cat(
  "Visual candidate genes not represented by these symbols:",
  length(visual_genes_missing),
  "\n"
)

print(
  visual_genes_found
)

print(
  visual_genes_missing
)



## SPECIES-VS-REST VOLCANO PLOT

make_species_vs_rest_volcano <- function(
    res_tbl,
    focal_species,
    alpha = 0.05,
    lfc_guide = 1,
    labels_per_side = 2,
    visual_genes = visual_genes_plot,
    x_limit = 20,
    y_limit = 55
) {
  
  plot_tbl <- res_tbl %>%
    dplyr::mutate(
      
      ## display gene symbol
      gene_label = tolower(
        trimws(
          gene_name
        )
      ),
      
      ## protect against infinite values
      neg_log10_padj = dplyr::case_when(
        is.na(padj) ~ NA_real_,
        
        TRUE ~ -log10(
          pmax(
            padj,
            .Machine$double.xmin
          )
        )
      ),
      
      ## display coordinates for common volcano plot axes
      plot_log2FC = pmax(
        pmin(log2FoldChange, x_limit),
        -x_limit
      ),
      
      plot_neg_log10_padj = pmin(
        neg_log10_padj,
        y_limit
      ),
      
      ## expression direction for plotting
      expression_change = dplyr::case_when(
        
        !is.na(padj) &
          padj < alpha &
          log2FoldChange > 0 ~
          "Focal species higher",
        
        !is.na(padj) &
          padj < alpha &
          log2FoldChange < 0 ~
          "Other species higher",
        
        TRUE ~
          "Not significant"
      ),
      
      ## final curated visual-system candidate
      is_visual_gene =
        gene_label %in% visual_genes &
        !is.na(padj) &
        padj < alpha
      )
  
  
  ## strongest genes with higher expression in focal species
  
  top_positive <- plot_tbl %>%
    dplyr::filter(
      expression_change == "Focal species higher",
      !is_visual_gene
    ) %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(log2FoldChange))
    ) %>%
    dplyr::slice_head(
      n = labels_per_side
    )
  
  
  ## strongest genes with higher expression in other species
  
  top_negative <- plot_tbl %>%
    dplyr::filter(
      expression_change == "Other species higher",
      !is_visual_gene
    ) %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(log2FoldChange))
    ) %>%
    dplyr::slice_head(
      n = labels_per_side
    )
  
  ## all significant final curated visual-system genes
  visual_points <- plot_tbl %>%
    dplyr::filter(
      is_visual_gene
    )
  
  ## label only the strongest visual-system candidates
  visual_labels <- visual_points %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(log2FoldChange))
    ) %>%
    dplyr::slice_head(n = 5)
  
  ## combine strongest general DE genes
  
  top_labels <- dplyr::bind_rows(
    top_positive,
    top_negative
  ) %>%
    dplyr::distinct(
      gene_name,
      .keep_all = TRUE
    ) %>%
    dplyr::filter(
      abs(log2FoldChange) <= x_limit,
      neg_log10_padj <= y_limit
    )
  
  ## move sowahd separately to avoid acer2 label
  sowahd_label <- top_labels %>%
    dplyr::filter(
      gene_label == "sowahd"
    )
  
  top_labels <- top_labels %>%
    dplyr::filter(
      gene_label != "sowahd"
    )
  
  ## colour levels
  
  sig_levels <- c(
    "Focal species higher",
    "Other species higher",
    "Not significant"
  )
  
  
  color_values <- setNames(
    c(
      "#9ecae1",
      "#fdae6b",
      "grey80"
    ),
    sig_levels
  )
  
  ## genes truncated by common display limits

  truncated_points <- plot_tbl %>%
    dplyr::filter(
      abs(log2FoldChange) > x_limit |
        neg_log10_padj > y_limit
    ) %>%
    dplyr::mutate(
      trunc_label = dplyr::case_when(
        
        neg_log10_padj > y_limit ~
          paste0(
            gene_label,
            " (",
            round(neg_log10_padj, 1),
            ")"
          ),
        
        abs(log2FoldChange) > x_limit ~
          paste0(
            gene_label,
            " (",
            round(log2FoldChange, 1),
            ")"
          ),
        
        TRUE ~ gene_label
      )
    )
  
  ## create plot
  
  p <- ggplot(
    plot_tbl,
    aes(
      x = plot_log2FC,
      y = plot_neg_log10_padj
    )
  ) +
    
    ## non-significant genes
    geom_point(
      data = dplyr::filter(
        plot_tbl,
        expression_change ==
          "Not significant"
      ),
      aes(
        color = expression_change
      ),
      alpha = 0.6,
      size = 1.7
    ) +
    
    ## significant genes
    geom_point(
      data = dplyr::filter(
        plot_tbl,
        expression_change !=
          "Not significant"
      ),
      aes(
        color = expression_change
      ),
      alpha = 0.9,
      size = 2.2
    ) +
    
    ## highlight significant visual-system genes
    geom_point(
      data = visual_points,
      aes(
        x = plot_log2FC,
        y = plot_neg_log10_padj,
        shape = "Visual-system candidate"
      ),
      inherit.aes = FALSE,
      size = 4,
      stroke = 1,
      fill = "white",
      color = "black"
    ) +
    
    scale_shape_manual(
      name = NULL,
      values = c(
        "Visual-system candidate" = 24
      )
    ) +
    
    ##legends easier to read
    guides(
      color = guide_legend(
        order = 1,
        title = "Expression change",
        override.aes = list(
          size = 3,
          alpha = 1
        )
      ),
      shape = guide_legend(
        order = 2,
        title = NULL,
        override.aes = list(
          size = 4,
          fill = "white",
          color = "black"
        )
      )
    ) +
    
    ## label strongest general DE genes
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(
        x = plot_log2FC,
        y = plot_neg_log10_padj,
        label = gene_label
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0,
      seed = 1
    ) +
    
    ## separately position sowahd
    ggrepel::geom_text_repel(
      data = sowahd_label,
      aes(
        x = plot_log2FC,
        y = plot_neg_log10_padj,
        label = gene_label
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.5,
      nudge_x = 1.5,
      nudge_y = -3,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 999
    ) +
    
    ## label significant visual-system genes
    ggrepel::geom_text_repel(
      data = visual_labels,
      aes(
        x = plot_log2FC,
        y = plot_neg_log10_padj,
        label = gene_label
      ),
      inherit.aes = FALSE,
      color = "black",
      fontface = "bold",
      size = 3.2,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 123
    ) +
    
    ## label genes truncated at x-axis limits
    ggrepel::geom_text_repel(
      data = truncated_points %>%
        dplyr::filter(
          abs(log2FoldChange) > x_limit
        ),
      aes(
        x = plot_log2FC,
        y = plot_neg_log10_padj,
        label = trunc_label
      ),
      inherit.aes = FALSE,
      color = "black",
      fontface = "italic",
      size = 3.2,
      nudge_x = 3,
      direction = "y",
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 456
    ) +
    
    ## label genes truncated at y-axis limit
    ggrepel::geom_text_repel(
      data = truncated_points %>%
        dplyr::filter(
          neg_log10_padj > y_limit,
          abs(log2FoldChange) <= x_limit
        ),
      aes(
        x = plot_log2FC,
        y = plot_neg_log10_padj,
        label = trunc_label
      ),
      inherit.aes = FALSE,
      color = "black",
      fontface = "italic",
      size = 3.2,
      nudge_y = -4,
      direction = "x",
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 789
    ) +
    
    ## FDR significance threshold
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    ## effect-size guides only
    geom_vline(
      xintercept = c(
        -lfc_guide,
        lfc_guide
      ),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    scale_color_manual(
      values = color_values,
      breaks = sig_levels,
      drop = FALSE,
      name = "Expression change"
    ) +
    
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = seq(-x_limit, x_limit, 10)
    ) +
    
    scale_y_continuous(
      limits = c(0, y_limit),
      breaks = seq(0, 50, 10),
      expand = expansion(
        mult = c(0, 0.02)
      )
    ) +
    
    labs(
      x = expression(
        log[2] * " fold change"
      ),
      y = expression(
        -log[10] * "(FDR)"
      )
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    theme(
      axis.title = element_text(
        size = 11
      ),
      
      axis.text = element_text(
        size = 9
      ),
      
      legend.position = "bottom",
      
      legend.box = "horizontal",
      
      legend.title = element_text(
        face = "bold",
        size = 10
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      plot.margin = margin(
        8, 8, 8, 8
      )
    )
  
  
  return(p)
}



## RUN ALL SPECIES VS REST COMPARISONS

species_vs_rest_results <- list()

species_vs_rest_degs <- list()

species_vs_rest_plots <- list()

species_vs_rest_summary <- list()


for (focal_species in species_levels) {
  
  message(
    "Running ",
    focal_species,
    " vs remaining species..."
  )
  
  species_res <- make_species_vs_rest_res(
    dds = dds_species,
    focal_species = focal_species,
    species_levels = species_levels,
    alpha = alpha
  )
  
  res_tbl <- species_res$results
  
  deg_tbl <- species_res$degs
  
  species_vs_rest_results[[focal_species]] <- res_tbl
  
  species_vs_rest_degs[[focal_species]] <- deg_tbl
  
  
  ## export complete results
  
  write.csv(
    res_tbl,
    file.path(
      species_vs_rest_output_dir,
      paste0(
        "DESeq2_",
        focal_species,
        "_vs_rest_all_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  ## export significant genes
  
  write.csv(
    deg_tbl,
    file.path(
      species_vs_rest_output_dir,
      paste0(
        "DESeq2_",
        focal_species,
        "_vs_rest_DEGs_FDR0.05.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  ## summary
  
  species_vs_rest_summary[[focal_species]] <- tibble::tibble(
    species = focal_species,
    
    n_tested = sum(
      !is.na(
        res_tbl$pvalue
      )
    ),
    
    n_DEGs = nrow(
      deg_tbl
    ),
    
    n_species_higher = sum(
      deg_tbl$sig ==
        paste(
          focal_species,
          "higher"
        )
    ),
    
    n_other_species_higher = sum(
      deg_tbl$sig ==
        "Other species higher"
    )
  )
  
  
  ## volcano plot
  volcano_plot <- make_species_vs_rest_volcano(
    res_tbl = res_tbl,
    focal_species = focal_species,
    alpha = alpha,
    lfc_guide = volcano_lfc_guide,
    labels_per_side = 3,
    visual_genes = visual_genes_plot
  )
  
  species_vs_rest_plots[[focal_species]] <- volcano_plot
  
  ggsave(
    filename = file.path(
      species_vs_rest_figure_dir,
      paste0(
        "Volcano_",
        focal_species,
        "_vs_rest.png"
      )
    ),
    plot = volcano_plot,
    width = 7,
    height = 6,
    dpi = 300
  )
}



## CHECK
## Investigate annotation, orthology and filtering?
## Examine final curated visual system candidates across all species vs rest DE comparisons

## VISUAL SYSTEM CANDIDATE RESULTS ACROSS SPECIES VS REST COMPARISONS

visual_candidate_results <- dplyr::bind_rows(
  lapply(
    names(
      species_vs_rest_results
    ),
    function(sp) {
      
      x <- species_vs_rest_results[[sp]]
      
      x %>%
        dplyr::mutate(
          species = sp,
          
          gene_match = tolower(
            trimws(
              gene_name
            )
          )
        ) %>%
        
        dplyr::filter(
          gene_match %in%
            visual_genes_plot
        ) %>%
        
        dplyr::select(
          species,
          gene_name,
          log2FoldChange,
          lfcSE,
          pvalue,
          padj
        )
    }
  )
)


## summary of significant visual-system candidates
## significance follows the main DESeq2 criterion:
## FDR-adjusted p-value < 0.05

visual_candidate_summary <- visual_candidate_results %>%
  dplyr::mutate(
    
    significant =
      !is.na(padj) &
      padj < alpha,
    
    direction = dplyr::case_when(
      
      significant &
        log2FoldChange > 0 ~
        paste(
          species,
          "higher"
        ),
      
      significant &
        log2FoldChange < 0 ~
        "Other species higher",
      
      TRUE ~
        "Not significant"
    )
  ) %>%
  
  dplyr::group_by(
    species
  ) %>%
  
  dplyr::summarise(
    
    visual_candidates_tested =
      dplyr::n(),
    
    significant_visual_candidates =
      sum(
        significant,
        na.rm = TRUE
      ),
    
    significant_genes = paste(
      gene_name[
        significant
      ],
      collapse = "; "
    ),
    
    .groups = "drop"
  )


## display

print(
  visual_candidate_results,
  n = Inf,
  width = Inf
)

print(
  visual_candidate_summary,
  n = Inf,
  width = Inf
)


## export full candidate results

write.csv(
  visual_candidate_results,
  file.path(
    species_vs_rest_output_dir,
    "DESeq2_species_vs_rest_visual_candidate_results.csv"
  ),
  row.names = FALSE
)


## export candidate summary

write.csv(
  visual_candidate_summary,
  file.path(
    species_vs_rest_output_dir,
    "DESeq2_species_vs_rest_visual_candidate_summary.csv"
  ),
  row.names = FALSE
)



## FIVE-PANEL SPECIES VS REST VOLCANO FIGURE

library(patchwork)


## order species in final figure
species_plot_order <- c(
  "Ab",
  "Mz",
  "Nb",
  "On",
  "Pn"
)


## retrieve plots in required order
final_volcano_plots <- species_vs_rest_plots[
  species_plot_order
]


## combine plots
## top row    = Ab, Mz
## middle row = Nb, On
## bottom row = Pn

volcano_multi_panel <- patchwork::wrap_plots(
  plotlist = final_volcano_plots,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )

volcano_multi_panel

ggsave(
  filename = file.path(
    species_vs_rest_figure_dir,
    "Figure2_species_vs_rest_volcano.png"
  ),
  plot = volcano_multi_panel,
  width = 8.2,
  height = 11.2,
  dpi = 600
)


