library(ggplot2)
library(stringr)
library(DESeq2)
library(dplyr)
library(tibble)
library(pheatmap)
library(grid)
library(tidyr)
library(purrr)


go_plot_dir <- "Figures/GO_plots"

dir.create(
  go_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

plot_go_terms <- function(
    go_tbl,
    comparison_name,
    group_name,
    top_n = 15
) {
  
  plot_data <- go_tbl %>%
    filter(
      comparison == comparison_name,
      group == group_name,
      !is.na(p_value),
      p_value < 0.05
    ) %>%
    arrange(p_value) %>%
    slice_head(n = top_n) %>%
    mutate(
      term_name = factor(
        term_name,
        levels = rev(unique(term_name))
      ),
      neg_log10_p = -log10(
        pmax(
          p_value,
          .Machine$double.xmin
        )
      )
    )
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  ggplot(
    plot_data,
    aes(
      x = neg_log10_p,
      y = term_name,
      size = intersection_size
    )
  ) +
    geom_point(
      aes(color = source),
      alpha = 0.85
    ) +
    scale_size_continuous(
      name = "Genes"
    ) +
    labs(
      title = paste(
        comparison_name,
        "-",
        group_name
      ),
      x = expression(
        -log[10]("FDR-adjusted p-value")
      ),
      y = "GO term",
      color = "GO source"
    ) +
    theme_classic(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 10
      )
    )
}

go_plot_groups <- all_GO_export %>%
  distinct(
    comparison,
    group
  )

go_plots <- pmap(
  go_plot_groups,
  function(comparison, group) {
    
    plot_object <- plot_go_terms(
      go_tbl = all_GO_export,
      comparison_name = comparison,
      group_name = group,
      top_n = 15
    )
    
    if (!is.null(plot_object)) {
      
      file_name <- paste0(
        "GO_dotplot_",
        comparison,
        "_",
        gsub(" ", "_", group),
        ".png"
      )
      
      ggsave(
        filename = file.path(
          go_plot_dir,
          file_name
        ),
        plot = plot_object,
        width = 10,
        height = 7,
        dpi = 300
      )
    }
    
    plot_object
  }
)

names(go_plots) <- paste(
  go_plot_groups$comparison,
  go_plot_groups$group,
  sep = "_"
)


vision_go_plot_data <- vision_GO_export %>%
  filter(
    !is.na(p_value),
    p_value < 0.05
  ) %>%
  group_by(
    comparison,
    group
  ) %>%
  arrange(
    p_value,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 15
  ) %>%
  ungroup() %>%
  mutate(
    comparison_group = paste(
      comparison,
      group,
      sep = ": "
    ),
    neg_log10_p = -log10(
      pmax(
        p_value,
        .Machine$double.xmin
      )
    )
  )

if (nrow(vision_go_plot_data) > 0) {
  
  p_vision_go <- ggplot(
    vision_go_plot_data,
    aes(
      x = neg_log10_p,
      y = reorder(
        term_name,
        neg_log10_p
      ),
      size = intersection_size
    )
  ) +
    geom_point(
      aes(color = comparison_group),
      alpha = 0.85
    ) +
    facet_wrap(
      ~ comparison_group,
      scales = "free_y"
    ) +
    labs(
      title = "Vision-related GO enrichment",
      x = expression(
        -log[10]("FDR-adjusted p-value")
      ),
      y = "GO term",
      size = "Genes",
      color = "Comparison and group"
    ) +
    theme_classic(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8
      ),
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(
      go_plot_dir,
      "GO_vision_all_pairwise_dotplot.png"
    ),
    plot = p_vision_go,
    width = 14,
    height = 10,
    dpi = 300
  )
}


visual_genes <- toupper(c(
  "RH1",
  "RHO",
  "SWS1",
  "SWS2A",
  "SWS2B",
  "RH2A",
  "RH2AALPHA",
  "RH2ABETA",
  "RH2B",
  "LWS",
  "OPN1SW1",
  "OPN1SW2",
  "OPN1MW",
  "OPN1LW",
  "GNAT1",
  "GNAT2",
  "GNB1",
  "GNGT1",
  "GNGT2",
  "PDE6A",
  "PDE6B",
  "PDE6C",
  "CNGA1",
  "CNGA3",
  "CNGB1",
  "CNGB3",
  "ARR3",
  "SAG",
  "RCVRN",
  "RGS9",
  "GUCA1A",
  "GUCA1B",
  "RPE65",
  "LRAT",
  "ABCA4",
  "RDH5",
  "RDH8",
  "CRX",
  "PAX6",
  "NRL",
  "NR2E3",
  "OTX2",
  "SIX6"
))

visual_DE_genes <- imap_dfr(
  all_deg_list,
  function(deg_tbl, comparison_name) {
    
    deg_tbl %>%
      mutate(
        gene_upper = toupper(
          trimws(gene_name)
        )
      ) %>%
      filter(
        !is.na(gene_upper),
        gene_upper != "",
        gene_upper %in% visual_genes
      ) %>%
      select(
        -gene_upper
      ) %>%
      mutate(
        comparison = comparison_name
      )
  }
)

visual_DE_genes_export <- visual_DE_genes %>%
  mutate(
    gene_name = toupper(trimws(gene_name))
  ) %>%
  select(
    comparison,
    gene_name,
    baseMean,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj,
    sig
  ) %>%
  arrange(
    comparison,
    padj,
    desc(abs(log2FoldChange))
  )

write.csv(
  visual_DE_genes_export,
  "DEGs_opsin_phototransduction_visual_genes.csv",
  row.names = FALSE
)


visual_DE_gene_summary <- visual_DE_genes_export %>%
  group_by(
    comparison,
    sig
  ) %>%
  summarise(
    n_visual_DE_genes = n_distinct(gene_name),
    genes = paste(
      sort(unique(gene_name)),
      collapse = "; "
    ),
    .groups = "drop"
  )


write.csv(
  visual_DE_gene_summary,
  "DEGs_opsin_phototransduction_visual_summary.csv",
  row.names = FALSE
)


visual_DE_gene_frequency <- visual_DE_genes_export %>%
  dplyr::count(
    gene_name,
    sort = TRUE,
    name = "n_comparisons"
  )


write.csv(
  visual_DE_gene_frequency,
  "DEGs_visual_gene_frequency.csv",
  row.names = FALSE
)


visual_DE_genes_export %>%
  select(
    comparison,
    gene_name,
    log2FoldChange,
    padj,
    sig
  )


visual_DE_gene_summary

visual_DE_gene_frequency


## complete comparisons table
all_visual_groups <- all_deg_list %>%
  imap_dfr(
    function(deg_tbl, comparison_name) {
      tibble(
        comparison = comparison_name,
        sig = unique(deg_tbl$sig)
      )
    }
  ) %>%
  filter(
    !is.na(sig),
    sig != "Not significant"
  ) %>%
  distinct()

visual_DE_gene_summary_complete <- all_visual_groups %>%
  left_join(
    visual_DE_gene_summary,
    by = c("comparison", "sig")
  ) %>%
  mutate(
    n_visual_DE_genes = tidyr::replace_na(
      n_visual_DE_genes,
      0L
    ),
    genes = tidyr::replace_na(
      genes,
      "None"
    )
  ) %>%
  arrange(
    comparison,
    sig
  )


##HEATMAP - candidate visual-gene expression
## producing a species-ordered heatmap

## create output directory
heatmap_dir <- "Figures/"

dir.create(
  heatmap_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## use the significant candidate genes identified across all comparisons
significant_visual_genes <- visual_DE_genes_export %>%
  mutate(
    gene_name = toupper(trimws(gene_name))
  ) %>%
  filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  distinct(gene_name) %>%
  pull(gene_name)

significant_visual_genes
length(significant_visual_genes)

## variance-stabilising transformation
## blind = FALSE retains differences associated with the species design
vsd <- vst(
  dds,
  blind = FALSE
)

vst_mat <- assay(vsd)

## match the DESeq2 matrix rows to gene names
## this assumes rownames(dds) are the shared gene names
vst_gene_names <- toupper(
  trimws(rownames(vst_mat))
)

keep_visual <- vst_gene_names %in% significant_visual_genes

visual_vst_mat <- vst_mat[
  keep_visual,
  ,
  drop = FALSE
]

## standardise row names
rownames(visual_vst_mat) <- vst_gene_names[
  keep_visual
]

## collapse duplicate gene-name rows if any remain
visual_vst_mat <- rowsum(
  visual_vst_mat,
  group = rownames(visual_vst_mat),
  reorder = FALSE
)

## check genes found in the expression matrix
genes_found <- rownames(visual_vst_mat)

genes_missing <- setdiff(
  significant_visual_genes,
  genes_found
)

genes_found
genes_missing

## arrange genes alphabetically
visual_vst_mat <- visual_vst_mat[
  order(rownames(visual_vst_mat)),
  ,
  drop = FALSE
]

##sample annotation
sample_annotation <- as.data.frame(
  colData(dds)
) %>%
  select(species)

sample_annotation$species <- factor(
  sample_annotation$species
)

rownames(sample_annotation) <- colnames(visual_vst_mat)

## confirm sample order matches
stopifnot(
  identical(
    colnames(visual_vst_mat),
    rownames(sample_annotation)
  )
)

##  save unscaled transformed expression values
write.csv(
  visual_vst_mat,
  file.path(
    heatmap_dir,
    "Visual_genes_VST_expression.csv"
  ),
  row.names = TRUE
)

## row scale each gene so that relative expression patterns are comparable
visual_vst_scaled <- t(
  scale(
    t(visual_vst_mat)
  )
)

## remove genes with zero variance, if present
visual_vst_scaled <- visual_vst_scaled[
  apply(
    visual_vst_scaled,
    1,
    function(x) all(is.finite(x))
  ),
  ,
  drop = FALSE
]

##save scaled values
write.csv(
  visual_vst_scaled,
  file.path(
    heatmap_dir,
    "Visual_genes_VST_row_scaled.csv"
  ),
  row.names = TRUE
)



## publication-style species-ordered heatmap
## heatmap palette
heatmap_colours <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(100)

## species annotation colours
species_colours <- c(
  "Ab" = "#E69F00",
  "Mz" = "#CC79A7",
  "Nb" = "#009E73",
  "On" = "#0072B2",
  "Pn" = "#D55E00"
)

annotation_colours <- list(
  species = species_colours
)

## order samples by species
sample_order <- order(
  sample_annotation$species
)

visual_vst_species_order <- visual_vst_scaled[
  ,
  sample_order,
  drop = FALSE
]

sample_annotation_species_order <- sample_annotation[
  sample_order,
  ,
  drop = FALSE
]

## add gaps between species
species_counts <- table(
  sample_annotation_species_order$species
)

gaps_col <- cumsum(
  as.numeric(species_counts)
)

gaps_col <- gaps_col[
  -length(gaps_col)
]

## create heatmap object
p_visual <- pheatmap(
  visual_vst_species_order,
  
  annotation_col = sample_annotation_species_order,
  annotation_colors = annotation_colours,
  
  color = heatmap_colours,
  breaks = seq(
    min(
      visual_vst_species_order,
      na.rm = TRUE
    ),
    max(
      visual_vst_species_order,
      na.rm = TRUE
    ),
    length.out = 101
  ),
  
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  
  gaps_col = gaps_col,
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  fontsize = 13,
  fontsize_row = 16,
  fontsize_col = 10,
  
  angle_col = 90,
  border_color = NA,
  
  annotation_names_col = FALSE,
  annotation_legend = TRUE,
  legend = TRUE,
  
  main = "Visual-gene expression grouped by species",
  
  silent = TRUE
)

## italicise gene names
row_label_index <- which(
  p_visual$gtable$layout$name == "row_names"
)

if (length(row_label_index) == 1) {
  
  p_visual$gtable$grobs[[row_label_index]]$gp <- grid::gpar(
    fontsize = 14,
    fontface = "italic"
  )
}

## save high-resolution PNG
png(
  filename = file.path(
    heatmap_dir,
    "Visual_genes_VST_heatmap_species_order_publication.png"
  ),
  width = 12,
  height = 8,
  units = "in",
  res = 600
)

grid::grid.newpage()
grid::grid.draw(
  p_visual$gtable
)

dev.off()

## save vector PDF
pdf(
  file = file.path(
    heatmap_dir,
    "Visual_genes_VST_heatmap_species_order_publication.pdf"
  ),
  width = 12,
  height = 8
)

grid::grid.newpage()
grid::grid.draw(
  p_visual$gtable
)

dev.off()



## summary table
## summary table
visual_DE_gene_frequency_table <- visual_DE_genes_export %>%
  dplyr::mutate(
    gene_name = toupper(trimws(gene_name))
  ) %>%
  dplyr::distinct(
    comparison,
    gene_name
  ) %>%
  dplyr::count(
    gene_name,
    name = "n_significant_comparisons",
    sort = TRUE
  ) %>%
  dplyr::rename(
    Gene = gene_name,
    `Significant comparisons` = n_significant_comparisons
  )

visual_DE_gene_frequency_table


##save
write.csv(
  visual_DE_gene_frequency_table,
  file.path(
    heatmap_dir,
    "Visual_gene_significant_comparison_frequency.csv"
  ),
  row.names = FALSE
)



## biological interpretation
## what biological processes are actually being enriched?

## 50 most sig GO terms across all comparisons
all_GO_export %>%
  arrange(p_value) %>%
  select(
    comparison,
    group,
    term_name,
    source,
    p_value,
    intersection_size
  ) %>%
  print(n = 50)


## most frequently enriched GO terms
all_GO_export %>%
  dplyr::count(
    term_name,
    sort = TRUE
  ) %>%
  print(n = 50)

## GO domains represented
all_GO_export %>%
  dplyr::count(source)

##save
go_summary <- all_GO_export %>%
  arrange(p_value)

write.csv(
  go_summary,
  "GO_all_terms_ranked.csv",
  row.names = FALSE
)



## table showing candidate visual genes across all comparisons
visual_tile_dir <- "Figures/Visual_gene_plots"

dir.create(
  visual_tile_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## prepare one row per gene per comparison
visual_lfc_plot_data <- visual_DE_genes_export %>%
  mutate(
    gene_name = toupper(trimws(gene_name))
  ) %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    !is.na(log2FoldChange)
  ) %>%
  distinct(
    comparison,
    gene_name,
    .keep_all = TRUE
  ) %>%
  select(
    comparison,
    gene_name,
    log2FoldChange,
    padj
  )

## define comparison order
comparison_order <- c(
  "Ab_vs_Mz",
  "Ab_vs_Nb",
  "Ab_vs_On",
  "Ab_vs_Pn",
  "Mz_vs_Nb",
  "Mz_vs_On",
  "Mz_vs_Pn",
  "Nb_vs_On",
  "Nb_vs_Pn",
  "On_vs_Pn"
)

## add missing gene-comparison combinations
## missing combinations are left as NA rather than treated as zero
visual_lfc_complete <- visual_lfc_plot_data %>%
  complete(
    gene_name,
    comparison = comparison_order
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = comparison_order
    )
  )

## order genes by how often they were significant
gene_order <- c(
  # Opsin
  "OPN1SW2",
  
  # Phototransduction
  "GNAT1",
  "GNAT2",
  "GNGT1",
  "PDE6A",
  "PDE6C",
  "GUCA1A",
  "GUCA1B",
  
  # Development
  "CRX",
  "NRL",
  "NR2E3",
  
  # Visual cycle
  "RDH5"
)

visual_lfc_complete <- visual_lfc_complete %>%
  mutate(
    gene_name = factor(
      gene_name,
      levels = rev(gene_order)
    )
  )

## symmetrical colour limits
lfc_limit <- max(
  abs(visual_lfc_complete$log2FoldChange),
  na.rm = TRUE
)

## create tile plot
p_visual_lfc <- ggplot(
  visual_lfc_complete,
  aes(
    x = comparison,
    y = gene_name,
    fill = log2FoldChange
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = case_when(
        is.na(log2FoldChange) ~ "",
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        TRUE ~ ""
      )
    ),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -lfc_limit,
      lfc_limit
    ),
    na.value = "grey90",
    name = expression(log[2]~fold~change)
  ) +
  labs(
    title = "Differential expression of candidate visual genes",
    x = "Pairwise comparison",
    y = "Candidate visual gene",
    caption = paste(
      "Positive values indicate higher expression in the first species;",
      "negative values indicate higher expression in the second species.",
      "* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001"
    )
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 11
    ),
    panel.grid = element_blank()
  )

## display plot
p_visual_lfc

## save
ggsave(
  filename = file.path(
    visual_tile_dir,
    "Visual_genes_log2FC_tileplot.png"
  ),
  plot = p_visual_lfc,
  width = 12,
  height = 7,
  dpi = 600
)

## save pdf
ggsave(
  filename = file.path(
    visual_tile_dir,
    "Visual_genes_log2FC_tileplot.pdf"
  ),
  plot = p_visual_lfc,
  width = 12,
  height = 7
)

## save plotted data
write.csv(
  visual_lfc_complete,
  file.path(
    visual_tile_dir,
    "Visual_genes_log2FC_tileplot_data.csv"
  ),
  row.names = FALSE
)



