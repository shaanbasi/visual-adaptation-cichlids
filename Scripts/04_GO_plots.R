library(ggplot2)
library(stringr)
library(DESeq2)
library(dplyr)
library(tibble)
library(pheatmap)
library(grid)
library(tidyr)
library(purrr)


go_plot_dir <- "Figures/GO"

visual_plot_dir <- "Figures/Visual_genes"

visual_result_dir <- "Results/Visual_genes"

go_result_dir <- "Results/GO"

dir.create(
  go_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  visual_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  visual_result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  go_result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


## load shared visual system candidate gene list

visual_genes <- readr::read_csv(
  "Data/Reference/visual_candidate_genes.csv",
  show_col_types = FALSE
) %>%
  dplyr::pull(
    gene_name
  ) %>%
  trimws() %>%
  toupper() %>%
  unique()

visual_genes_lower <- tolower(
  visual_genes
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


## SPECIES-VS-REST GO ENRICHMENT PLOTS

species_vs_rest_go_dir <- file.path(
  "Results",
  "GO",
  "species_vs_rest"
)

species_vs_rest_go_plot_dir <- file.path(
  "Figures",
  "GO",
  "species_vs_rest"
)

dir.create(
  species_vs_rest_go_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


## read species-vs-rest GO results

species_vs_rest_GO <- readr::read_csv(
  file.path(
    species_vs_rest_go_dir,
    "GO_species_vs_rest_all_clean.csv"
  ),
  show_col_types = FALSE
)


## Very generic GO terms excluded from representative term selection
## Broader but biologically interpretable terms such as translation, ribosome, ATP metabolic process are retained

uninformative_GO_terms <- c(
  
  ## V broad molecular functions
  "binding",
  "protein binding",
  "small molecule binding",
  "ion binding",
  "catalytic activity",
  "hydrolase activity",
  "transferase activity",
  
  ## V broad cellular components
  "cytoplasm",
  "cytosol",
  "membrane",
  "organelle",
  "intracellular organelle",
  "intracellular membrane-bounded organelle",
  "organelle membrane",
  "organelle subcompartment",
  "endomembrane system",
  "intracellular",
  "intracellular anatomical structure",
  
  ## V broad biological processes
  "cellular process",
  "metabolic process",
  "primary metabolic process",
  "macromolecule metabolic process",
  "protein metabolic process",
  "catabolic process",
  "biological regulation",
  "regulation of biological process",
  "regulation of cellular process",
  "developmental process",
  "multicellular organismal process"
)


species_vs_rest_GO_plot_data <- species_vs_rest_GO %>%
  dplyr::filter(
    !term_name %in%
      uninformative_GO_terms,
    !is.na(p_value),
    p_value < 0.05
  )


## plotting function

plot_species_vs_rest_GO <- function(
    go_tbl,
    focal_species,
    direction,
    top_n = 15
) {
  
  plot_data <- go_tbl %>%
    dplyr::filter(
      species == focal_species,
      group == direction
    ) %>%
    dplyr::arrange(
      p_value
    ) %>%
    dplyr::slice_head(
      n = top_n
    )
  
  
  ## return nothing if no significant terms are available
  
  if (nrow(plot_data) == 0) {
    
    message(
      "No significant GO terms for ",
      focal_species,
      " | ",
      direction
    )
    
    return(NULL)
  }
  
  
  plot_data <- plot_data %>%
    dplyr::mutate(
      
      neg_log10_p =
        -log10(
          pmax(
            p_value,
            .Machine$double.xmin
          )
        ),
      
      term_name = factor(
        term_name,
        levels = rev(
          unique(
            term_name
          )
        )
      )
    )
  
  
  ## figure subtitle
  
  if (
    direction ==
    paste(
      focal_species,
      "higher"
    )
  ) {
    
    plot_subtitle <- paste(
      "Genes with higher expression in",
      focal_species,
      "relative to the remaining species"
    )
    
  } else {
    
    plot_subtitle <- paste(
      "Genes with lower expression in",
      focal_species,
      "relative to the remaining species"
    )
  }
  
  
  ggplot(
    plot_data,
    aes(
      x = neg_log10_p,
      y = term_name,
      size = intersection_size,
      color = source
    )
  ) +
    
    geom_point(
      alpha = 0.85
    ) +
    
    scale_size_continuous(
      name = "Genes"
    ) +
    
    labs(
      title = paste(
        focal_species,
        "GO enrichment"
      ),
      
      subtitle = plot_subtitle,
      
      x = expression(
        -log[10](
          "FDR-adjusted p-value"
        )
      ),
      
      y = NULL,
      
      color = "GO source"
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      
      plot.subtitle = element_text(
        size = 11
      ),
      
      axis.text.y = element_text(
        size = 9
      )
    )
}


## species order

species_order <- c(
  "Ab",
  "Mz",
  "Nb",
  "On",
  "Pn"
)



## PLOTS FOR GENES HIGHER IN THE FOCAL SPECIES

species_vs_rest_GO_plots <- list()


for (focal_species in species_order) {
  
  focal_group <- paste(
    focal_species,
    "higher"
  )
  
  
  plot_object <- plot_species_vs_rest_GO(
    go_tbl = species_vs_rest_GO_plot_data,
    focal_species = focal_species,
    direction = focal_group,
    top_n = 15
  )
  
  
  species_vs_rest_GO_plots[[focal_species]] <-
    plot_object
  
  
  if (!is.null(plot_object)) {
    
    ggsave(
      filename = file.path(
        species_vs_rest_go_plot_dir,
        paste0(
          "GO_",
          focal_species,
          "_higher_vs_rest.png"
        )
      ),
      plot = plot_object,
      width = 10,
      height = 7,
      dpi = 600
    )
    
    
    ggsave(
      filename = file.path(
        species_vs_rest_go_plot_dir,
        paste0(
          "GO_",
          focal_species,
          "_higher_vs_rest.pdf"
        )
      ),
      plot = plot_object,
      width = 10,
      height = 7
    )
  }
}




## PLOTS FOR GENES LOWER IN THE FOCAL SPECIES

species_vs_rest_lower_GO_plots <- list()


for (focal_species in species_order) {
  
  plot_object <- plot_species_vs_rest_GO(
    go_tbl = species_vs_rest_GO_plot_data,
    focal_species = focal_species,
    direction = "Other species higher",
    top_n = 15
  )
  
  
  species_vs_rest_lower_GO_plots[[focal_species]] <-
    plot_object
  
  
  if (!is.null(plot_object)) {
    
    ggsave(
      filename = file.path(
        species_vs_rest_go_plot_dir,
        paste0(
          "GO_",
          focal_species,
          "_lower_vs_rest.png"
        )
      ),
      plot = plot_object,
      width = 10,
      height = 7,
      dpi = 600
    )
    
    
    ggsave(
      filename = file.path(
        species_vs_rest_go_plot_dir,
        paste0(
          "GO_",
          focal_species,
          "_lower_vs_rest.pdf"
        )
      ),
      plot = plot_object,
      width = 10,
      height = 7
    )
  }
}




## VISION-RELATED SPECIES-VS-REST GO TERMS

vision_terms_species_vs_rest <- c(
  "visual perception",
  "phototransduction",
  "opsin",
  "rhodopsin",
  "retina",
  "retinal",
  "neural retina",
  "photoreceptor",
  "camera-type eye",
  "eye development",
  "lens development",
  "pigment",
  "pigmentation",
  "melanin",
  "melanogenesis",
  "chromatophore",
  "retinoid"
)


vision_pattern_species_vs_rest <- paste(
  vision_terms_species_vs_rest,
  collapse = "|"
)


species_vs_rest_visual_GO <-
  species_vs_rest_GO_plot_data %>%
  dplyr::filter(
    grepl(
      vision_pattern_species_vs_rest,
      term_name,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::arrange(
    species,
    p_value
  )


## display vision related terms

print(
  tibble::as_tibble(
    species_vs_rest_visual_GO
  ),
  n = Inf,
  width = Inf
)


## export vision related terms

readr::write_csv(
  species_vs_rest_visual_GO,
  file.path(
    species_vs_rest_go_dir,
    "GO_species_vs_rest_visual_terms.csv"
  )
)


## summarise vision related terms

species_vs_rest_visual_GO_summary <-
  species_vs_rest_visual_GO %>%
  dplyr::count(
    species,
    group,
    name = "n_visual_GO_terms"
  ) %>%
  dplyr::arrange(
    species,
    group
  )


print(
  species_vs_rest_visual_GO_summary,
  n = Inf
)


readr::write_csv(
  species_vs_rest_visual_GO_summary,
  file.path(
    species_vs_rest_go_dir,
    "GO_species_vs_rest_visual_summary.csv"
  )
)



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
  file.path(
    visual_result_dir,
    "DEGs_opsin_phototransduction_visual_genes.csv"
  ),
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
  file.path(
    visual_result_dir,
    "DEGs_opsin_phototransduction_visual_summary.csv"
  ),
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
  file.path(
    visual_result_dir,
    "DEGs_visual_gene_frequency.csv"
  ),
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
heatmap_dir <- visual_plot_dir

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
    visual_result_dir,
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
    visual_result_dir,
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
    visual_result_dir,
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
  file.path(
    go_result_dir,
    "GO_all_terms_ranked.csv"
  ),
  row.names = FALSE
)



## table showing candidate visual genes across all comparisons
visual_tile_dir <- visual_plot_dir

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

## add missing gene comp combinations
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
  ## opsin
  "OPN1SW2",
  
  ## Phototransduction
  "GNAT1",
  "GNAT2",
  "GNGT1",
  "PDE6A",
  "PDE6C",
  "GUCA1A",
  "GUCA1B",
  
  ## development
  "CRX",
  "NRL",
  "NR2E3",
  
  ## visual cycle
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
    visual_result_dir,
    "Visual_genes_log2FC_tileplot_data.csv"
  ),
  row.names = FALSE
)





## COMPARATIVE SPECIES-VS-REST GO BUBBLE PLOT (REDUNDANCY-REDUCED)

library(simplifyEnrichment)

species_vs_rest_GO <- readr::read_csv(
  file.path(
    "Results",
    "GO",
    "species_vs_rest",
    "GO_species_vs_rest_all_clean.csv"
  ),
  show_col_types = FALSE
)


## prepare significant GO results

species_order <- c(
  "Ab",
  "Mz",
  "Nb",
  "On",
  "Pn"
)

GO_root_terms <- c(
  "biological_process",
  "molecular_function",
  "cellular_component"
)

comparative_GO_data <- species_vs_rest_GO %>%
  dplyr::filter(
    !is.na(p_value),
    p_value < 0.05,
    !term_name %in% GO_root_terms
  ) %>%
  dplyr::mutate(
    
    species = factor(
      species,
      levels = species_order
    ),
    
    direction = dplyr::case_when(
      
      group == paste(
        as.character(species),
        "higher"
      ) ~ "Up-regulated",
      
      group == "Other species higher" ~
        "Down-regulated",
      
      TRUE ~ NA_character_
    ),
    
    direction = factor(
      direction,
      levels = c(
        "Up-regulated",
        "Down-regulated"
      )
    ),
    
    ontology = dplyr::recode(
      source,
      "GO:BP" = "BP",
      "GO:MF" = "MF",
      "GO:CC" = "CC"
    ),
    
    neg_log10_FDR = -log10(
      pmax(
        p_value,
        .Machine$double.xmin
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(direction),
    !is.na(fold_enrichment),
    is.finite(fold_enrichment),
    fold_enrichment > 0
  )


## Check significant enrichment by species and direction

GO_enrichment_check <- comparative_GO_data %>%
  dplyr::group_by(
    species,
    direction
  ) %>%
  dplyr::summarise(
    n_significant_GO_terms =
      dplyr::n_distinct(term_id),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    species,
    direction,
    fill = list(
      n_significant_GO_terms = 0
    )
  ) %>%
  dplyr::arrange(
    species,
    direction
  )

print(
  GO_enrichment_check,
  n = Inf
)


## Very generic GO terms excluded from representative term selection
## Broader but biologically interpretable terms such as translation, ribosome, ATP metabolic process are retained

uninformative_GO_terms <- c(
  
  ## V broad molecular functions
  "binding",
  "protein binding",
  "small molecule binding",
  "ion binding",
  "catalytic activity",
  "hydrolase activity",
  "transferase activity",
  
  ## V broad cellular components
  "cytoplasm",
  "cytosol",
  "membrane",
  "organelle",
  "intracellular organelle",
  "intracellular membrane-bounded organelle",
  "membrane-bounded organelle",
  "organelle membrane",
  "organelle subcompartment",
  "endomembrane system",
  "intracellular",
  "intracellular anatomical structure",
  
  ## V broad biological processes
  "cellular process",
  "metabolic process",
  "primary metabolic process",
  "macromolecule metabolic process",
  "protein metabolic process",
  "catabolic process",
  "biological regulation",
  "regulation of biological process",
  "regulation of cellular process",
  "developmental process",
  "multicellular organismal process"
)



## function to reduce semantic redundancy within one species x direction 
## x ontology combination

reduce_GO_subset <- function(
    go_subset,
    ontology_name
) {
  
  ## no terms
  if (nrow(go_subset) == 0) {
    return(tibble::tibble())
  }
  
  ## one term only
  if (nrow(go_subset) == 1) {
    
    return(
      go_subset %>%
        dplyr::mutate(
          semantic_cluster = 1L,
          broad_term =
            stringr::str_to_lower(
              stringr::str_squish(term_name)
            ) %in%
            stringr::str_to_lower(
              stringr::str_squish(uninformative_GO_terms)
            )
        )
    )
  }
  
  ## semantic similarity matrix
  similarity_matrix <-
    simplifyEnrichment::GO_similarity(
      go_subset$term_id,
      ont = ontology_name,
      db = "org.Dr.eg.db"
    )
  
  ## cluster semantically related terms
  cluster_table <-
    simplifyEnrichment::simplifyGO(
      similarity_matrix,
      method = "binary_cut",
      plot = FALSE,
      verbose = FALSE
    ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      term_id = id,
      semantic_cluster = cluster
    )
  
  ## attach semantic cluster identities
  go_subset %>%
    dplyr::inner_join(
      cluster_table,
      by = "term_id"
    ) %>%
    dplyr::mutate(
      broad_term =
        stringr::str_to_lower(
          stringr::str_squish(term_name)
        ) %in%
        stringr::str_to_lower(
          stringr::str_squish(uninformative_GO_terms)
        )
    )
}


## perform semantic reduction separately for each

GO_semantic_clusters <- comparative_GO_data %>%
  dplyr::group_by(
    species,
    direction,
    ontology
  ) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    reduced = purrr::map2(
      data,
      as.character(ontology),
      reduce_GO_subset
    )
  ) %>%
  dplyr::select(
    -data
  ) %>%
  tidyr::unnest(
    cols = reduced
  ) %>%
  dplyr::ungroup()


## select one informative representative per semantic cluster
## broad terms are deprioritised and of the remaining terms the most significant term is
## selected

GO_representatives <- GO_semantic_clusters %>%
  dplyr::filter(
    !broad_term
  ) %>%
  dplyr::group_by(
    species,
    direction,
    ontology,
    semantic_cluster
  ) %>%
  dplyr::arrange(
    p_value,
    dplyr::desc(fold_enrichment),
    dplyr::desc(intersection_size)
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup()



## select a manageable number of representatives

##up to 3 representative terms per ontology for each species x direction

GO_plot_terms <- GO_representatives %>%
  dplyr::group_by(
    species,
    direction,
    ontology
  ) %>%
  dplyr::arrange(
    p_value,
    dplyr::desc(fold_enrichment)
  ) %>%
  dplyr::slice_head(
    n = 3
  ) %>%
  dplyr::ungroup()


visual_GO_pattern <- paste0(
  "\\b(",
  paste(
    c(
      "visual perception",
      "phototransduction",
      "opsin",
      "rhodopsin",
      "retina",
      "retinal",
      "neural retina",
      "photoreceptor",
      "camera-type eye",
      "eye development",
      "lens development",
      "pigment",
      "pigmentation",
      "melanin",
      "melanogenesis",
      "chromatophore",
      "retinoid"
    ),
    collapse = "|"
  ),
  ")\\b"
)

visual_GO_plot_terms <- comparative_GO_data %>%
  dplyr::filter(
    stringr::str_detect(
      stringr::str_to_lower(term_name),
      stringr::regex(
        visual_GO_pattern,
        ignore_case = TRUE
      )
    )
  )

GO_plot_terms <- dplyr::bind_rows(
  GO_plot_terms,
  visual_GO_plot_terms
) %>%
  dplyr::distinct(
    species,
    direction,
    term_id,
    .keep_all = TRUE
  )


##Add plotting labels

GO_plot_data <- GO_plot_terms %>%
  dplyr::mutate(
    
    term_label = paste0(
      term_name,
      " (",
      ontology,
      ")"
    ),
    
    ## unique identifier prevents terms with the same name
    ## in different facets from interfering with ordering
    plot_term = paste(
      species,
      direction,
      term_id,
      sep = "___"
    )
  )


## Check what existed before reduction versus what is plotted

GO_significant_check <- comparative_GO_data %>%
  dplyr::count(
    species,
    direction,
    name = "n_significant_terms"
  ) %>%
  tidyr::complete(
    species = factor(
      species_order,
      levels = species_order
    ),
    direction = factor(
      c("Up-regulated", "Down-regulated"),
      levels = c(
        "Up-regulated",
        "Down-regulated"
      )
    ),
    fill = list(
      n_significant_terms = 0
    )
  )


GO_plotted_check <- GO_plot_data %>%
  dplyr::count(
    species,
    direction,
    name = "n_plotted_terms"
  ) %>%
  tidyr::complete(
    species = factor(
      species_order,
      levels = species_order
    ),
    direction = factor(
      c("Up-regulated", "Down-regulated"),
      levels = c(
        "Up-regulated",
        "Down-regulated"
      )
    ),
    fill = list(
      n_plotted_terms = 0
    )
  )

GO_selection_check <- GO_significant_check %>%
  dplyr::left_join(
    GO_plotted_check,
    by = c(
      "species",
      "direction"
    )
  ) %>%
  dplyr::arrange(
    species,
    direction
  )

print(
  GO_selection_check,
  n = Inf
)


## Order GO terms by fold enrichment

GO_plot_order <- GO_plot_data %>%
  dplyr::arrange(
    species,
    direction,
    fold_enrichment
  ) %>%
  dplyr::pull(
    plot_term
  ) %>%
  unique()

GO_plot_data <- GO_plot_data %>%
  dplyr::mutate(
    plot_term = factor(
      plot_term,
      levels = GO_plot_order
    )
  )



## FINAL COMPARATIVE GO BUBBLE PLOT

## factor order for final comparatice figure
## Species = columns
## DE direction = rows

GO_plot_data <- GO_plot_data %>%
  dplyr::mutate(
    species = factor(
      species,
      levels = c(
        "Ab",
        "Mz",
        "Nb",
        "On",
        "Pn"
      )
    ),
    
    direction = factor(
      direction,
      levels = c(
        "Up-regulated",
        "Down-regulated"
      )
    )
  )



## FIGURE 3
## Comparative species vs rest GO enrichment
##
## columns = focal species
## rows = DE direction
## x = fold enrichment
## y = representative GO term
## size = number of contributing DE genes
## colour = -log10(FDR)
##
## Mz is omitted because no significant GO enrichment detected in either DE direction


## PREPARE FINAL PLOT DATA

GO_plot_data_final <- GO_plot_data %>%
  dplyr::filter(
    species %in% c(
      "Ab",
      "Nb",
      "On",
      "Pn"
    )
  ) %>%
  dplyr::mutate(
    
    species = factor(
      species,
      levels = c(
        "Ab",
        "Nb",
        "On",
        "Pn"
      )
    ),
    
    direction = factor(
      direction,
      levels = c(
        "Up-regulated",
        "Down-regulated"
      )
    )
  )


## ORDER TERMS WITHIN EACH DIRECTION

term_order <- GO_plot_data_final %>%
  dplyr::group_by(
    direction,
    term_label
  ) %>%
  dplyr::summarise(
    max_fold_enrichment = max(
      fold_enrichment,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    direction,
    max_fold_enrichment
  ) %>%
  dplyr::pull(
    term_label
  ) %>%
  unique()


GO_plot_data_final <- GO_plot_data_final %>%
  dplyr::mutate(
    term_label = factor(
      term_label,
      levels = term_order
    )
  )


## CREATE FIGURE 3

p_GO_comparative_reduced <- ggplot(
  GO_plot_data_final,
  aes(
    x = fold_enrichment,
    y = term_label
  )
) +
  
  geom_point(
    aes(
      size = intersection_size,
      fill = neg_log10_FDR
    ),
    shape = 21,
    colour = "black",
    stroke = 0.3,
    alpha = 0.9
  ) +
  
  facet_grid(
    rows = vars(direction),
    cols = vars(species),
    scales = "free_y",
    space = "free_y",
    drop = TRUE
  ) +
  
  ## common x-axis makes comparisons meaningful
  scale_x_continuous(
    limits = c(
      0,
      35
    ),
    breaks = c(
      0,
      10,
      20,
      30
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +
  
  scale_fill_viridis_c(
    name = expression(
      -log[10]*"(FDR)"
    )
  ) +
  
  scale_size_continuous(
    name = "Genes",
    range = c(
      3,
      9
    )
  ) +
  
  labs(
    x = "Fold enrichment",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    ## Facet strips
    strip.background = element_rect(
      fill = "white",
      colour = "grey50",
      linewidth = 0.5
    ),
    
    strip.text.x = element_text(
      face = "bold",
      size = 11
    ),
    
    strip.text.y = element_text(
      face = "bold",
      size = 10
    ),
    
    ## Axis text
    axis.text.y = element_text(
      size = 8.5
    ),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(
        t = 8
      )
    ),
    
    ## box each facet
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.4
    ),
    
    ## Space between facets
    panel.spacing.x = grid::unit(
      0.35,
      "lines"
    ),
    
    panel.spacing.y = grid::unit(
      0.35,
      "lines"
    ),
    
    ## Legend
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    ## Margins
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


p_GO_comparative_reduced


##check
GO_plot_data_final %>%
  dplyr::filter(
    tolower(term_name) %in%
      tolower(uninformative_GO_terms)
  )


##SAVE
ggsave(
  filename = file.path(
    species_vs_rest_go_plot_dir,
    "Figure3_GO_comparative_bubbleplot.png"
  ),
  plot = p_GO_comparative_reduced,
  width = 12,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

## EXPORT REDUNDANCY REDUCTION RESULTS
readr::write_csv(
  GO_semantic_clusters,
  file.path(
    species_vs_rest_go_dir,
    "GO_semantic_clusters.csv"
  )
)

readr::write_csv(
  GO_representatives,
  file.path(
    species_vs_rest_go_dir,
    "GO_semantic_cluster_representatives.csv"
  )
)

readr::write_csv(
  GO_plot_data %>%
    dplyr::select(
      species,
      direction,
      ontology,
      term_id,
      term_name,
      p_value,
      neg_log10_FDR,
      intersection_size,
      query_size,
      term_size,
      effective_domain_size,
      fold_enrichment,
      semantic_cluster,
      broad_term
    ),
  file.path(
    species_vs_rest_go_dir,
    "GO_Figure3_plot_data.csv"
  )
  )
