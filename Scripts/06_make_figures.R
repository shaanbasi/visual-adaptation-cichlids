
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


##Volcano plots - editted for poster
make_volcano_plot <- function(res_tbl,
                              species_a,
                              species_b,
                              alpha = 0.05,
                              labels_per_side = 5) {
  
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  volcano_df <- res_tbl %>%
    mutate(
      neg_log10_padj = -log10(padj)
    )
  
  label_genes <- volcano_df %>%
    filter(!is.na(padj), padj < alpha, sig != "Not significant") %>%
    group_by(sig) %>%
    arrange(padj) %>%
    slice_head(n = labels_per_side) %>%
    ungroup()
  
  if (nrow(label_genes) > 0) {
    x_max <- max(abs(label_genes$log2FoldChange), na.rm = TRUE)
    y_max <- max(label_genes$neg_log10_padj, na.rm = TRUE)
    x_lim <- ceiling(x_max) + 1
    y_lim <- ceiling(y_max) + 2
  } else {
    x_lim <- 5
    y_lim <- 5
  }
  
  color_values <- setNames(
    c("#9ecae1", "#fdae6b", "grey80"),
    c(
      paste(species_a, "higher"),
      paste(species_b, "higher"),
      "Not significant"
    )
  )
  
  p <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    geom_point(aes(color = sig), alpha = 0.9, size = 2.5) +
    
    geom_text_repel(
      data = label_genes,
      aes(label = gene_name),
      size = 6,
      max.overlaps = 50,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "black",
      min.segment.length = 0
    ) +
    
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    scale_color_manual(
      values = color_values,
      name = "Expression change"
    ) +
    
    coord_cartesian(
      xlim = c(-x_lim, x_lim),
      ylim = c(0, y_lim)
    ) +
    
    theme_minimal(base_size = 16) +
    
    labs(
      title = paste("Volcano Plot:", species_a, "vs", species_b),
      x = "Log2 fold change",
      y = "−log10(adj. p)"
    ) +
    
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      
      legend.key.size = unit(0.6, "cm")
    )
  
  list(
    plot = p,
    volcano_df = volcano_df,
    label_genes = label_genes
  )
}



## Mz vs Ab
mz_ab <- make_pairwise_res(dds, species_a = "Mz", species_b = "Ab")

res_tbl_Mz_vs_Ab <- mz_ab$results
deg_tbl_Mz_vs_Ab <- mz_ab$degs

res_tbl_Mz_vs_Ab %>%
  select(gene_name, log2FoldChange, lfcSE, pvalue, padj, sig) %>%
  head(20)

summary(res_tbl_Mz_vs_Ab$padj)
table(res_tbl_Mz_vs_Ab$sig)

## clean DEG table for export
deg_tbl_Mz_vs_Ab_clean <- deg_tbl_Mz_vs_Ab %>%
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

write.csv(
  deg_tbl_Mz_vs_Ab_clean,
  "DESeq2_Mz_vs_Ab_DEGs_FDR0.05.csv",
  row.names = FALSE
)

## volcano plot
volcano_mz_ab <- make_volcano_plot(
  res_tbl = res_tbl_Mz_vs_Ab,
  species_a = "Mz",
  species_b = "Ab",
  alpha = 0.05,
  labels_per_side = 5
)

p_volcano_mz_ab <- volcano_mz_ab$plot
p_volcano_mz_ab

ggsave(
  filename = "Volcano_Mz_vs_Ab_clean.png",
  plot = p_volcano_mz_ab,
  width = 8,
  height = 6,
  dpi = 300
)


## Mz vs On
## this now works because releveling is automatic
mz_on <- make_pairwise_res(dds, species_a = "Mz", species_b = "On")

res_tbl_Mz_vs_On <- mz_on$results
deg_tbl_Mz_vs_On <- mz_on$degs

res_tbl_Mz_vs_On %>%
  select(gene_name, log2FoldChange, lfcSE, pvalue, padj, sig) %>%
  head(20)

summary(res_tbl_Mz_vs_On$padj)
table(res_tbl_Mz_vs_On$sig)

deg_tbl_Mz_vs_On_clean <- deg_tbl_Mz_vs_On %>%
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

write.csv(
  deg_tbl_Mz_vs_On_clean,
  "DESeq2_Mz_vs_On_DEGs_FDR0.05.csv",
  row.names = FALSE
)

volcano_mz_on <- make_volcano_plot(
  res_tbl = res_tbl_Mz_vs_On,
  species_a = "Mz",
  species_b = "On",
  alpha = 0.05,
  labels_per_side = 5
)

p_volcano_mz_on <- volcano_mz_on$plot
p_volcano_mz_on

ggsave(
  filename = "Volcano_Mz_vs_On_clean.png",
  plot = p_volcano_mz_on,
  width = 8,
  height = 6,
  dpi = 300
)



## LOOP THROUGH ALL PAIRWISE COMPARISONS
species_levels <- c("Ab", "Mz", "Nb", "On", "Pn")
comparisons <- combn(species_levels, 2, simplify = FALSE)

all_results_list <- list()
all_deg_list <- list()
all_volcano_list <- list()

for (comp in comparisons) {
  
  species_a <- comp[1]
  species_b <- comp[2]
  comparison_name <- paste0(species_a, "_vs_", species_b)
  
  message("Running comparison: ", comparison_name)
  
  pair_res <- make_pairwise_res(
    dds = dds,
    species_a = species_a,
    species_b = species_b,
    alpha = 0.05
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
  
  ## save all results and DEGs
  write.csv(
    res_tbl,
    paste0("DESeq2_", comparison_name, "_all_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    deg_tbl_clean,
    paste0("DESeq2_", comparison_name, "_DEGs_FDR0.05.csv"),
    row.names = FALSE
  )
  
  ## make volcano plot
  volcano_out <- make_volcano_plot(
    res_tbl = res_tbl,
    species_a = species_a,
    species_b = species_b,
    alpha = 0.05,
    labels_per_side = 5
  )
  
  p <- volcano_out$plot
  all_volcano_list[[comparison_name]] <- p
  
  ggsave(
    filename = paste0("Volcano_", comparison_name, "_clean.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

## combine all DEGs into one master table
all_pairwise_DEGs <- bind_rows(all_deg_list)

write.csv(
  all_pairwise_DEGs,
  "DESeq2_all_pairwise_DEGs_FDR0.05.csv",
  row.names = FALSE
)


##BAR PLOT - how many genes 
mz_on <- make_pairwise_res(dds, "Mz", "On")

res_tbl_Mz_vs_On <- mz_on$results
deg_tbl_Mz_vs_On <- mz_on$degs

table(res_tbl_Mz_vs_On$sig)

deg_counts <- res_tbl_Mz_vs_On %>%
  dplyr::filter(sig != "Not significant") %>%
  dplyr::count(sig)

deg_counts <- deg_counts %>%
  mutate(
    sig = recode(sig,
                 "Mz higher" = "Mz",
                 "On higher" = "On")
  )

library(ggplot2)

bar_mz_on <- ggplot(deg_counts, aes(x = sig, y = n, fill = sig)) +
  geom_bar(stat = "identity", width = 0.6) +
  
  scale_fill_manual(values = c(
    "Mz" = "#9ecae1",   # light blue
    "On" = "#fdae6b"    # light orange
  )) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Differentially Expressed Genes (Mz vs On)",
    x = "Upregulation",
    y = "Number of genes"
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )

bar_mz_on



## check
class(res_tbl_Mz_vs_On)
## not list

find("count")

unique(deg_counts$sig)
