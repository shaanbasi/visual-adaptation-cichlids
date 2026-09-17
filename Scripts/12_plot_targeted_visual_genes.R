library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

options(stringsAsFactors = FALSE)


##INPUT/OUTPOUT
input_dir <- file.path(
  "Results",
  "targeted_visual_genes"
)

figure_dir <- file.path(
  "Figures",
  "targeted_visual_genes"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


species_levels <- c(
  "Ab",
  "Mz",
  "Nb",
  "On",
  "Pn"
)



## LOAD CURATED VISUAL-GENE DE RESULTS
visual_DE <- readr::read_csv(
  file.path(
    input_dir,
    "curated_visual_genes_all_species_vs_rest.csv"
  ),
  show_col_types = FALSE
)


## CLEAN NAMES
visual_DE <- visual_DE %>%
  dplyr::mutate(
    
    gene_name =
      dplyr::coalesce(
        gene_name.y,
        gene_name.x
      ),
    
    species = factor(
      species,
      levels = species_levels
    ),
    
    significant =
      !is.na(padj) &
      padj < 0.05
  )



## ORDER VISUAL CATEGORIES

visual_category_levels <- c(
  "Opsin / photopigment",
  "Phototransduction",
  "CNG channel",
  "Photoreceptor recovery",
  "Visual cycle",
  "Photoreceptor regulation",
  "Retinal regulation"
)


visual_DE <- visual_DE %>%
  dplyr::mutate(
    visual_category = factor(
      visual_category,
      levels = visual_category_levels
    )
  )




## ORDER GENES WITHIN FUNCTIONAL CATEGORY

gene_order <- visual_DE %>%
  
  dplyr::group_by(
    visual_category,
    gene_name
  ) %>%
  
  dplyr::summarise(
    strongest_abs_log2FC =
      max(
        abs(
          log2FoldChange
        ),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    visual_category,
    dplyr::desc(
      strongest_abs_log2FC
    )
  )


visual_DE <- visual_DE %>%
  dplyr::mutate(
    gene_name = factor(
      gene_name,
      levels = rev(
        gene_order$gene_name
      )
    )
  )



## TARGETED VISUAL SYSTEM HEATMAP

p_visual_heatmap <- ggplot(
  visual_DE,
  aes(
    x = species,
    y = gene_name,
    fill = log2FoldChange
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  
  ## mark significant species-vs-rest results
  geom_point(
    data = visual_DE %>%
      dplyr::filter(
        significant
      ),
    shape = 21,
    size = 2.3,
    stroke = 0.5,
    fill = "white",
    color = "black"
  ) +
  
  facet_grid(
    visual_category ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_gradient2(
      name = "Shrunken\nlog2 fold change",
      low = "#3B4CC0",
      mid = "grey",
      high = "#B40426",
      midpoint = 0,
      limits = c(-6, 6),
      oob = scales::squish
  ) +
  
  labs(
    title =
      "Species-associated expression shifts in visual-system genes",
    
    subtitle =
      paste(
        "Curated visual-system genes",
        "(circles indicate FDR < 0.05)"
      ),
    
    x = "Focal species",
    
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),
    
    strip.placement = "outside",
    
    panel.spacing.y = grid::unit(
      0.5,
      "lines"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )


p_visual_heatmap

## SAVE
ggsave(
  filename = file.path(
    figure_dir,
    "Visual_system_species_vs_rest_heatmap.png"
  ),
  plot = p_visual_heatmap,
  width = 10,
  height = 13,
  dpi = 600
)