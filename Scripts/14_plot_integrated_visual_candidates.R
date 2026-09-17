library(dplyr)
library(readr)
library(ggplot2)

options(stringsAsFactors = FALSE)


## INPUT/OUTPUT
input_file <- file.path(
  "Results",
  "targeted_visual_genes",
  "WGCNA_integration",
  "visual_gene_final_DE_WGCNA_hub_integration.csv"
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


visual_category_levels <- c(
  "Opsin / photopigment",
  "Phototransduction",
  "CNG channel",
  "Photoreceptor recovery",
  "Visual cycle",
  "Photoreceptor regulation",
  "Retinal regulation"
)



## LOAD FINAL INTEGRATED RESULTS

visual_integrated <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)



## PREPARE DATA

visual_integrated <- visual_integrated %>%
  
  dplyr::mutate(
    
    species = factor(
      species,
      levels = species_levels
    ),
    
    visual_category = factor(
      visual_category,
      levels = visual_category_levels
    ),
    
    significant =
      !is.na(padj) &
      padj < 0.05,
    
    strongest_integrated_support =
      significant &
      dplyr::coalesce(
        species_associated_module,
        FALSE
      ) &
      dplyr::coalesce(
        supported_hub,
        FALSE
      ),
    
    evidence_symbol =
      dplyr::case_when(
        
        strongest_integrated_support ~
          "DE + WGCNA module + hub",
        
        significant ~
          "DE (FDR < 0.05)",
        
        TRUE ~
          NA_character_
      )
  )



## ORDER GENES WITHIN FUNCTIONAL GROUPS

gene_order <- visual_integrated %>%
  
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


visual_integrated <- visual_integrated %>%
  
  dplyr::mutate(
    
    gene_name = factor(
      gene_name,
      levels = rev(
        gene_order$gene_name
      )
    )
  )



## INTEGRATED HEATMAP

p_visual_integrated <- ggplot(
  visual_integrated,
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
  
  ## significant DE result
  geom_point(
    data = visual_integrated %>%
      dplyr::filter(
        significant &
          !strongest_integrated_support
      ),
    
    aes(
      shape = "DE (FDR < 0.05)"
    ),
    
    size = 2.5,
    stroke = 0.6,
    fill = "white",
    color = "black"
  ) +
  
  ## strongest integrated evidence
  geom_point(
    data = visual_integrated %>%
      dplyr::filter(
        strongest_integrated_support
      ),
    
    aes(
      shape = "DE + WGCNA module + hub"
    ),
    
    size = 4.2,
    stroke = 0.8,
    fill = "black",
    color = "black"
  ) +
  
  ## functional groups
  facet_grid(
    visual_category ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  ## fold-change colour scale
  scale_fill_gradient2(
    name = expression(
      atop("Shrunken", log[2] * " fold change")
    ),
    low = "#3B4CC0",
    mid = "grey85",
    high = "#B40426",
    midpoint = 0,
    limits = c(-6, 6),
    oob = scales::squish
  ) +
  
  ## evidence symbols
  scale_shape_manual(
    name = "Evidence",
    values = c(
      "DE (FDR < 0.05)" = 21,
      "DE + WGCNA module + hub" = 23
    )
  ) +
  
  ## italicise gene symbols
  scale_y_discrete(
    labels = function(x) {
      parse(text = paste0("italic('", x, "')"))
    }
  ) +
  
  ## make legend match plotted symbols
  guides(
    shape = guide_legend(
      override.aes = list(
        fill = c(
          "white",
          "black"
        ),
        size = c(
          2.5,
          4.2
        )
      )
    )
  ) +
  
  labs(
    x = "Focal species",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
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
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )


p_visual_integrated



## SAVE

ggsave(
  filename = file.path(
    figure_dir,
    "Figure5_visual_system_integrated_DE_WGCNA_heatmap.png"
  ),
  plot = p_visual_integrated,
  width = 11,
  height = 13,
  dpi = 600
)
