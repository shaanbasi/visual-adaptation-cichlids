library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)
library(patchwork)

options(stringsAsFactors = FALSE)


## INPUT AND OUTPUT SETTINGS

wgcna_dir <- "Results/WGCNA"

hub_dir <- file.path(
  wgcna_dir,
  "hub_genes"
)

output_dir <- file.path(
  wgcna_dir,
  "module_expression"
)

figure_dir <- "Figures/WGCNA"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## LOAD WGCNA OBJECTS

load(
  file.path(
    wgcna_dir,
    "WGCNA_network_objects.RData"
  )
)

load(
  file.path(
    hub_dir,
    "WGCNA_hub_gene_objects.RData"
  )
)


## CHECK REQUIRED OBJECTS

required_objects <- c(
  "MEs",
  "speciesTraits",
  "selectedModuleSpeciesPairs"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}


## MATCH SAMPLE ORDER

speciesTraits <- speciesTraits[
  rownames(MEs),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(MEs),
    rownames(speciesTraits)
  )
)


## IDENTIFY SELECTED MODULE EIGENGENES

selected_modules <- unique(
  selectedModuleSpeciesPairs$module
)

selected_ME_columns <- paste0(
  "ME",
  selected_modules
)

missing_ME_columns <- setdiff(
  selected_ME_columns,
  colnames(MEs)
)

if (length(missing_ME_columns) > 0) {
  stop(
    "Missing module eigengene columns: ",
    paste(
      missing_ME_columns,
      collapse = ", "
    )
  )
}


## BUILD SAMPLE-LEVEL MODULE EIGENGENE TABLE

selected_ME_data <- MEs[
  ,
  selected_ME_columns,
  drop = FALSE
] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "sample"
  ) %>%
  dplyr::mutate(
    species = as.character(
      speciesTraits[
        sample,
        "species"
      ]
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      selected_ME_columns
    ),
    names_to = "module",
    values_to = "eigengene"
  ) %>%
  dplyr::mutate(
    module = sub(
      "^ME",
      "",
      module
    ),
    species = factor(
      species,
      levels = c(
        "Ab",
        "Mz",
        "Nb",
        "On",
        "Pn"
      )
    )
  )


## ADD ASSOCIATED SPECIES INFORMATION

selected_ME_data <- selected_ME_data %>%
  dplyr::left_join(
    selectedModuleSpeciesPairs %>%
      dplyr::select(
        associated_species = species,
        module,
        expression_pattern,
        module_trait_correlation = correlation,
        module_trait_FDR = FDR
      ),
    by = "module"
  )


## SUMMARISE MODULE EIGENGENES BY SPECIES

selected_ME_summary <- selected_ME_data %>%
  dplyr::group_by(
    module,
    associated_species,
    expression_pattern,
    species
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_eigengene = mean(
      eigengene,
      na.rm = TRUE
    ),
    sd_eigengene = sd(
      eigengene,
      na.rm = TRUE
    ),
    se_eigengene = sd_eigengene /
      sqrt(n),
    .groups = "drop"
  )


## EXPORT SAMPLE-LEVEL AND SUMMARY TABLES

readr::write_csv(
  selected_ME_data,
  file.path(
    output_dir,
    "selected_module_eigengenes_by_sample.csv"
  )
)

readr::write_csv(
  selected_ME_summary,
  file.path(
    output_dir,
    "selected_module_eigengenes_by_species.csv"
  )
)


## PLOT SAMPLE-LEVEL MODULE EIGENGENES

## helper function for individual module plots
plot_module_eigengene <- function(
    module_name
) {

  module_data <- selected_ME_data %>%
    dplyr::filter(
      module == module_name
    )

  ggplot(
    module_data,
    aes(
      x = species,
      y = eigengene
    )
  ) +

    ## individual biological samples
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.8
    ) +

    ## species mean
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 3.8
    ) +

    ## mean ± SE
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15,
      linewidth = 0.7
    ) +

    ## zero reference
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +

    labs(
      x = NULL,
      y = "Module eigengene"
    ) +

    theme_classic(
      base_size = 12
    ) +

    theme(
      axis.text.x = element_text(
        face = "bold",
        size = 10
      ),

      axis.text.y = element_text(
        size = 9
      ),

      axis.title.y = element_text(
        size = 11
      ),

      plot.margin = margin(
        8, 8, 8, 8
      )
    )
}

## panel B
p_magenta <- plot_module_eigengene(
  "magenta"
)


## panel C
p_lightgreen <- plot_module_eigengene(
  "lightgreen"
)


p_magenta <- p_magenta +
  facet_wrap(
    ~ module,
    labeller = labeller(
      module = c(
        "magenta" = "Magenta module"
      )
    )
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.title.y = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      6, 6, 6, 6
    )
  )


p_lightgreen <- p_lightgreen +
  facet_wrap(
    ~ module,
    labeller = labeller(
      module = c(
        "lightgreen" = "Lightgreen module"
      )
    )
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.title.y = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      6, 6, 6, 6
    )
  )


p_magenta
p_lightgreen


ggsave(
  filename = file.path(
    figure_dir,
    "Figure4B_magenta_module_eigengene.png"
  ),
  plot = p_magenta,
  width = 4,
  height = 3.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure4C_lightgreen_module_eigengene.png"
  ),
  plot = p_lightgreen,
  width = 4,
  height = 3.5,
  dpi = 600,
  bg = "white"
)


## FINAL FIGURE 4

## Reconstruct panel A from saved WGCNA results

## Hierarchically cluster non-grey modules according to their species-correl profiles
module_cor_mat <- moduleTraitCor[
  rownames(moduleTraitCor) != "MEgrey",
  c("Ab", "Mz", "Nb", "On", "Pn"),
  drop = FALSE
]


## Cluster modules using Euclidean distance between their five-species correlation profiles

module_row_clustering <- hclust(
  dist(module_cor_mat),
  method = "complete"
)


## Extract clustered module order
module_order <- rownames(module_cor_mat)[
  module_row_clustering$order
]


## Remove ME prefix so names match moduleSpeciesResults$module
module_order <- sub(
  "^ME",
  "",
  module_order
)


## Prepare heatmap data

heatmap_df <- moduleSpeciesResults %>%
  
  dplyr::filter(
    module != "grey"
  ) %>%
  
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
    
    module = factor(
      module,
      levels = rev(
        module_order
      )
    ),
    
    label = sprintf(
      "%.2f",
      correlation
    ),
    
    significant =
      !is.na(FDR) &
      FDR < 0.05,
    
    label_sig = ifelse(
      significant,
      paste0(
        label,
        "*"
      ),
      label
    )
  )

p_module_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = species,
    y = module,
    fill = correlation
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label = label_sig,
      fontface = ifelse(
        significant,
        "bold",
        "plain"
      )
    ),
    size = 3
  ) +
  
  scale_fill_gradient2(
    name = "Correlation",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 8
    ),
    
    axis.ticks = element_blank(),
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      6, 6, 6, 6
    )
  )


## Combine A above B and C
figure4 <- p_module_heatmap /
  (
    p_magenta |
      p_lightgreen
  ) +
  
  patchwork::plot_layout(
    heights = c(
      1.8,
      1
    )
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


## display
figure4


## save final Figure 4
ggsave(
  filename = file.path(
    figure_dir,
    "Figure4_WGCNA_module_species.png"
  ),
  plot = figure4,
  width = 8.2,
  height = 10.5,
  dpi = 600,
  bg = "white"
)





## FINAL SUMMARY

cat(
  "\n============================================\n"
)

cat(
  "Selected-module expression analysis complete\n"
)

cat(
  "Modules analysed:",
  paste(
    selected_modules,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Results saved to:",
  output_dir,
  "\n"
)

cat(
  "Figures saved to:",
  figure_dir,
  "\n"
)

cat(
  "============================================\n\n"
)


