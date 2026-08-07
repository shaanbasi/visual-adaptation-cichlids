library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)

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

p_selected_modules <- ggplot(
  selected_ME_data,
  aes(
    x = species,
    y = eigengene
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 3
  ) +
  facet_wrap(
    ~ module,
    scales = "free_y",
    labeller = labeller(
      module = function(x) {
        paste(
          x,
          "module"
        )
      }
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    title = "Species-associated WGCNA module eigengenes",
    x = "Species",
    y = "Module eigengene"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )
## Pink moduel eigengene is strongly positive in On, negative
## in every other species.
## Lightyellow module - only Pn has consistently positive
## eigengene values.
## This is what WGCNA module-trait analysis predicted.



## DISPLAY PLOT

p_selected_modules


## SAVE PLOT

ggsave(
  filename = file.path(
    figure_dir,
    "WGCNA_selected_module_eigengenes.png"
  ),
  plot = p_selected_modules,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  filename = file.path(
    figure_dir,
    "WGCNA_selected_module_eigengenes.pdf"
  ),
  plot = p_selected_modules,
  width = 10,
  height = 6
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