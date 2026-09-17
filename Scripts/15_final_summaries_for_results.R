##Total no. of significant DE genes - to use in write-up

species_vs_rest_dir <- file.path(output_dir, "species_vs_rest")

species <- c("Ab", "Mz", "Nb", "On", "Pn")

DE_summary <- lapply(species, function(sp) {
  
  x <- readr::read_csv(
    file.path(
      species_vs_rest_dir,
      paste0("DESeq2_", sp, "_vs_rest_DEGs_FDR0.05.csv")
    ),
    show_col_types = FALSE
  )
  
  data.frame(
    species = sp,
    total_DE = nrow(x),
    higher_in_focal = sum(x$log2FoldChange > 0, na.rm = TRUE),
    lower_in_focal  = sum(x$log2FoldChange < 0, na.rm = TRUE)
  )
  
}) %>%
  dplyr::bind_rows()

print(DE_summary)


## SAVE
readr::write_csv(
  DE_summary,
  file.path(species_vs_rest_dir, "DESeq2_species_vs_rest_DE_counts.csv")
)




## Read module sizes
module_sizes <- readr::read_csv(
  "Results/WGCNA/module_sizes.csv",
  show_col_types = FALSE
)

print(module_sizes)
names(module_sizes)

print(module_sizes, n = Inf)


## Total no. of genes
sum(module_sizes$module_size)

datExpr <- readRDS("Data/processed/datExpr.rds")
dim(datExpr)


## Two significant relationships

module_species <- readr::read_csv(
  "Results/WGCNA/module_species_results.csv",
  show_col_types = FALSE
)

names(module_species)

module_species %>%
  dplyr::filter(
    (module == "pink" & species == "On") |
      (module == "lightyellow" & species == "Pn")
  ) %>%
  dplyr::select(
    module,
    species,
    correlation,
    pvalue,
    FDR,
    module_size
  )

##check
module_species %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  dplyr::select(
    module,
    species,
    correlation,
    pvalue,
    FDR,
    module_size
  )




## Major GO function in the On-Pink module

module_GO_all <- readr::read_csv(
  "Results/WGCNA/module_GO/WGCNA_selected_modules_GO_all.csv",
  show_col_types = FALSE
)

names(module_GO_all)

module_GO_all %>%
  dplyr::filter(
    species == "On",
    module == "pink"
  ) %>%
  print(n = Inf)




## Strong visual gene patterns for figure 5

visual_DE <- readr::read_csv(
  "Results/DESeq2/species_vs_rest/DESeq2_species_vs_rest_visual_candidate_results.csv",
  show_col_types = FALSE
)

names(visual_DE)


visual_summary <- readr::read_csv(
  "Results/DESeq2/species_vs_rest/DESeq2_species_vs_rest_visual_candidate_summary.csv",
  show_col_types = FALSE
)

print(visual_summary, n = Inf)





## Confirm evidence for rdh12 and rdh8a

visual_hubs <- readr::read_csv(
  "Results/WGCNA/hub_genes/12_visual_candidate_hub_overlap.csv",
  show_col_types = FALSE
)

names(visual_hubs)

visual_hubs %>%
  dplyr::filter(
    gene_id %in% c("rdh12", "rdh8a")
  ) %>%
  print(n = Inf)





