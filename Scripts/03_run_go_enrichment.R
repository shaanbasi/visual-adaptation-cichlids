library(dplyr)
library(purrr)
library(gprofiler2)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)



output_dir <- "Results/GO"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## GO enrichment for all pairwise comparisons

run_go_for_pair <- function(res_tbl, deg_tbl, comparison_name) {
  
  species <- strsplit(
    comparison_name,
    "_vs_",
    fixed = TRUE
  )[[1]]
  
  if (length(species) != 2) {
    stop(
      "Comparison name must follow the format 'speciesA_vs_speciesB': ",
      comparison_name
    )
  }
  
  species_a <- species[1]
  species_b <- species[2]
  
  group_a <- paste(species_a, "higher")
  group_b <- paste(species_b, "higher")
  
  
  ## genes with higher expression in species A
  genes_a_up <- deg_tbl %>%
    filter(
      sig == group_a,
      !is.na(gene_name),
      trimws(gene_name) != ""
    ) %>%
    pull(gene_name) %>%
    unique()
  
  
  ## genes with higher expression in species B
  genes_b_up <- deg_tbl %>%
    filter(
      sig == group_b,
      !is.na(gene_name),
      trimws(gene_name) != ""
    ) %>%
    pull(gene_name) %>%
    unique()
  
  
  ## rna-seq background: all genes with a deseq2 p-value
  background_genes <- res_tbl %>%
    filter(
      !is.na(gene_name),
      trimws(gene_name) != "",
      !is.na(pvalue)
    ) %>%
    pull(gene_name) %>%
    unique()
  
  
  ## ensure query genes are present in the background
  genes_a_up <- intersect(
    genes_a_up,
    background_genes
  )
  
  genes_b_up <- intersect(
    genes_b_up,
    background_genes
  )
  
  
  ##helper function for one expression direction
  run_single_gost <- function(query_genes, group_name) {
    
    if (length(query_genes) == 0) {
      
      message(
        "No DE genes for ",
        comparison_name,
        " | ",
        group_name
      )
      
      return(tibble())
    }
    
    gost_result <- tryCatch(
      
      gost(
        query = query_genes,
        organism = "drerio",
        custom_bg = background_genes,
        domain_scope = "custom",
        sources = c(
          "GO:BP",
          "GO:MF",
          "GO:CC"
        ),
        correction_method = "fdr",
        evcodes = TRUE,
        user_threshold = 0.05,
        significant = TRUE
      ),
      
      error = function(e) {
        
        warning(
          "g:Profiler failed for ",
          comparison_name,
          " | ",
          group_name,
          ": ",
          conditionMessage(e)
        )
        
        return(NULL)
      }
    )
    
    
    ## handle comparisons with no significant enrichment
    if (
      is.null(gost_result) ||
      is.null(gost_result$result) ||
      nrow(gost_result$result) == 0
    ) {
      
      message(
        "No significant GO terms for ",
        comparison_name,
        " | ",
        group_name
      )
      
      return(tibble())
    }
    
    
    gost_result$result %>%
      mutate(
        group = group_name,
        comparison = comparison_name
      )
  }
  
  
  ## run the two directions separately
  go_a <- run_single_gost(
    query_genes = genes_a_up,
    group_name = group_a
  )
  
  go_b <- run_single_gost(
    query_genes = genes_b_up,
    group_name = group_b
  )
  
  
  bind_rows(
    go_a,
    go_b
  )
}


## run GO enrichment across all comparisons

all_GO <- map_dfr(
  names(all_results_list),
  ~ run_go_for_pair(
    res_tbl = all_results_list[[.x]],
    deg_tbl = all_deg_list[[.x]],
    comparison_name = .x
  )
)
## no statistically significant GO enrichment was detected after FDR correction


##keep only GO categories
all_GO_clean <- all_GO %>%
  filter(source %in% c("GO:BP", "GO:MF", "GO:CC"))

## export
all_GO_export <- all_GO_clean %>%
  select(
    comparison,
    group,
    term_id,
    term_name,
    source,
    p_value,
    intersection_size,
    query_size
  )

write.csv(
  all_GO_export,
  file.path(
    output_dir,
    "GO_all_pairwise_clean.csv"
  ),
  row.names = FALSE
)


##vision and pigment related GO terms
vision_terms <- c(
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

vision_pattern <- paste(
  vision_terms,
  collapse = "|"
)

vision_GO_all <- all_GO_clean %>%
  filter(
    grepl(
      vision_pattern,
      term_name,
      ignore.case = TRUE
    )
  )

vision_GO_export <- vision_GO_all %>%
  select(
    comparison,
    group,
    term_id,
    term_name,
    source,
    p_value,
    intersection_size,
    query_size
  )

write.csv(
  vision_GO_export,
  file.path(
    output_dir,
    "GO_vision_all_pairwise.csv"
  ),
  row.names = FALSE
)


##summarise no of GO terms 
vision_go_by_comparison <- vision_GO_export %>%
  dplyr::count(
    comparison,
    group,
    name = "n_GO_terms"
  )

print(vision_go_by_comparison)


## summarise most frequent vision-related GO terms
vision_go_frequency <- vision_GO_export %>%
  dplyr::count(
    term_name,
    sort = TRUE,
    name = "n_comparisons"
  )

print(vision_go_frequency)


## check total number of enriched GO terms
print(
  nrow(all_GO_export)
)


## summarise all enriched terms by comparison and direction
all_go_by_comparison <- all_GO_export %>%
  dplyr::count(
    comparison,
    group,
    name = "n_GO_terms"
  )

print(all_go_by_comparison)



## GO ENRICHMENT FOR SPECIES VS REST COMPARISONS

species_vs_rest_go_dir <- file.path(
  output_dir,
  "species_vs_rest"
)

dir.create(
  species_vs_rest_go_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


## helper function for one species-vs-rest comparison

run_go_for_species_vs_rest <- function(
    res_tbl,
    deg_tbl,
    focal_species
) {
  
  focal_group <- paste(
    focal_species,
    "higher"
  )
  
  rest_group <- "Other species higher"
  
  
  ## genes with higher expression in the focal species
  
  genes_focal_up <- deg_tbl %>%
    dplyr::filter(
      sig == focal_group,
      !is.na(gene_name),
      trimws(gene_name) != ""
    ) %>%
    dplyr::pull(
      gene_name
    ) %>%
    unique()
  
  
  ## genes with higher expression in the remaining species
  
  genes_rest_up <- deg_tbl %>%
    dplyr::filter(
      sig == rest_group,
      !is.na(gene_name),
      trimws(gene_name) != ""
    ) %>%
    dplyr::pull(
      gene_name
    ) %>%
    unique()
  
  
  ## RNA-seq background - all genes tested for this species vs rest contrast
  
  background_genes <- res_tbl %>%
    dplyr::filter(
      !is.na(gene_name),
      trimws(gene_name) != "",
      !is.na(pvalue)
    ) %>%
    dplyr::pull(
      gene_name
    ) %>%
    unique()
  
  
  ## ensure query genes are represented in the background
  
  genes_focal_up <- intersect(
    genes_focal_up,
    background_genes
  )
  
  genes_rest_up <- intersect(
    genes_rest_up,
    background_genes
  )
  
  
  ## helper function for one expression direction
  
  run_single_gost <- function(
    query_genes,
    group_name
  ) {
    
    if (length(query_genes) == 0) {
      
      message(
        "No DE genes for ",
        focal_species,
        " vs rest | ",
        group_name
      )
      
      return(
        tibble::tibble()
      )
    }
    
    
    message(
      "Running GO enrichment for ",
      focal_species,
      " vs rest | ",
      group_name,
      " | ",
      length(query_genes),
      " genes"
    )
    
    
    gost_result <- tryCatch(
      
      gprofiler2::gost(
        query = query_genes,
        organism = "drerio",
        custom_bg = background_genes,
        domain_scope = "custom",
        sources = c(
          "GO:BP",
          "GO:MF",
          "GO:CC"
        ),
        correction_method = "fdr",
        evcodes = TRUE,
        user_threshold = 0.05,
        significant = TRUE
      ),
      
      error = function(e) {
        
        warning(
          "g:Profiler failed for ",
          focal_species,
          " vs rest | ",
          group_name,
          ": ",
          conditionMessage(e)
        )
        
        return(NULL)
      }
    )
    
    
    if (
      is.null(gost_result) ||
      is.null(gost_result$result) ||
      nrow(gost_result$result) == 0
    ) {
      
      message(
        "No significant GO terms for ",
        focal_species,
        " vs rest | ",
        group_name
      )
      
      return(
        tibble::tibble()
      )
    }
    
    
    gost_result$result %>%
      dplyr::mutate(
        species = focal_species,
        
        comparison = paste0(
          focal_species,
          "_vs_rest"
        ),
        
        group = group_name
      )
  }
  
  
  ## run focal-species-higher genes
  
  go_focal <- run_single_gost(
    query_genes = genes_focal_up,
    group_name = focal_group
  )
  
  
  ## run genes higher in the remaining species
  
  go_rest <- run_single_gost(
    query_genes = genes_rest_up,
    group_name = rest_group
  )
  
  
  dplyr::bind_rows(
    go_focal,
    go_rest
  )
}


## run GO enrichment across all five species-vs-rest comparisons

species_vs_rest_GO <- purrr::map_dfr(
  names(
    species_vs_rest_results
  ),
  function(focal_species) {
    
    run_go_for_species_vs_rest(
      res_tbl =
        species_vs_rest_results[[focal_species]],
      
      deg_tbl =
        species_vs_rest_degs[[focal_species]],
      
      focal_species =
        focal_species
    )
  }
)


## keep only GO categories

species_vs_rest_GO_clean <- species_vs_rest_GO %>%
  dplyr::filter(
    source %in% c(
      "GO:BP",
      "GO:MF",
      "GO:CC"
    )
  )


## calculate fold enrichment and export

species_vs_rest_GO_export <- species_vs_rest_GO_clean %>%
  dplyr::mutate(
    fold_enrichment =
      (intersection_size / query_size) /
      (term_size / effective_domain_size)
  ) %>%
  dplyr::select(
    species,
    comparison,
    group,
    term_id,
    term_name,
    source,
    p_value,
    intersection,
    intersection_size,
    query_size,
    term_size,
    effective_domain_size,
    fold_enrichment
  )


write.csv(
  species_vs_rest_GO_export,
  file.path(
    species_vs_rest_go_dir,
    "GO_species_vs_rest_all_clean.csv"
  ),
  row.names = FALSE
)



## Figure 3 terms to investigate further

figure3_gene_details <- species_vs_rest_GO_export %>%
  dplyr::filter(
    (species == "Ab" &
       group == "Other species higher" &
       term_name %in% c(
         "translation",
         "ribosome",
         "structural constituent of ribosome",
         "eukaryotic translation initiation factor 3 complex"
       )) |
      
      (species == "On" &
         group == "On higher" &
         term_name %in% c(
           "translation",
           "eukaryotic translation initiation factor 3 complex"
         )) |
      
      (species == "On" &
         group == "Other species higher" &
         term_name == "DNA endonuclease activity") |
      
      (species == "Pn" &
         group == "Other species higher" &
         term_name == "catalytic activity, acting on DNA")
  ) %>%
  dplyr::select(
    species,
    group,
    term_id,
    term_name,
    source,
    p_value,
    fold_enrichment,
    intersection_size,
    intersection
  )

figure3_gene_details


## split GO-term intersections into one gene per row
figure3_genes_long <- figure3_gene_details %>%
  tidyr::separate_rows(
    intersection,
    sep = ","
  ) %>%
  dplyr::mutate(
    intersection = trimws(intersection)
  ) %>%
  dplyr::rename(
    gene_name = intersection
  )



## TRANSLATION

## Ab-lower translation genes
ab_translation <- figure3_genes_long %>%
  dplyr::filter(
    species == "Ab",
    group == "Other species higher",
    term_name == "translation"
  ) %>%
  dplyr::pull(gene_name) %>%
  unique()


## On-higher translation genes
on_translation <- figure3_genes_long %>%
  dplyr::filter(
    species == "On",
    group == "On higher",
    term_name == "translation"
  ) %>%
  dplyr::pull(gene_name) %>%
  unique()


## shared genes
translation_shared <- intersect(
  ab_translation,
  on_translation
)


## genes unique to Ab-lower enrichment
translation_ab_only <- setdiff(
  ab_translation,
  on_translation
)


## genes unique to On-higher enrichment
translation_on_only <- setdiff(
  on_translation,
  ab_translation
)


## union of both gene sets
translation_union <- union(
  ab_translation,
  on_translation
)


## overlap statistics
translation_jaccard <-
  length(translation_shared) /
  length(translation_union)

translation_percent_ab_shared <-
  100 *
  length(translation_shared) /
  length(ab_translation)

translation_percent_on_shared <-
  100 *
  length(translation_shared) /
  length(on_translation)


## inspect translation results
length(ab_translation)
length(on_translation)
length(translation_shared)
length(translation_ab_only)
length(translation_on_only)

translation_shared
translation_ab_only
translation_on_only

translation_jaccard
translation_percent_ab_shared
translation_percent_on_shared



## EIF3 COMPLEX

## Ab-lower eIF3 complex genes
ab_eif3 <- figure3_genes_long %>%
  dplyr::filter(
    species == "Ab",
    group == "Other species higher",
    term_name ==
      "eukaryotic translation initiation factor 3 complex"
  ) %>%
  dplyr::pull(gene_name) %>%
  unique()


## On-higher eIF3-complex genes
on_eif3 <- figure3_genes_long %>%
  dplyr::filter(
    species == "On",
    group == "On higher",
    term_name ==
      "eukaryotic translation initiation factor 3 complex"
  ) %>%
  dplyr::pull(gene_name) %>%
  unique()


## shared genes
eif3_shared <- intersect(
  ab_eif3,
  on_eif3
)


## genes unique to Ab-lower enrichment
eif3_ab_only <- setdiff(
  ab_eif3,
  on_eif3
)


## genes unique to On-higher enrichment
eif3_on_only <- setdiff(
  on_eif3,
  ab_eif3
)


## union of both gene sets
eif3_union <- union(
  ab_eif3,
  on_eif3
)


## overlap statistiscs
eif3_jaccard <-
  length(eif3_shared) /
  length(eif3_union)

eif3_percent_ab_shared <-
  100 *
  length(eif3_shared) /
  length(ab_eif3)

eif3_percent_on_shared <-
  100 *
  length(eif3_shared) /
  length(on_eif3)


## inspect eIF3 results
length(ab_eif3)
length(on_eif3)
length(eif3_shared)
length(eif3_ab_only)
length(eif3_on_only)

eif3_shared
eif3_ab_only
eif3_on_only

eif3_jaccard
eif3_percent_ab_shared
eif3_percent_on_shared




## SUMMARY TABLE

go_overlap_summary <- tibble::tibble(
  GO_term = c(
    "translation",
    "eukaryotic translation initiation factor 3 complex"
  ),
  
  Ab_lower_n = c(
    length(ab_translation),
    length(ab_eif3)
  ),
  
  On_higher_n = c(
    length(on_translation),
    length(on_eif3)
  ),
  
  shared_n = c(
    length(translation_shared),
    length(eif3_shared)
  ),
  
  Ab_only_n = c(
    length(translation_ab_only),
    length(eif3_ab_only)
  ),
  
  On_only_n = c(
    length(translation_on_only),
    length(eif3_on_only)
  ),
  
  percent_Ab_shared = c(
    translation_percent_ab_shared,
    eif3_percent_ab_shared
  ),
  
  percent_On_shared = c(
    translation_percent_on_shared,
    eif3_percent_on_shared
  ),
  
  Jaccard = c(
    translation_jaccard,
    eif3_jaccard
  )
)

go_overlap_summary



## Compare shrunken LFCs across species for shared translation genes

## COMBINE SPECIES-VS-REST DE RESULTS

de_all_species <- purrr::map_dfr(
  names(species_vs_rest_results),
  function(sp) {
    
    species_vs_rest_results[[sp]] %>%
      dplyr::mutate(
        species = sp
      )
  }
)


## check available columns
names(de_all_species)



## SHARED TRANSLATION GENES

## extract LFC and significance across all five species
translation_lfc_long <- de_all_species %>%
  dplyr::filter(
    gene_name %in% translation_shared
  ) %>%
  dplyr::select(
    species,
    gene_name,
    log2FoldChange,
    padj
  )


## make LFC table with one row per gene
translation_lfc_wide <- translation_lfc_long %>%
  dplyr::select(
    species,
    gene_name,
    log2FoldChange
  ) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = log2FoldChange
  )


## make significance table
translation_sig_wide <- translation_lfc_long %>%
  dplyr::mutate(
    significant =
      !is.na(padj) &
      padj < 0.05
  ) %>%
  dplyr::select(
    species,
    gene_name,
    significant
  ) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = significant,
    names_prefix = "sig_"
  )


## combine LFC and significance information
translation_summary <- translation_lfc_wide %>%
  dplyr::left_join(
    translation_sig_wide,
    by = "gene_name"
  )


## identify whether Ab and On occupy opposite extremes
translation_extremes <- translation_lfc_wide %>%
  dplyr::mutate(
    
    Ab_is_lowest =
      Ab == pmin(
        Ab,
        Mz,
        Nb,
        On,
        Pn,
        na.rm = TRUE
      ),
    
    On_is_highest =
      On == pmax(
        Ab,
        Mz,
        Nb,
        On,
        Pn,
        na.rm = TRUE
      ),
    
    opposite_extremes =
      Ab_is_lowest &
      On_is_highest
  )


## summarise translation pattern
translation_extreme_summary <- translation_extremes %>%
  dplyr::summarise(
    
    n_genes = dplyr::n(),
    
    Ab_lowest_n =
      sum(
        Ab_is_lowest,
        na.rm = TRUE
      ),
    
    On_highest_n =
      sum(
        On_is_highest,
        na.rm = TRUE
      ),
    
    opposite_extremes_n =
      sum(
        opposite_extremes,
        na.rm = TRUE
      ),
    
    opposite_extremes_percent =
      100 *
      mean(
        opposite_extremes,
        na.rm = TRUE
      )
  )



## SHARED EIF3-COMPLEX GENES

## extract LFC and significance across all five species
eif3_lfc_long <- de_all_species %>%
  dplyr::filter(
    gene_name %in% eif3_shared
  ) %>%
  dplyr::select(
    species,
    gene_name,
    log2FoldChange,
    padj
  )


## make LFC table with one row per gene
eif3_lfc_wide <- eif3_lfc_long %>%
  dplyr::select(
    species,
    gene_name,
    log2FoldChange
  ) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = log2FoldChange
  )


## make significance table
eif3_sig_wide <- eif3_lfc_long %>%
  dplyr::mutate(
    significant =
      !is.na(padj) &
      padj < 0.05
  ) %>%
  dplyr::select(
    species,
    gene_name,
    significant
  ) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = significant,
    names_prefix = "sig_"
  )


## combine LFC and significance information
eif3_summary <- eif3_lfc_wide %>%
  dplyr::left_join(
    eif3_sig_wide,
    by = "gene_name"
  )


## identify whether Ab and On occupy opposite extremes
eif3_extremes <- eif3_lfc_wide %>%
  dplyr::mutate(
    
    Ab_is_lowest =
      Ab == pmin(
        Ab,
        Mz,
        Nb,
        On,
        Pn,
        na.rm = TRUE
      ),
    
    On_is_highest =
      On == pmax(
        Ab,
        Mz,
        Nb,
        On,
        Pn,
        na.rm = TRUE
      ),
    
    opposite_extremes =
      Ab_is_lowest &
      On_is_highest
  )


## summarise eIF3 pattern
eif3_extreme_summary <- eif3_extremes %>%
  dplyr::summarise(
    
    n_genes = dplyr::n(),
    
    Ab_lowest_n =
      sum(
        Ab_is_lowest,
        na.rm = TRUE
      ),
    
    On_highest_n =
      sum(
        On_is_highest,
        na.rm = TRUE
      ),
    
    opposite_extremes_n =
      sum(
        opposite_extremes,
        na.rm = TRUE
      ),
    
    opposite_extremes_percent =
      100 *
      mean(
        opposite_extremes,
        na.rm = TRUE
      )
  )


## PRINT RESULTS

cat("\n AVAILABLE DE COLUMNS: \n")
print(names(de_all_species))

cat("\n TRANSLATION LFCs: \n")
print(
  translation_lfc_wide,
  n = Inf
)

cat("\n TRANSLATION EXTREME SUMMARY: \n")
print(
  translation_extreme_summary
)

cat("\n EIF3 LFCs: \n")
print(
  eif3_lfc_wide,
  n = Inf
)

cat("\n EIF3 EXTREME SUMMARY: \n")
print(
  eif3_extreme_summary
)



## Additional GO follow up?

## ATP METABOLIC PROCESS IN Ab

atp_details <- species_vs_rest_GO_export %>%
  dplyr::filter(
    species == "Ab",
    group == "Other species higher",
    term_name == "ATP metabolic process"
  ) %>%
  dplyr::select(
    species,
    group,
    term_id,
    term_name,
    source,
    p_value,
    fold_enrichment,
    intersection_size,
    intersection
  )

cat("\n ATP METABOLIC PROCESS DETAILS: \n")
print(atp_details)


## extract individual genes
atp_genes <- atp_details %>%
  tidyr::separate_rows(
    intersection,
    sep = ","
  ) %>%
  dplyr::mutate(
    intersection = trimws(intersection)
  ) %>%
  dplyr::pull(intersection) %>%
  unique()

cat("\n ATP METABOLIC PROCESS GENES: \n")
print(atp_genes)

cat("\nNumber of ATP-related genes:\n")
print(length(atp_genes))


## examine these genes across all five species
atp_lfc_long <- de_all_species %>%
  dplyr::filter(
    gene_name %in% atp_genes
  ) %>%
  dplyr::select(
    species,
    gene_name,
    log2FoldChange,
    padj
  )


## LFC table
atp_lfc_wide <- atp_lfc_long %>%
  dplyr::select(
    species,
    gene_name,
    log2FoldChange
  ) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = log2FoldChange
  )


## significance counts across species
atp_significance_counts <- atp_lfc_long %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    n_significant = sum(
      !is.na(padj) & padj < 0.05
    ),
    percent_significant =
      100 * n_significant / n_genes,
    .groups = "drop"
  )


cat("\n ATP GENES: LFCs ACROSS SPECIES: \n")
print(
  atp_lfc_wide,
  n = Inf
)

cat("\n ATP GENES: SIGNIFICANCE ACROSS SPECIES: \n")
print(atp_significance_counts)


## DNA-RELATED TERMS IN On AND Pn

dna_details <- species_vs_rest_GO_export %>%
  dplyr::filter(
    (
      species == "On" &
        group == "Other species higher" &
        term_name == "DNA endonuclease activity"
    ) |
      (
        species == "Pn" &
          group == "Other species higher" &
          term_name == "catalytic activity, acting on DNA"
      )
  ) %>%
  dplyr::select(
    species,
    group,
    term_id,
    term_name,
    source,
    p_value,
    fold_enrichment,
    intersection_size,
    intersection
  )

cat("\n DNA-RELATED GO DETAILS: \n")
print(dna_details)


## On-lower DNA endonuclease genes
on_dna <- dna_details %>%
  dplyr::filter(
    species == "On",
    term_name == "DNA endonuclease activity"
  ) %>%
  tidyr::separate_rows(
    intersection,
    sep = ","
  ) %>%
  dplyr::mutate(
    intersection = trimws(intersection)
  ) %>%
  dplyr::pull(intersection) %>%
  unique()


## Pn-lower catalytic activity acting on DNA genes
pn_dna <- dna_details %>%
  dplyr::filter(
    species == "Pn",
    term_name == "catalytic activity, acting on DNA"
  ) %>%
  tidyr::separate_rows(
    intersection,
    sep = ","
  ) %>%
  dplyr::mutate(
    intersection = trimws(intersection)
  ) %>%
  dplyr::pull(intersection) %>%
  unique()


## compare gene sets
dna_shared <- intersect(
  on_dna,
  pn_dna
)

dna_on_only <- setdiff(
  on_dna,
  pn_dna
)

dna_pn_only <- setdiff(
  pn_dna,
  on_dna
)

dna_union <- union(
  on_dna,
  pn_dna
)


## overlap statistics
dna_jaccard <-
  ifelse(
    length(dna_union) > 0,
    length(dna_shared) /
      length(dna_union),
    NA_real_
  )

dna_percent_on_shared <-
  ifelse(
    length(on_dna) > 0,
    100 *
      length(dna_shared) /
      length(on_dna),
    NA_real_
  )

dna_percent_pn_shared <-
  ifelse(
    length(pn_dna) > 0,
    100 *
      length(dna_shared) /
      length(pn_dna),
    NA_real_
  )


## summary table
dna_overlap_summary <- tibble::tibble(
  On_lower_n = length(on_dna),
  Pn_lower_n = length(pn_dna),
  shared_n = length(dna_shared),
  On_only_n = length(dna_on_only),
  Pn_only_n = length(dna_pn_only),
  percent_On_shared = dna_percent_on_shared,
  percent_Pn_shared = dna_percent_pn_shared,
  Jaccard = dna_jaccard
)


cat("\n On DNA ENDONUCLEASE GENES: \n")
print(on_dna)

cat("\n Pn DNA-ACTIVITY GENES: \n")
print(pn_dna)

cat("\n SHARED DNA-RELATED GENES: \n")
print(dna_shared)

cat("\n On-ONLY DNA-RELATED GENES: \n")
print(dna_on_only)

cat("\n Pn-ONLY DNA-RELATED GENES: \n")
print(dna_pn_only)

cat("\n DNA OVERLAP SUMMARY: \n")
print(dna_overlap_summary)



## HEATMAP:
## Shared translation-related genes across all five species


## Get ribosome-related genes from the Ab-lower GO results

ab_ribosome <- figure3_genes_long %>%
  dplyr::filter(
    species == "Ab",
    group == "Other species higher",
    term_name == "ribosome"
  ) %>%
  dplyr::pull(gene_name) %>%
  unique()


ab_ribosome_mf <- figure3_genes_long %>%
  dplyr::filter(
    species == "Ab",
    group == "Other species higher",
    term_name == "structural constituent of ribosome"
  ) %>%
  dplyr::pull(gene_name) %>%
  unique()


ribosome_genes <- union(
  ab_ribosome,
  ab_ribosome_mf
)


## Prep heatmap data

translation_heatmap_data <- translation_lfc_long %>%
  dplyr::mutate(
    
    ## fixed species order
    species = factor(
      species,
      levels = c("Ab", "Mz", "Nb", "On", "Pn")
    ),
    
    ## DE significance
    significant =
      !is.na(padj) & padj < 0.05,
    
    ## assign every gene to ONE functional category
    ## eIF3 takes priority over ribosome membership
    functional_group = dplyr::case_when(
      
      gene_name %in% eif3_shared ~
        "eIF3-related",
      
      gene_name %in% ribosome_genes ~
        "Ribosome-related",
      
      TRUE ~
        "Other translation-related"
    ),
    
    functional_group = factor(
      functional_group,
      levels = c(
        "eIF3-related",
        "Ribosome-related",
        "Other translation-related"
      )
    ),
    
    ## cap only for colour display
    ## original LFC values are unchanged
    lfc_plot = pmax(
      pmin(log2FoldChange, 8),
      -8
    )
  )


## Check functional assignments

translation_gene_groups <- translation_heatmap_data %>%
  dplyr::distinct(
    gene_name,
    functional_group
  ) %>%
  dplyr::arrange(
    functional_group,
    gene_name
  )


cat("\n FUNCTIONAL GROUP ASSIGNMENTS: \n")

print(
  translation_gene_groups,
  n = Inf
)


cat("\n NUMBER OF GENES PER GROUP: \n")

print(
  table(
    translation_gene_groups$functional_group
  )
)


cat("\n TOTAL UNIQUE GENES: \n")

print(
  dplyr::n_distinct(
    translation_gene_groups$gene_name
  )
)


## This should equal 32
stopifnot(
  dplyr::n_distinct(
    translation_gene_groups$gene_name
  ) == 32
)


## cluster genes within each functional group

cluster_gene_group <- function(group_name) {
  
  group_matrix <- translation_heatmap_data %>%
    
    dplyr::filter(
      functional_group == group_name
    ) %>%
    
    dplyr::select(
      gene_name,
      species,
      log2FoldChange
    ) %>%
    
    tidyr::pivot_wider(
      names_from = species,
      values_from = log2FoldChange
    ) %>%
    
    tibble::column_to_rownames(
      "gene_name"
    ) %>%
    
    as.matrix()
  
  
  ## enforce species order
  group_matrix <- group_matrix[
    ,
    c("Ab", "Mz", "Nb", "On", "Pn"),
    drop = FALSE
  ]
  
  
  ## cluster if group contains >1 gene
  if (nrow(group_matrix) > 1) {
    
    hc <- hclust(
      dist(
        group_matrix,
        method = "euclidean"
      ),
      method = "complete"
    )
    
    rownames(group_matrix)[hc$order]
    
  } else {
    
    rownames(group_matrix)
  }
}


## obtain gene order separately for each category
eif3_order <- cluster_gene_group(
  "eIF3-related"
)

ribosome_order <- cluster_gene_group(
  "Ribosome-related"
)

other_translation_order <- cluster_gene_group(
  "Other translation-related"
)


## combined ordering
translation_gene_order <- c(
  eif3_order,
  ribosome_order,
  other_translation_order
)




## Set factor order for plotting
translation_heatmap_data <- translation_heatmap_data %>%
  dplyr::mutate(
    
    gene_name = factor(
      gene_name,
      levels = rev(
        translation_gene_order
      )
    )
  )


## Generate heatmap

translation_heatmap <- ggplot(
  translation_heatmap_data,
  aes(
    x = species,
    y = gene_name,
    fill = lfc_plot
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +
  
  ## open circle = significant DE
  geom_point(
    data = translation_heatmap_data %>%
      dplyr::filter(
        significant
      ),
    aes(
      x = species,
      y = gene_name
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = NA,
    colour = "black",
    stroke = 0.65,
    size = 2.2
  ) +
  
  ## separate functional groups
  facet_grid(
    functional_group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  ## same interpretation as other figures - colour coding
  scale_fill_gradient2(
    low = "#D95F02",
    mid = "grey95",
    high = "#1F78B4",
    midpoint = 0,
    limits = c(-8, 8),
    breaks = c(-8, -4, 0, 4, 8),
    name = expression(
      "Log"[2] * " fold change"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    
    ## gene symbols italicised
    axis.text.y = element_text(
      size = 8,
      face = "italic"
    ),
    
    axis.ticks = element_blank(),
    
    ## put functional group labels on left
    strip.placement = "outside",
    
    strip.background = element_blank(),
    
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 9
    ),
    
    panel.spacing.y = grid::unit(
      0.45,
      "lines"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 9
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      5.5,
      5.5,
      5.5,
      5.5
    )
  )


## display
translation_heatmap


dir.create(
  "Figures/GO/species_vs_rest",
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename =
    "Figures/GO/species_vs_rest/Figure4_translation_shared_genes_grouped.png",
  plot = translation_heatmap,
  width = 7.2,
  height = 8.5,
  units = "in",
  dpi = 600
)

