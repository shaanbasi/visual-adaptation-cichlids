library(dplyr)
library(purrr)
library(gprofiler2)
library(ggplot2)
library(stringr)


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
  "GO_all_pairwise_clean.csv",
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
  "GO_vision_all_pairwise.csv",
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