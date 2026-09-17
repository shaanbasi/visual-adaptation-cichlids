library(dplyr)
library(readr)
library(tibble)
library(tidyr)

options(stringsAsFactors = FALSE)


## INPUT/OUTPUT
visual_dir <- file.path(
  "Results",
  "targeted_visual_genes"
)

wgcna_dir <- file.path(
  "Results",
  "WGCNA"
)

output_dir <- file.path(
  visual_dir,
  "WGCNA_integration"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



## LOAD TARGETED VISUAL-GENE DE RESULTS

visual_DE <- readr::read_csv(
  file.path(
    visual_dir,
    "curated_visual_genes_all_species_vs_rest.csv"
  ),
  show_col_types = FALSE
)


## standardise gene name

visual_DE <- visual_DE %>%
  dplyr::mutate(
    
    gene_name =
      dplyr::coalesce(
        gene_name.y,
        gene_name.x
      ),
    
    gene_match =
      tolower(
        trimws(
          gene_match
        )
      ),
    
    significant =
      !is.na(padj) &
      padj < 0.05
  )




## LOAD FINAL WGCNA NETWORK OBJECT

load(
  file.path(
    wgcna_dir,
    "WGCNA_network_objects.RData"
  )
)


required_objects <- c(
  "geneModuleAssignments",
  "moduleSpeciesResults"
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
    "Missing WGCNA objects: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}



## GENE TO WGCNA MODULE LOOKUP

visual_module_lookup <- geneModuleAssignments %>%
  
  dplyr::transmute(
    
    dataset_gene_name =
      as.character(
        Gene
      ),
    
    gene_match =
      tolower(
        trimws(
          dataset_gene_name
        )
      ),
    
    WGCNA_module =
      sub(
        "^ME",
        "",
        as.character(
          Module
        )
      )
  ) %>%
  
  dplyr::distinct(
    gene_match,
    .keep_all = TRUE
  )



## MODULE TO SPECIES ASSOCIATION LOOKUP

module_species_lookup <- moduleSpeciesResults %>%
  
  dplyr::transmute(
    
    species =
      as.character(
        species
      ),
    
    WGCNA_module =
      as.character(
        module
      ),
    
    module_trait_correlation =
      correlation,
    
    module_trait_FDR =
      FDR,
    
    module_expression_pattern =
      expression_pattern,
    
    species_associated_module =
      !is.na(FDR) &
      FDR < 0.05
  )



## INTEGRATE DE AND WGCNA
visual_DE_WGCNA <- visual_DE %>%
  
  dplyr::left_join(
    visual_module_lookup,
    by = "gene_match"
  ) %>%
  
  dplyr::left_join(
    module_species_lookup,
    by = c(
      "species",
      "WGCNA_module"
    )
  ) %>%
  
  dplyr::mutate(
    
    ## explicit evidence categories
    evidence_class =
      dplyr::case_when(
        
        significant &
          species_associated_module ~
          "DE + species-associated module",
        
        significant ~
          "DE only",
        
        species_associated_module ~
          "Species-associated module only",
        
        TRUE ~
          "No significant species-associated evidence"
      )
  )


readr::write_csv(
  visual_DE_WGCNA,
  file.path(
    output_dir,
    "visual_gene_DE_WGCNA_integration.csv"
  )
)



## SIGNIFICANT VISUAL DE + NETWORK CONTEXT

significant_visual_integrated <- visual_DE_WGCNA %>%
  
  dplyr::filter(
    significant
  ) %>%
  
  dplyr::arrange(
    species,
    dplyr::desc(
      species_associated_module
    ),
    padj
  ) %>%
  
  dplyr::select(
    species,
    gene_name,
    visual_category,
    log2FoldChange,
    padj,
    direction,
    WGCNA_module,
    module_trait_correlation,
    module_trait_FDR,
    module_expression_pattern,
    species_associated_module,
    evidence_class
  )


print(
  significant_visual_integrated,
  n = Inf,
  width = Inf
)


readr::write_csv(
  significant_visual_integrated,
  file.path(
    output_dir,
    "significant_visual_genes_with_WGCNA_context.csv"
  )
)



## VISUAL GENES WITH BOTH DE AND MODULE SUPPORT

visual_DE_module_overlap <- visual_DE_WGCNA %>%
  
  dplyr::filter(
    significant,
    species_associated_module
  ) %>%
  
  dplyr::arrange(
    species,
    module_trait_FDR,
    padj
  )


print(
  visual_DE_module_overlap,
  n = Inf,
  width = Inf
)


readr::write_csv(
  visual_DE_module_overlap,
  file.path(
    output_dir,
    "visual_genes_DE_and_species_module_support.csv"
  )
)


##SUMMARY

visual_evidence_summary <- visual_DE_WGCNA %>%
  
  dplyr::filter(
    significant
  ) %>%
  
  dplyr::count(
    species,
    evidence_class,
    name = "n_visual_gene_results"
  )


print(
  visual_evidence_summary,
  n = Inf
)


readr::write_csv(
  visual_evidence_summary,
  file.path(
    output_dir,
    "visual_gene_evidence_summary.csv"
  )
)




## GENELEVEL PRIORTISATION

visual_gene_prioritisation <- visual_DE_WGCNA %>%
  
  dplyr::group_by(
    gene_match,
    gene_name,
    visual_category
  ) %>%
  
  dplyr::summarise(
    
    n_species_DE =
      sum(
        significant,
        na.rm = TRUE
      ),
    
    n_species_module_supported =
      sum(
        significant &
          species_associated_module,
        na.rm = TRUE
      ),
    
    strongest_abs_log2FC =
      max(
        abs(
          log2FoldChange
        ),
        na.rm = TRUE
      ),
    
    best_DE_FDR =
      if (
        all(
          is.na(
            padj
          )
        )
      ) {
        NA_real_
      } else {
        min(
          padj,
          na.rm = TRUE
        )
      },
    
    DE_patterns =
      paste(
        direction[
          significant
        ],
        collapse = "; "
      ),
    
    modules =
      paste(
        unique(
          WGCNA_module[
            !is.na(
              WGCNA_module
            )
          ]
        ),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      n_species_module_supported
    ),
    dplyr::desc(
      n_species_DE
    ),
    best_DE_FDR
  )


print(
  visual_gene_prioritisation,
  n = Inf,
  width = Inf
)


readr::write_csv(
  visual_gene_prioritisation,
  file.path(
    output_dir,
    "visual_gene_prioritisation.csv"
  )
)






## ADD GENE-LEVEL WGCNA MM / GS / HUB SUPPORT

hub_dir <- file.path(
  wgcna_dir,
  "hub_genes"
)


## complete gene-level MM and GS metrics

all_gene_WGCNA_metrics <- readr::read_csv(
  file.path(
    hub_dir,
    "05_all_gene_MM_and_GS_metrics.csv"
  ),
  show_col_types = FALSE
)


## genes passing the final supported-hub criteria

supported_hubs <- readr::read_csv(
  file.path(
    hub_dir,
    "08_supported_hub_genes.csv"
  ),
  show_col_types = FALSE
)


## inspect

print(
  names(
    all_gene_WGCNA_metrics
  )
)

print(
  names(
    supported_hubs
  )
)




## EXTRACT OWN-MODULE MM FOR EACH GENE

MM_long <- all_gene_WGCNA_metrics %>%
  
  dplyr::select(
    gene_id,
    module,
    dplyr::starts_with("MM.")
  ) %>%
  
  tidyr::pivot_longer(
    cols = dplyr::starts_with("MM."),
    names_to = "MM_module",
    values_to = "MM"
  ) %>%
  
  dplyr::mutate(
    MM_module = sub(
      "^MM\\.",
      "",
      MM_module
    )
  ) %>%
  
  dplyr::filter(
    module == MM_module
  ) %>%
  
  dplyr::transmute(
    gene_match =
      tolower(
        trimws(
          gene_id
        )
      ),
    
    WGCNA_module =
      module,
    
    MM
  )



## EXTRACT OWN MODULE MM P-VALUE

MM_p_long <- all_gene_WGCNA_metrics %>%
  
  dplyr::select(
    gene_id,
    module,
    dplyr::starts_with("p.MM.")
  ) %>%
  
  tidyr::pivot_longer(
    cols = dplyr::starts_with("p.MM."),
    names_to = "MM_module",
    values_to = "MM_pvalue"
  ) %>%
  
  dplyr::mutate(
    MM_module = sub(
      "^p\\.MM\\.",
      "",
      MM_module
    )
  ) %>%
  
  dplyr::filter(
    module == MM_module
  ) %>%
  
  dplyr::transmute(
    gene_match =
      tolower(
        trimws(
          gene_id
        )
      ),
    
    WGCNA_module =
      module,
    
    MM_pvalue
  )


## EXTRACT SPECIES-SPECIFIC GS

GS_long <- all_gene_WGCNA_metrics %>%
  
  dplyr::select(
    gene_id,
    dplyr::starts_with("GS.")
  ) %>%
  
  tidyr::pivot_longer(
    cols = dplyr::starts_with("GS."),
    names_to = "species",
    values_to = "GS"
  ) %>%
  
  dplyr::mutate(
    species = sub(
      "^GS\\.",
      "",
      species
    ),
    
    gene_match =
      tolower(
        trimws(
          gene_id
        )
      )
  ) %>%
  
  dplyr::select(
    gene_match,
    species,
    GS
  )



## EXTRACT SPECIES-SPECIFIC GS P-VALUES

GS_p_long <- all_gene_WGCNA_metrics %>%
  
  dplyr::select(
    gene_id,
    dplyr::starts_with("p.GS.")
  ) %>%
  
  tidyr::pivot_longer(
    cols = dplyr::starts_with("p.GS."),
    names_to = "species",
    values_to = "GS_pvalue"
  ) %>%
  
  dplyr::mutate(
    
    species = sub(
      "^p\\.GS\\.",
      "",
      species
    ),
    
    gene_match =
      tolower(
        trimws(
          gene_id
        )
      )
  ) %>%
  
  dplyr::select(
    gene_match,
    species,
    GS_pvalue
  )



## SUPPORTED HUB LOOKUP

supported_hub_lookup <- supported_hubs %>%
  
  dplyr::mutate(
    gene_match =
      tolower(
        trimws(
          gene_id
        )
      ),
    
    supported_hub = TRUE
  ) %>%
  
  dplyr::select(
    gene_match,
    species,
    module,
    supported_hub,
    passes_MM,
    passes_aligned_GS,
    passes_GS_pvalue,
    passes_all_support_criteria,
    exploratory_joint_score,
    supported_hub_rank
  )



## FINAL VISUAL GENE DE + WGCNA + HUB INTEGRATION

visual_gene_final <- visual_DE_WGCNA %>%
  
  dplyr::left_join(
    MM_long,
    by = c(
      "gene_match",
      "WGCNA_module"
    )
  ) %>%
  
  dplyr::left_join(
    MM_p_long,
    by = c(
      "gene_match",
      "WGCNA_module"
    )
  ) %>%
  
  dplyr::left_join(
    GS_long,
    by = c(
      "gene_match",
      "species"
    )
  ) %>%
  
  dplyr::left_join(
    GS_p_long,
    by = c(
      "gene_match",
      "species"
    )
  ) %>%
  
  dplyr::left_join(
    supported_hub_lookup,
    by = c(
      "gene_match",
      "species",
      "WGCNA_module" = "module"
    )
  ) %>%
  
  dplyr::mutate(
    
    supported_hub =
      dplyr::coalesce(
        supported_hub,
        FALSE
      ),
    
    evidence_class =
      dplyr::case_when(
        
        significant &
          species_associated_module &
          supported_hub ~
          "DE + module + supported hub",
        
        significant &
          species_associated_module ~
          "DE + species-associated module",
        
        significant &
          supported_hub ~
          "DE + supported hub",
        
        significant ~
          "DE only",
        
        species_associated_module &
          supported_hub ~
          "Module + supported hub",
        
        species_associated_module ~
          "Species-associated module only",
        
        supported_hub ~
          "Supported hub only",
        
        TRUE ~
          "No significant integrated evidence"
      )
  )

readr::write_csv(
  visual_gene_final,
  file.path(
    output_dir,
    "visual_gene_final_DE_WGCNA_hub_integration.csv"
  )
)


## FINAL SIGNIFICANT VISUAL GENE CANDIDATES

final_visual_candidates <- visual_gene_final %>%
  
  dplyr::filter(
    significant
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      supported_hub
    ),
    dplyr::desc(
      species_associated_module
    ),
    species,
    padj
  ) %>%
  
  dplyr::select(
    species,
    gene_name,
    visual_category,
    
    log2FoldChange,
    padj,
    direction,
    
    WGCNA_module,
    
    module_trait_correlation,
    module_trait_FDR,
    
    MM,
    MM_pvalue,
    
    GS,
    GS_pvalue,
    
    species_associated_module,
    supported_hub,
    
    evidence_class
  )


print(
  final_visual_candidates,
  n = Inf,
  width = Inf
)


readr::write_csv(
  final_visual_candidates,
  file.path(
    output_dir,
    "final_visual_candidate_evidence.csv"
  )
)



## MULTI-LAYER SUPPORTED VISUAL CANDIDATES
multi_layer_visual_candidates <- final_visual_candidates %>%
  
  dplyr::filter(
    species_associated_module |
      supported_hub
  ) %>%
  
  dplyr::arrange(
    species,
    dplyr::desc(
      supported_hub
    ),
    module_trait_FDR,
    padj
  )


print(
  multi_layer_visual_candidates,
  n = Inf,
  width = Inf
)


readr::write_csv(
  multi_layer_visual_candidates,
  file.path(
    output_dir,
    "multi_layer_supported_visual_candidates.csv"
  )
)



