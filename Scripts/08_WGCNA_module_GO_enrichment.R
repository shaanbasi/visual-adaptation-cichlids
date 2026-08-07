library(WGCNA)
library(dplyr)
library(purrr)
library(gprofiler2)
library(readr)
library(tibble)

options(stringsAsFactors = FALSE)


## INPUT AND OUTPUT SETTINGS

wgcna_dir <- "Results/WGCNA"

hub_dir <- file.path(
  wgcna_dir,
  "hub_genes"
)

output_dir <- file.path(
  wgcna_dir,
  "module_GO"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}



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
  "datExpr",
  "geneModuleAssignments",
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



## STANDARDISE GENE MODULE TABLE

geneModuleTable <- geneModuleAssignments %>%
  dplyr::transmute(
    gene_id = as.character(Gene),
    module = sub(
      "^ME",
      "",
      as.character(Module)
    )
  )

stopifnot(
  setequal(
    geneModuleTable$gene_id,
    colnames(datExpr)
  )
)



## RNA-SEQ BACKGROUND
## all genes retained in the WGCNA network

background_genes <- geneModuleTable %>%
  dplyr::filter(
    !is.na(gene_id),
    trimws(gene_id) != ""
  ) %>%
  dplyr::pull(gene_id) %>%
  unique()



## RUN ONE MODULE GO ANALYSIS

run_module_GO <- function(
    module_name,
    species_name,
    background_genes
) {
  
  module_genes <- geneModuleTable %>%
    dplyr::filter(
      module == module_name,
      !is.na(gene_id),
      trimws(gene_id) != ""
    ) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  module_genes <- intersect(
    module_genes,
    background_genes
  )
  
  message(
    "Running GO enrichment for ",
    species_name,
    " | ",
    module_name,
    " module | ",
    length(module_genes),
    " genes"
  )
  
  if (length(module_genes) == 0) {
    
    message(
      "No genes found for ",
      species_name,
      " | ",
      module_name
    )
    
    return(tibble())
  }
  
  gost_result <- tryCatch(
    
    gost(
      query = module_genes,
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
        species_name,
        " | ",
        module_name,
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
      species_name,
      " | ",
      module_name
    )
    
    return(tibble())
  }
  
  gost_result$result %>%
    dplyr::mutate(
      species = species_name,
      module = module_name,
      module_gene_count = length(module_genes),
      background_gene_count = length(background_genes)
    )
}



## RUN GO ENRICHMENT FOR SELECTED MODULE-SPECIES PAIRS

module_GO_results <- purrr::map_dfr(
  seq_len(
    nrow(selectedModuleSpeciesPairs)
  ),
  function(i) {
    
    run_module_GO(
      module_name =
        selectedModuleSpeciesPairs$module[i],
      
      species_name =
        selectedModuleSpeciesPairs$species[i],
      
      background_genes =
        background_genes
    )
  }
)



## KEEP GO RESULTS
module_GO_clean <- module_GO_results %>%
  dplyr::filter(
    source %in% c(
      "GO:BP",
      "GO:MF",
      "GO:CC"
    )
  ) %>%
  dplyr::arrange(
    species,
    module,
    p_value
  )


## EXPORT ALL SIGNIFICANT MODULE GO RESULTS

module_GO_export <- module_GO_clean %>%
  dplyr::select(
    species,
    module,
    term_id,
    term_name,
    source,
    p_value,
    intersection_size,
    query_size,
    term_size,
    effective_domain_size,
    module_gene_count,
    background_gene_count,
    dplyr::everything()
  )

readr::write_csv(
  module_GO_export,
  file.path(
    output_dir,
    "WGCNA_selected_modules_GO_all.csv"
  )
)




## VISION- AND PIGMENT-RELATED TERMS

vision_terms <- c(
  "visual perception",
  "phototransduction",
  "light stimulus",
  "response to light",
  "opsin",
  "rhodopsin",
  "retina",
  "retinal",
  "neural retina",
  "photoreceptor",
  "camera-type eye",
  "eye development",
  "eye morphogenesis",
  "lens development",
  "pigment",
  "pigmentation",
  "melanin",
  "melanogenesis",
  "melanocyte",
  "chromatophore",
  "retinoid",
  "retinol",
  "visual cycle"
)

vision_pattern <- paste(
  vision_terms,
  collapse = "|"
)

module_vision_GO <- module_GO_clean %>%
  dplyr::filter(
    grepl(
      vision_pattern,
      term_name,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::arrange(
    species,
    module,
    p_value
  )



## EXPORT VISION-RELATED MODULE GO RESULTS

readr::write_csv(
  module_vision_GO,
  file.path(
    output_dir,
    "WGCNA_selected_modules_GO_visual_terms.csv"
  )
)



## SUMMARISE GO TERMS BY MODULE

module_GO_summary <- module_GO_clean %>%
  dplyr::count(
    species,
    module,
    source,
    name = "n_significant_GO_terms"
  ) %>%
  dplyr::arrange(
    species,
    module,
    source
  )

readr::write_csv(
  module_GO_summary,
  file.path(
    output_dir,
    "WGCNA_selected_modules_GO_summary.csv"
  )
)




## SUMMARISE VISION-RELATED GO TERMS

module_vision_GO_summary <- module_vision_GO %>%
  dplyr::group_by(
    species,
    module,
    term_id,
    term_name,
    source
  ) %>%
  dplyr::summarise(
    p_value = min(
      p_value,
      na.rm = TRUE
    ),
    intersection_size = max(
      intersection_size,
      na.rm = TRUE
    ),
    n_occurrences = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    species,
    module,
    p_value
  )

readr::write_csv(
  module_vision_GO_summary,
  file.path(
    output_dir,
    "WGCNA_selected_modules_visual_GO_summary.csv"
  )
)



## SEE RESULTS

print(
  tibble::as_tibble(module_GO_summary),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(module_vision_GO),
  n = Inf,
  width = Inf
)


##CHECK On-PINK ENRICHED TERMS

module_GO_clean %>%
  dplyr::select(
    species,
    module,
    term_id,
    term_name,
    source,
    p_value,
    intersection_size
  ) %>%
  dplyr::arrange(
    p_value
  ) %>%
  tibble::as_tibble() %>%
  print(
    n = Inf,
    width = Inf
  )


##any vision related?
module_vision_GO %>%
  dplyr::select(
    species,
    module,
    term_id,
    term_name,
    source,
    p_value,
    intersection_size
  ) %>%
  tibble::as_tibble() %>%
  print(
    n = Inf,
    width = Inf
  )
## no vision related in either module
## GO:0007399 – nervous system development
## not eye specific - but directly related ot neural development!


## FINAL SUMMARY

cat(
  "\n============================================\n"
)

cat(
  "WGCNA module GO enrichment complete\n"
)

cat(
  "Selected modules analysed:",
  nrow(selectedModuleSpeciesPairs),
  "\n"
)

cat(
  "Significant GO terms:",
  nrow(module_GO_clean),
  "\n"
)

cat(
  "Vision-related GO terms:",
  nrow(module_vision_GO),
  "\n"
)

cat(
  "Results saved to:",
  output_dir,
  "\n"
)

cat(
  "============================================\n\n"
)



## VISUAL-GENE OVERLAP WITH WGCNA MODULES AND HUB GENES

## candidate genes involved in opsin function, phototransduction,
## retinal development and the visual cycle
visual_genes <- toupper(c(
  "RH1",
  "RHO",
  "SWS1",
  "SWS2A",
  "SWS2B",
  "RH2A",
  "RH2AALPHA",
  "RH2ABETA",
  "RH2B",
  "LWS",
  "OPN1SW1",
  "OPN1SW2",
  "OPN1MW",
  "OPN1LW",
  "GNAT1",
  "GNAT2",
  "GNB1",
  "GNGT1",
  "GNGT2",
  "PDE6A",
  "PDE6B",
  "PDE6C",
  "CNGA1",
  "CNGA3",
  "CNGB1",
  "CNGB3",
  "ARR3",
  "SAG",
  "RCVRN",
  "RGS9",
  "GUCA1A",
  "GUCA1B",
  "RPE65",
  "LRAT",
  "ABCA4",
  "RDH5",
  "RDH8",
  "CRX",
  "PAX6",
  "NRL",
  "NR2E3",
  "OTX2",
  "SIX6"
))



## visual genes present anywhere in the selected WGCNA modules

visual_genes_in_selected_modules <- selectedModuleSpeciesPairs %>%
  dplyr::select(
    species,
    module,
    module_size,
    correlation,
    FDR
  ) %>%
  dplyr::left_join(
    geneModuleTable,
    by = "module"
  ) %>%
  dplyr::mutate(
    gene_upper = toupper(
      trimws(gene_id)
    )
  ) %>%
  dplyr::filter(
    !is.na(gene_upper),
    gene_upper != "",
    gene_upper %in% visual_genes
  ) %>%
  dplyr::select(
    species,
    module,
    gene_id,
    module_size,
    correlation,
    FDR
  ) %>%
  dplyr::arrange(
    species,
    module,
    gene_id
  )



## supported visual hub genes

visual_supported_hubs <- supportedHubs %>%
  dplyr::mutate(
    gene_upper = toupper(
      trimws(gene_id)
    )
  ) %>%
  dplyr::filter(
    !is.na(gene_upper),
    gene_upper != "",
    gene_upper %in% visual_genes
  ) %>%
  dplyr::select(
    species,
    module,
    gene_id,
    module_size,
    module_trait_correlation,
    module_trait_FDR,
    MM,
    MM_pvalue,
    GS,
    aligned_GS,
    GS_pvalue,
    passes_all_support_criteria
  ) %>%
  dplyr::arrange(
    species,
    module,
    dplyr::desc(MM)
  )


## visual genes among the top MM-ranked genes

visual_top_ranked_hubs <- topRankedHubs %>%
  dplyr::mutate(
    gene_upper = toupper(
      trimws(gene_id)
    )
  ) %>%
  dplyr::filter(
    !is.na(gene_upper),
    gene_upper != "",
    gene_upper %in% visual_genes
  ) %>%
  dplyr::select(
    species,
    module,
    gene_id,
    MM,
    MM_rank,
    GS,
    aligned_GS,
    GS_pvalue,
    module_trait_correlation,
    module_trait_FDR
  ) %>%
  dplyr::arrange(
    species,
    module,
    MM_rank
  )



## summary by module

visual_hub_overlap_summary <- selectedModuleSpeciesPairs %>%
  dplyr::select(
    species,
    module,
    module_size
  ) %>%
  dplyr::left_join(
    visual_genes_in_selected_modules %>%
      dplyr::group_by(
        species,
        module
      ) %>%
      dplyr::summarise(
        n_visual_genes_in_module =
          dplyr::n_distinct(gene_id),
        
        visual_genes_in_module = paste(
          sort(
            unique(gene_id)
          ),
          collapse = "; "
        ),
        
        .groups = "drop"
      ),
    
    by = c(
      "species",
      "module"
    )
  ) %>%
  dplyr::left_join(
    visual_supported_hubs %>%
      dplyr::group_by(
        species,
        module
      ) %>%
      dplyr::summarise(
        n_supported_visual_hubs =
          dplyr::n_distinct(gene_id),
        
        supported_visual_hubs = paste(
          sort(
            unique(gene_id)
          ),
          collapse = "; "
        ),
        
        .groups = "drop"
      ),
    
    by = c(
      "species",
      "module"
    )
  ) %>%
  dplyr::mutate(
    n_visual_genes_in_module =
      tidyr::replace_na(
        n_visual_genes_in_module,
        0L
      ),
    
    visual_genes_in_module =
      tidyr::replace_na(
        visual_genes_in_module,
        "None"
      ),
    
    n_supported_visual_hubs =
      tidyr::replace_na(
        n_supported_visual_hubs,
        0L
      ),
    
    supported_visual_hubs =
      tidyr::replace_na(
        supported_visual_hubs,
        "None"
      )
  ) %>%
  dplyr::arrange(
    species,
    module
  )



## EXPORT RESULTS

readr::write_csv(
  visual_genes_in_selected_modules,
  file.path(
    output_dir,
    "WGCNA_visual_genes_in_selected_modules.csv"
  )
)

readr::write_csv(
  visual_supported_hubs,
  file.path(
    output_dir,
    "WGCNA_supported_visual_hub_genes.csv"
  )
)

readr::write_csv(
  visual_top_ranked_hubs,
  file.path(
    output_dir,
    "WGCNA_visual_genes_among_top_ranked_hubs.csv"
  )
)

readr::write_csv(
  visual_hub_overlap_summary,
  file.path(
    output_dir,
    "WGCNA_visual_hub_overlap_summary.csv"
  )
)



## DISPLAY RESULTS

print(
  tibble::as_tibble(
    visual_hub_overlap_summary
  ),
  n = Inf,
  width = Inf
)

print(
  tibble::as_tibble(
    visual_supported_hubs
  ),
  n = Inf,
  width = Inf
)


## FINAL OVERLAP SUMMARY

cat(
  "\n============================================\n"
)

cat(
  "Visual-gene and WGCNA hub overlap complete\n"
)

cat(
  "Visual genes in selected modules:",
  nrow(
    visual_genes_in_selected_modules
  ),
  "\n"
)

cat(
  "Supported visual hub genes:",
  nrow(
    visual_supported_hubs
  ),
  "\n"
)

cat(
  "Visual genes among top-ranked hubs:",
  nrow(
    visual_top_ranked_hubs
  ),
  "\n"
)

cat(
  "Results saved to:",
  output_dir,
  "\n"
)

cat(
  "============================================\n\n"
)