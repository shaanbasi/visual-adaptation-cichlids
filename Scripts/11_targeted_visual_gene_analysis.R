library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(purrr) 

options(stringsAsFactors = FALSE)


## INPUT/ OUTPUT

deseq_dir <- "Results/DESeq2"

species_vs_rest_dir <- file.path(
  deseq_dir,
  "species_vs_rest"
)

wgcna_dir <- "Results/WGCNA"

output_dir <- file.path(
  "Results",
  "targeted_visual_genes"
)

figure_dir <- file.path(
  "Figures",
  "targeted_visual_genes"
)

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


species_levels <- c(
  "Ab",
  "Mz",
  "Nb",
  "On",
  "Pn"
)



## LOAD VISUAL-GENE REFERENCE

visual_gene_reference <- readr::read_csv(
  "Data/Reference/visual_candidate_genes.csv",
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    gene_match = tolower(
      trimws(
        gene_name
      )
    )
  )


## LOAD NORMALIZED EXPRESSION
normalized_counts <- read.csv(
  file.path(
    deseq_dir,
    "DESeq2_normalized_counts.csv"
  ),
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()


## LOAD COMPLETE SPECIES VS REST DE RESULTS

species_vs_rest_results <- lapply(
  species_levels,
  function(sp) {
    
    readr::read_csv(
      file.path(
        species_vs_rest_dir,
        paste0(
          "DESeq2_",
          sp,
          "_vs_rest_all_results.csv"
        )
      ),
      show_col_types = FALSE
    ) %>%
      dplyr::mutate(
        species = sp,
        gene_match = tolower(
          trimws(
            gene_name
          )
        )
      )
  }
)

names(
  species_vs_rest_results
) <- species_levels


species_vs_rest_all <- dplyr::bind_rows(
  species_vs_rest_results
)


## CHECK

cat(
  "Genes in normalized-count matrix:",
  nrow(normalized_counts),
  "\n"
)

cat(
  "Species-vs-rest result rows:",
  nrow(species_vs_rest_all),
  "\n"
)

cat(
  "Visual genes in reference:",
  nrow(visual_gene_reference),
  "\n"
)



## VISUAL-GENE REPRESENTATION AUDIT

dataset_genes <- rownames(
  normalized_counts
)

dataset_gene_match <- tolower(
  trimws(
    dataset_genes
  )
)


visual_gene_audit <- visual_gene_reference %>%
  dplyr::mutate(
    
    present_in_DESeq2 =
      gene_match %in%
      dataset_gene_match,
    
    dataset_gene_name =
      dataset_genes[
        match(
          gene_match,
          dataset_gene_match
        )
      ]
  )


print(
  visual_gene_audit,
  n = Inf,
  width = Inf
)


readr::write_csv(
  visual_gene_audit,
  file.path(
    output_dir,
    "visual_gene_representation_audit.csv"
  )
)



## EXTRACT VISUAL-GENE DE RESULTS

visual_DE_results <- species_vs_rest_all %>%
  
  dplyr::filter(
    gene_match %in%
      visual_gene_reference$gene_match
  ) %>%
  
  dplyr::left_join(
    visual_gene_reference,
    by = "gene_match",
    suffix = c(
      "_DESeq2",
      "_reference"
    )
  ) %>%
  
  dplyr::mutate(
    
    significant =
      !is.na(padj) &
      padj < 0.05,
    
    direction =
      dplyr::case_when(
        
        significant &
          log2FoldChange > 0 ~
          paste(
            species,
            "higher"
          ),
        
        significant &
          log2FoldChange < 0 ~
          "Other species higher",
        
        TRUE ~
          "Not significant"
      )
  )


readr::write_csv(
  visual_DE_results,
  file.path(
    output_dir,
    "visual_gene_species_vs_rest_results.csv"
  )
)


## significant visual-system genes only

significant_visual_DE <- visual_DE_results %>%
  dplyr::filter(
    significant
  ) %>%
  dplyr::arrange(
    species,
    padj,
    dplyr::desc(
      abs(
        log2FoldChange
      )
    )
  )


print(
  significant_visual_DE,
  n = Inf,
  width = Inf
)


readr::write_csv(
  significant_visual_DE,
  file.path(
    output_dir,
    "significant_visual_gene_species_vs_rest_results.csv"
  )
)




## ANNOTATION/PARALOGUE AUDIT FOR VISUAL-SYSTEM GENES
## all gene names represented in the DESeq2 dataset

all_deseq_genes <- tibble::tibble(
  gene_name = rownames(
    normalized_counts
  )
) %>%
  dplyr::mutate(
    gene_match = tolower(
      trimws(
        gene_name
      )
    )
  )


## broad search patterns for visual-system gene families

visual_family_patterns <- tibble::tribble(
  ~visual_family,              ~pattern,
  
  "Opsins / photopigments",    "^(opn|rho|rh1|rh2|sws|lws)",
  "Transducin / G proteins",   "^(gnat|gngt|gnb)",
  "PDE6 family",               "^pde6",
  "CNG channels",              "^(cnga|cngb)",
  "Arrestin / recoverin",      "^(arr|sag|rcvrn)",
  "RGS proteins",              "^rgs",
  "Guanylate cyclase activators", "^guca",
  "Retinoid cycle",            "^(rpe|lrat|abca4|rdh|rbp|stra6|bco)",
  "Photoreceptor regulators",  "^(crx|nrl|nr2e3|otx|pax|six|rx)",
  "Other retinal candidates",  "^(vsx|foxn4|neurod|prdm)"
)


## search every DESeq2 gene name against each family

visual_family_hits <- purrr::map_dfr(
  seq_len(
    nrow(
      visual_family_patterns
    )
  ),
  function(i) {
    
    family_name <- visual_family_patterns$visual_family[i]
    family_pattern <- visual_family_patterns$pattern[i]
    
    all_deseq_genes %>%
      dplyr::filter(
        grepl(
          family_pattern,
          gene_match,
          ignore.case = TRUE
        )
      ) %>%
      dplyr::mutate(
        visual_family = family_name
      )
  }
) %>%
  dplyr::distinct(
    gene_name,
    visual_family,
    .keep_all = TRUE
  ) %>%
  dplyr::arrange(
    visual_family,
    gene_name
  )



## FLAG WHETHER EACH HIT WAS ALREADY IN THE PREDEFINED LIST

visual_family_hits <- visual_family_hits %>%
  dplyr::mutate(
    in_predefined_visual_list =
      gene_match %in%
      visual_gene_reference$gene_match
  )


## display
print(
  visual_family_hits,
  n = Inf,
  width = Inf
)


## export

readr::write_csv(
  visual_family_hits,
  file.path(
    output_dir,
    "visual_gene_family_annotation_audit.csv"
  )
)


## SUMMARY BY FAMILY

visual_family_summary <- visual_family_hits %>%
  dplyr::group_by(
    visual_family
  ) %>%
  dplyr::summarise(
    n_genes_found =
      dplyr::n(),
    
    n_already_in_predefined_list =
      sum(
        in_predefined_visual_list
      ),
    
    n_additional_candidates =
      sum(
        !in_predefined_visual_list
      ),
    
    additional_gene_names =
      paste(
        gene_name[
          !in_predefined_visual_list
        ],
        collapse = "; "
      ),
    
    .groups = "drop"
  )


print(
  visual_family_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  visual_family_summary,
  file.path(
    output_dir,
    "visual_gene_family_annotation_summary.csv"
  )
)




## CURATED VISUAL-SYSTEM GENE PANEL

## Candidate genes retained after inspection of the broad
## annotation/paralogue audit

curated_visual_genes <- tibble::tribble(
  ~gene_name, ~visual_category,
  
  ## classical visual opsins
  "opn1lw1",  "Opsin / photopigment",
  "opn1sw2",  "Opsin / photopigment",
  
  ## phototransduction G proteins
  "gnat1",    "Phototransduction",
  "gnat2",    "Phototransduction",
  "gngt1",    "Phototransduction",
  
  ## phosphodiesterase components
  "pde6a",    "Phototransduction",
  "pde6b",    "Phototransduction",
  "pde6c",    "Phototransduction",
  "pde6ga",   "Phototransduction",
  
  ## cyclic nucleotide-gated channels
  "cnga1a",   "CNG channel",
  "cnga1b",   "CNG channel",
  "cnga3a",   "CNG channel",
  "cngb1b",   "CNG channel",
  
  ## photoreceptor response/recovery
  "guca1a",   "Photoreceptor recovery",
  "guca1b",   "Photoreceptor recovery",
  "guca1c",   "Photoreceptor recovery",
  "guca1d",   "Photoreceptor recovery",
  "arr3a",    "Photoreceptor recovery",
  "rcvrn2",   "Photoreceptor recovery",
  "rcvrn3",   "Photoreceptor recovery",
  "rgs9a",    "Photoreceptor recovery",
  "rgs9b",    "Photoreceptor recovery",
  
  ## visual /retinoid cycle
  "rdh5",     "Visual cycle",
  "rdh8a",    "Visual cycle",
  "rdh12",    "Visual cycle",
  "rpe65a",   "Visual cycle",
  "rpe65b",   "Visual cycle",
  "abca4a",   "Visual cycle",
  "abca4b",   "Visual cycle",
  
  ## photoreceptor/retinal regulation
  "crx",      "Photoreceptor regulation",
  "nrl",      "Photoreceptor regulation",
  "nr2e3",    "Photoreceptor regulation",
  "otx2b",    "Photoreceptor regulation",
  "pax6b",    "Retinal regulation",
  "rx1",      "Retinal regulation",
  "rx2",      "Retinal regulation",
  "six6a",    "Retinal regulation",
  "six7",     "Retinal regulation"
)


## create matching key explicitly

curated_visual_genes <- curated_visual_genes %>%
  dplyr::mutate(
    gene_match = tolower(
      trimws(
        gene_name
      )
    )
  )


## check

print(
  names(
    curated_visual_genes
  )
)

print(
  curated_visual_genes,
  n = Inf
)


## CHECK CURATED PANEL REPRESENTATION

curated_visual_gene_audit <- curated_visual_genes %>%
  dplyr::mutate(
    
    present_in_DESeq2 =
      gene_match %in%
      dataset_gene_match,
    
    dataset_gene_name =
      dataset_genes[
        match(
          gene_match,
          dataset_gene_match
        )
      ]
  )


print(
  curated_visual_gene_audit,
  n = Inf,
  width = Inf
)


readr::write_csv(
  curated_visual_gene_audit,
  file.path(
    output_dir,
    "curated_visual_gene_audit.csv"
  )
)



## CURATED VISUAL GENES × SPECIES VS REST DESEQ2

curated_visual_DE <- species_vs_rest_all %>%
  
  dplyr::filter(
    gene_match %in%
      curated_visual_genes$gene_match
  ) %>%
  
  dplyr::left_join(
    curated_visual_genes,
    by = "gene_match"
  ) %>%
  
  dplyr::mutate(
    
    ## significance follows final DESeq2 workflow
    ## FDR < 0.05 only
    significant =
      !is.na(padj) &
      padj < 0.05,
    
    direction =
      dplyr::case_when(
        
        significant &
          log2FoldChange > 0 ~
          paste(
            species,
            "higher"
          ),
        
        significant &
          log2FoldChange < 0 ~
          paste(
            species,
            "lower"
          ),
        
        TRUE ~
          "Not significant"
      )
  )


## check object exists

print(
  dim(
    curated_visual_DE
  )
)

print(
  names(
    curated_visual_DE
  )
)


## export full table

readr::write_csv(
  curated_visual_DE,
  file.path(
    output_dir,
    "curated_visual_genes_all_species_vs_rest.csv"
  )
)



## SIGNIFICANT CURATED VISUAL-SYSTEM RESULTS
curated_visual_DE_significant <- curated_visual_DE %>%
  
  dplyr::filter(
    significant
  ) %>%
  
  dplyr::arrange(
    species,
    padj,
    dplyr::desc(
      abs(
        log2FoldChange
      )
    )
  )


print(
  curated_visual_DE_significant,
  n = Inf,
  width = Inf
)


readr::write_csv(
  curated_visual_DE_significant,
  file.path(
    output_dir,
    "curated_visual_genes_significant_species_vs_rest.csv"
  )
)


## GENE LEVEL EVIDENCE SUMMARY

curated_visual_gene_summary <- curated_visual_DE %>%
  
  dplyr::group_by(
    gene_match,
    gene_name.y,
    visual_category
  ) %>%
  
  dplyr::summarise(
    
    n_species_significant =
      sum(
        significant,
        na.rm = TRUE
      ),
    
    strongest_abs_log2FC =
      max(
        abs(
          log2FoldChange
        ),
        na.rm = TRUE
      ),
    
    best_FDR = if (
      all(is.na(padj))
    ) {
      NA_real_
    } else {
      min(
        padj,
        na.rm = TRUE
      )
    },
    
    significant_patterns =
      paste(
        unique(
          direction[
            significant
          ]
        ),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      n_species_significant
    ),
    best_FDR
  )


print(
  curated_visual_gene_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  curated_visual_gene_summary,
  file.path(
    output_dir,
    "curated_visual_gene_DE_evidence_summary.csv"
  )
)

