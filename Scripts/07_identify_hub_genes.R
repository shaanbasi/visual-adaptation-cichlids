library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)

options(stringsAsFactors = FALSE)


## INPUT AND OUTPUT SETTINGS

wgcna_dir <- "Results/WGCNA"

output_dir <- file.path(
  wgcna_dir,
  "hub_genes"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}



## LOAD THE COMPLETED NETWORK

load(
  file.path(
    wgcna_dir,
    "WGCNA_network_objects.RData"
  )
)

required_objects <- c(
  "datExpr",
  "speciesTraits",
  "MEs",
  "geneModuleAssignments"
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


## check the metadata
colnames(speciesTraits)
table(speciesTraits$species)

## make sure sample order matches
speciesTraits <- speciesTraits[
  rownames(datExpr),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(datExpr),
    rownames(speciesTraits)
  )
)

## ensure species is a factor
speciesTraits$species <- factor(
  speciesTraits$species,
  levels = c(
    "Ab",
    "Mz",
    "Nb",
    "On",
    "Pn"
  )
)

## create one binary indicator column per species
speciesTraitMatrix <- model.matrix(
  ~ 0 + species,
  data = speciesTraits
)

colnames(speciesTraitMatrix) <- sub(
  "^species",
  "",
  colnames(speciesTraitMatrix)
)

speciesTraitMatrix <- as.data.frame(
  speciesTraitMatrix
)

rownames(speciesTraitMatrix) <- rownames(
  speciesTraits
)

## check the resulting indicators
head(speciesTraitMatrix)

colSums(speciesTraitMatrix)

stopifnot(
  identical(
    rownames(MEs),
    rownames(speciesTraitMatrix)
  )
)

stopifnot(
  all(
    rowSums(speciesTraitMatrix) == 1
  )
)


## RECALCULATE MODULE-SPECIES ASSOCIATIONS

moduleTraitCor <- WGCNA::bicor(
  MEs,
  speciesTraitMatrix,
  use = "pairwise.complete.obs",
  maxPOutliers = 0.10,
  robustY = FALSE
)

moduleTraitPvalue <- WGCNA::corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(datExpr)
)

moduleTraitFDR <- matrix(
  p.adjust(
    as.vector(moduleTraitPvalue),
    method = "BH"
  ),
  nrow = nrow(moduleTraitPvalue),
  ncol = ncol(moduleTraitPvalue),
  dimnames = dimnames(moduleTraitPvalue)
)

moduleSizes <- geneModuleAssignments %>%
  dplyr::transmute(
    Module = sub(
      "^ME",
      "",
      as.character(Module)
    )
  ) %>%
  dplyr::count(
    Module,
    name = "module_size"
  )

moduleSpeciesResults <- expand.grid(
  module = rownames(moduleTraitCor),
  species = colnames(moduleTraitCor),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    correlation = as.vector(
      moduleTraitCor
    ),
    absolute_correlation = abs(
      correlation
    ),
    pvalue = as.vector(
      moduleTraitPvalue
    ),
    FDR = as.vector(
      moduleTraitFDR
    ),
    direction = dplyr::case_when(
      correlation > 0 ~ "Positive",
      correlation < 0 ~ "Negative",
      TRUE ~ "None"
    ),
    module = sub(
      "^ME",
      "",
      module
    )
  ) %>%
  dplyr::left_join(
    moduleSizes,
    by = c(
      "module" = "Module"
    )
  ) %>%
  dplyr::arrange(
    species,
    FDR,
    dplyr::desc(
      absolute_correlation
    )
  )

stopifnot(
  setequal(
    unique(moduleSpeciesResults$species),
    c(
      "Ab",
      "Mz",
      "Nb",
      "On",
      "Pn"
    )
  )
)


## ANALYSIS THRESHOLDS

## module-species pair selection
module_FDR_threshold <- 0.10
module_correlation_threshold <- 0.50

## hub support thresholds
MM_threshold <- 0.80
aligned_GS_threshold <- 0.70
GS_pvalue_threshold <- 0.05

## export the highest-MM genes regardless of whether they pass all thresholds
top_n_genes <- 10



## STANDARDISE GENE-MODULE TABLE

geneModuleTable <- geneModuleAssignments %>%
  transmute(
    gene_id = as.character(Gene),
    module = sub(
      "^ME",
      "",
      as.character(Module)
    )
  )

stopifnot(
  setequal(
    colnames(datExpr),
    geneModuleTable$gene_id
  ),
  !anyDuplicated(
    geneModuleTable$gene_id
  ),
  identical(
    rownames(datExpr),
    rownames(MEs)
  ),
  identical(
    rownames(datExpr),
    rownames(speciesTraits)
  )
)



## SPECIES BY SPECIES MODULE SUMMARY

moduleSpeciesResults_nonGrey <- moduleSpeciesResults %>%
  dplyr::filter(
    module != "grey"
  ) %>%
  dplyr::mutate(
    correlation = as.numeric(correlation),
    absolute_correlation = as.numeric(absolute_correlation),
    pvalue = as.numeric(pvalue),
    FDR = as.numeric(FDR),
    
    passes_FDR =
      !is.na(FDR) &
      FDR <= module_FDR_threshold,
    
    passes_correlation =
      !is.na(absolute_correlation) &
      absolute_correlation >=
      module_correlation_threshold,
    
    selected_for_hub_analysis =
      passes_FDR &
      passes_correlation
  ) %>%
  dplyr::arrange(
    species,
    FDR,
    dplyr::desc(absolute_correlation)
  )


speciesModuleSummary <- moduleSpeciesResults_nonGrey %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    total_modules = dplyr::n(),
    
    modules_with_valid_results =
      sum(
        is.finite(absolute_correlation) &
          !is.na(FDR)
      ),
    
    modules_FDR_005 =
      sum(
        !is.na(FDR) &
          FDR < 0.05
      ),
    
    modules_FDR_010 =
      sum(
        !is.na(FDR) &
          FDR < 0.10
      ),
    
    modules_selected =
      sum(
        selected_for_hub_analysis,
        na.rm = TRUE
      ),
    
    strongest_module = {
      valid <- which(
        is.finite(absolute_correlation)
      )
      
      if (length(valid) == 0) {
        NA_character_
      } else {
        module[
          valid[
            which.max(
              absolute_correlation[valid]
            )
          ]
        ]
      }
    },
    
    strongest_correlation = {
      valid <- which(
        is.finite(absolute_correlation)
      )
      
      if (length(valid) == 0) {
        NA_real_
      } else {
        correlation[
          valid[
            which.max(
              absolute_correlation[valid]
            )
          ]
        ]
      }
    },
    
    strongest_FDR = {
      valid <- which(
        is.finite(absolute_correlation)
      )
      
      if (length(valid) == 0) {
        NA_real_
      } else {
        FDR[
          valid[
            which.max(
              absolute_correlation[valid]
            )
          ]
        ]
      }
    },
    
    .groups = "drop"
  )

print(
  speciesModuleSummary,
  width = Inf
)



## SELECT MODULE-SPECIES PAIRS

selectedModuleSpeciesPairs <-
  moduleSpeciesResults_nonGrey %>%
  filter(
    selected_for_hub_analysis
  ) %>%
  select(
    species,
    module,
    correlation,
    absolute_correlation,
    direction,
    pvalue,
    FDR,
    module_size
  )

print(
  tibble::as_tibble(selectedModuleSpeciesPairs),
  n = Inf,
  width = Inf
)

if (nrow(selectedModuleSpeciesPairs) == 0) {
  stop(
    paste0(
      "No module-species pairs passed the current thresholds. ",
      "Inspect moduleSpeciesResults_nonGrey before changing the criteria."
    )
  )
}


##CHECK - NO MODULES PASS THRESHOLDS
moduleSpeciesResults_nonGrey %>%
  dplyr::select(
    species,
    module,
    correlation,
    absolute_correlation,
    pvalue,
    FDR,
    passes_FDR,
    passes_correlation,
    selected_for_hub_analysis
  ) %>%
  dplyr::arrange(
    FDR,
    dplyr::desc(absolute_correlation)
  ) %>%
  tibble::as_tibble() %>%
  print(
    n = 50,
    width = Inf
  )

moduleSpeciesResults_nonGrey %>%
  dplyr::summarise(
    minimum_FDR = min(
      FDR,
      na.rm = TRUE
    ),
    maximum_absolute_correlation = max(
      absolute_correlation,
      na.rm = TRUE
    ),
    n_passing_FDR = sum(
      FDR <= module_FDR_threshold,
      na.rm = TRUE
    ),
    n_passing_correlation = sum(
      absolute_correlation >=
        module_correlation_threshold,
      na.rm = TRUE
    ),
    n_passing_both = sum(
      selected_for_hub_analysis,
      na.rm = TRUE
    )
  )


## CALCULATE MODULE MEMBERSHIP USING BICOR

geneMMStats <- bicorAndPvalue(
  x = datExpr,
  y = MEs,
  use = "pairwise.complete.obs",
  maxPOutliers = 0.10,
  robustX = TRUE,
  robustY = TRUE
)

geneMM <- as.data.frame(
  geneMMStats$bicor
)

geneMMpvalue <- as.data.frame(
  geneMMStats$p
)

rownames(geneMM) <- colnames(datExpr)
rownames(geneMMpvalue) <- colnames(datExpr)

colnames(geneMM) <- paste0(
  "MM.",
  sub(
    "^ME",
    "",
    colnames(MEs)
  )
)

colnames(geneMMpvalue) <- paste0(
  "p.MM.",
  sub(
    "^ME",
    "",
    colnames(MEs)
  )
)



## CALCULATE GENE SIGNIFICANCE FOR EACH SPECIES

geneTraitStats <- bicorAndPvalue(
  x = datExpr,
  y = speciesTraitMatrix,
  use = "pairwise.complete.obs",
  maxPOutliers = 0.10,
  robustX = TRUE,
  
  ## Species traits are binary
  robustY = FALSE
)

geneGS <- as.data.frame(
  geneTraitStats$bicor
)

geneGSpvalue <- as.data.frame(
  geneTraitStats$p
)

rownames(geneGS) <- colnames(datExpr)
rownames(geneGSpvalue) <- colnames(datExpr)

colnames(geneGS) <- paste0(
  "GS.",
  colnames(speciesTraitMatrix)
)

colnames(geneGSpvalue) <- paste0(
  "p.GS.",
  colnames(speciesTraitMatrix)
)



## BUILD THE MASTER GENE METRICS TABLE

hubMetrics <- geneModuleTable %>%
  left_join(
    geneMM %>%
      rownames_to_column(
        "gene_id"
      ),
    by = "gene_id"
  ) %>%
  left_join(
    geneMMpvalue %>%
      rownames_to_column(
        "gene_id"
      ),
    by = "gene_id"
  ) %>%
  left_join(
    geneGS %>%
      rownames_to_column(
        "gene_id"
      ),
    by = "gene_id"
  ) %>%
  left_join(
    geneGSpvalue %>%
      rownames_to_column(
        "gene_id"
      ),
    by = "gene_id"
  )

stopifnot(
  nrow(hubMetrics) == ncol(datExpr),
  !anyDuplicated(
    hubMetrics$gene_id
  )
)



## EXTRACT METRICS FOR EACH MODULE-SPECIES PAIR

focalHubMetrics <- map_dfr(
  seq_len(
    nrow(selectedModuleSpeciesPairs)
  ),
  function(i) {
    
    focal_module <-
      selectedModuleSpeciesPairs$module[i]
    
    focal_species <-
      selectedModuleSpeciesPairs$species[i]
    
    module_correlation <-
      selectedModuleSpeciesPairs$correlation[i]
    
    MM_column <- paste0(
      "MM.",
      focal_module
    )
    
    MM_pvalue_column <- paste0(
      "p.MM.",
      focal_module
    )
    
    GS_column <- paste0(
      "GS.",
      focal_species
    )
    
    GS_pvalue_column <- paste0(
      "p.GS.",
      focal_species
    )
    
    required_columns <- c(
      MM_column,
      MM_pvalue_column,
      GS_column,
      GS_pvalue_column
    )
    
    missing_columns <- setdiff(
      required_columns,
      colnames(hubMetrics)
    )
    
    if (length(missing_columns) > 0) {
      stop(
        "Missing columns for ",
        focal_module,
        " and ",
        focal_species,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    }
    
    hubMetrics %>%
      filter(
        module == focal_module
      ) %>%
      transmute(
        gene_id,
        module,
        species = focal_species,
        
        module_size =
          selectedModuleSpeciesPairs$module_size[i],
        
        module_trait_correlation =
          module_correlation,
        
        module_trait_direction =
          selectedModuleSpeciesPairs$direction[i],
        
        module_trait_pvalue =
          selectedModuleSpeciesPairs$pvalue[i],
        
        module_trait_FDR =
          selectedModuleSpeciesPairs$FDR[i],
        
        MM =
          .data[[MM_column]],
        
        MM_pvalue =
          .data[[MM_pvalue_column]],
        
        GS =
          .data[[GS_column]],
        
        GS_pvalue =
          .data[[GS_pvalue_column]],
        
        ## positive value = the gene follows the direction
        ## of the module-species relationship
        aligned_GS =
          GS *
          sign(
            module_correlation
          ),
        
        absolute_GS =
          abs(GS)
      )
  }
)



## RANK HUB GENES

rankedHubs <- focalHubMetrics %>%
  mutate(
    passes_MM =
      MM >= MM_threshold,
    
    passes_aligned_GS =
      aligned_GS >=
      aligned_GS_threshold,
    
    passes_GS_pvalue =
      GS_pvalue <=
      GS_pvalue_threshold,
    
    passes_all_support_criteria =
      passes_MM &
      passes_aligned_GS &
      passes_GS_pvalue,
    
    ## exploratory only - MM remains the primary ranking measure
    exploratory_joint_score =
      MM * aligned_GS
  ) %>%
  group_by(
    species,
    module
  ) %>%
  arrange(
    desc(MM),
    desc(aligned_GS),
    GS_pvalue,
    gene_id,
    .by_group = TRUE
  ) %>%
  mutate(
    MM_rank = row_number()
  ) %>%
  ungroup()



##CREATE HUB-GENE OUTPUT SETS

##highest MM genes for each module-species pair
topRankedHubs <- rankedHubs %>%
  group_by(
    species,
    module
  ) %>%
  slice_head(
    n = top_n_genes
  ) %>%
  ungroup()

## genes passing all supporting thresholds
supportedHubs <- rankedHubs %>%
  filter(
    passes_all_support_criteria
  ) %>%
  group_by(
    species,
    module
  ) %>%
  arrange(
    desc(MM),
    desc(aligned_GS),
    GS_pvalue,
    .by_group = TRUE
  ) %>%
  mutate(
    supported_hub_rank =
      row_number()
  ) %>%
  ungroup()



## COUNT SUPPORTED HUBS FOR EACH MODULE

moduleHubCounts <- selectedModuleSpeciesPairs %>%
  select(
    species,
    module,
    module_size,
    correlation,
    FDR
  ) %>%
  left_join(
    rankedHubs %>%
      group_by(
        species,
        module
      ) %>%
      summarise(
        genes_passing_MM =
          sum(
            passes_MM,
            na.rm = TRUE
          ),
        
        genes_passing_aligned_GS =
          sum(
            passes_aligned_GS,
            na.rm = TRUE
          ),
        
        genes_passing_GS_pvalue =
          sum(
            passes_GS_pvalue,
            na.rm = TRUE
          ),
        
        supported_hub_count =
          sum(
            passes_all_support_criteria,
            na.rm = TRUE
          ),
        
        maximum_MM =
          max(
            MM,
            na.rm = TRUE
          ),
        
        median_MM =
          median(
            MM,
            na.rm = TRUE
          ),
        
        maximum_aligned_GS =
          max(
            aligned_GS,
            na.rm = TRUE
          ),
        
        .groups = "drop"
      ),
    by = c(
      "species",
      "module"
    )
  ) %>%
  arrange(
    species,
    desc(supported_hub_count),
    FDR,
    desc(abs(correlation))
  )


print(
  tibble::as_tibble(moduleHubCounts),
  n = Inf,
  width = Inf
)



## SUMMARISE HUB COUNTS BY SPECIES

speciesHubSummary <- moduleHubCounts %>%
  group_by(species) %>%
  summarise(
    selected_modules = n(),
    
    modules_with_supported_hubs =
      sum(
        supported_hub_count > 0,
        na.rm = TRUE
      ),
    
    total_supported_module_gene_pairs =
      sum(
        supported_hub_count,
        na.rm = TRUE
      ),
    
    minimum_supported_hubs =
      min(
        supported_hub_count,
        na.rm = TRUE
      ),
    
    median_supported_hubs =
      median(
        supported_hub_count,
        na.rm = TRUE
      ),
    
    maximum_supported_hubs =
      max(
        supported_hub_count,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  speciesHubSummary,
  width = Inf
)



## THRESHOLD SENSITIVITY

hubThresholdSensitivity <- rankedHubs %>%
  group_by(
    species,
    module
  ) %>%
  summarise(
    module_size = n(),
    
    n_MM_070_GS_060 =
      sum(
        MM >= 0.70 &
          aligned_GS >= 0.60,
        na.rm = TRUE
      ),
    
    n_MM_075_GS_065 =
      sum(
        MM >= 0.75 &
          aligned_GS >= 0.65,
        na.rm = TRUE
      ),
    
    n_MM_080_GS_070 =
      sum(
        MM >= 0.80 &
          aligned_GS >= 0.70,
        na.rm = TRUE
      ),
    
    n_MM_080_GS_080 =
      sum(
        MM >= 0.80 &
          aligned_GS >= 0.80,
        na.rm = TRUE
      ),
    
    n_MM_090_GS_080 =
      sum(
        MM >= 0.90 &
          aligned_GS >= 0.80,
        na.rm = TRUE
      ),
    
    n_full_current_criteria =
      sum(
        MM >= MM_threshold &
          aligned_GS >= aligned_GS_threshold &
          GS_pvalue <= GS_pvalue_threshold,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    species,
    module
  )



## TOP MODULES FOR EVERY SPECIES

## this table includes exploratory modules that may not meet the selected FDR
## threshold - useful for seeing the strongest relationships for each
## species without claiming that all are significant

topModulesPerSpecies <-
  moduleSpeciesResults_nonGrey %>%
  group_by(species) %>%
  arrange(
    FDR,
    desc(absolute_correlation),
    .by_group = TRUE
  ) %>%
  mutate(
    module_rank_for_species =
      row_number()
  ) %>%
  filter(
    module_rank_for_species <= 5
  ) %>%
  ungroup()



##EXPORT TABLES

write_csv(
  moduleSpeciesResults_nonGrey,
  file.path(
    output_dir,
    "01_all_nonGrey_module_species_results.csv"
  )
)

write_csv(
  speciesModuleSummary,
  file.path(
    output_dir,
    "02_species_module_summary.csv"
  )
)

write_csv(
  topModulesPerSpecies,
  file.path(
    output_dir,
    "03_top_five_modules_per_species.csv"
  )
)

write_csv(
  selectedModuleSpeciesPairs,
  file.path(
    output_dir,
    "04_selected_module_species_pairs.csv"
  )
)

write_csv(
  hubMetrics,
  file.path(
    output_dir,
    "05_all_gene_MM_and_GS_metrics.csv"
  )
)

write_csv(
  rankedHubs,
  file.path(
    output_dir,
    "06_all_ranked_genes_selected_pairs.csv"
  )
)

write_csv(
  topRankedHubs,
  file.path(
    output_dir,
    "07_top_MM_ranked_genes_per_pair.csv"
  )
)

write_csv(
  supportedHubs,
  file.path(
    output_dir,
    "08_supported_hub_genes.csv"
  )
)

write_csv(
  moduleHubCounts,
  file.path(
    output_dir,
    "09_supported_hub_counts_by_module.csv"
  )
)

write_csv(
  speciesHubSummary,
  file.path(
    output_dir,
    "10_supported_hub_counts_by_species.csv"
  )
)

write_csv(
  hubThresholdSensitivity,
  file.path(
    output_dir,
    "11_hub_threshold_sensitivity.csv"
  )
)



## SAVE HUB-GENE OBJECTS
save(
  speciesTraitMatrix,
  moduleTraitCor,
  moduleTraitPvalue,
  moduleTraitFDR,
  moduleSpeciesResults,
  geneMM,
  geneMMpvalue,
  geneGS,
  geneGSpvalue,
  hubMetrics,
  selectedModuleSpeciesPairs,
  focalHubMetrics,
  rankedHubs,
  topRankedHubs,
  supportedHubs,
  moduleHubCounts,
  speciesHubSummary,
  hubThresholdSensitivity,
  
  file = file.path(
    output_dir,
    "WGCNA_hub_gene_objects.RData"
  )
)

saveRDS(
  supportedHubs,
  file.path(
    output_dir,
    "supported_hub_genes.rds"
  )
)



##how does this compare?
## inspect top-ranked hub genes
topRankedHubs %>%
  dplyr::select(
    species,
    module,
    gene_id,
    MM,
    aligned_GS,
    GS_pvalue,
    module_trait_FDR,
    MM_rank
  ) %>%
  dplyr::arrange(
    species,
    module,
    MM_rank
  ) %>%
  print(
    n = 30,
    width = Inf
  )

##overlap with visual-gene list
visual_hub_overlap <- supportedHubs %>%
  dplyr::mutate(
    gene_upper = toupper(trimws(gene_id))
  ) %>%
  dplyr::filter(
    gene_upper %in% visual_genes
  ) %>%
  dplyr::select(
    species,
    module,
    gene_id,
    MM,
    aligned_GS,
    GS_pvalue,
    module_trait_FDR
  ) %>%
  dplyr::arrange(
    species,
    module,
    dplyr::desc(MM)
  )

visual_hub_overlap

## save
write.csv(
  visual_hub_overlap,
  "Results/WGCNA/hub_genes/12_visual_candidate_hub_overlap.csv",
  row.names = FALSE
)


##FINAL SUMMARY
cat(
  "\n============================================\n"
)

cat(
  "Hub-gene analysis complete\n"
)

cat(
  "Selected module-species pairs:",
  nrow(selectedModuleSpeciesPairs),
  "\n"
)

cat(
  "Supported hub-gene rows:",
  nrow(supportedHubs),
  "\n"
)

cat(
  "MM threshold:",
  MM_threshold,
  "\n"
)

cat(
  "Aligned GS threshold:",
  aligned_GS_threshold,
  "\n"
)

cat(
  "GS P-value threshold:",
  GS_pvalue_threshold,
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


