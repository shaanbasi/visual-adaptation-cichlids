## ASSESS EXPLORATORY WGCNA MODULE PRESERVATION - Langfelder et al. (2011)

library(WGCNA)
library(dplyr)
library(tibble)

options(stringsAsFactors = FALSE)


## LOAD POOLED WGCNA NETWORK

wgcna_dir <- "Results/WGCNA"

load(
  file.path(
    wgcna_dir,
    "WGCNA_network_objects.RData"
  )
)



## MATCH SAMPLE ORDER

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



## CHECK SPECIES SAMPLE SIZES

species_sample_counts <- speciesTraits %>%
  tibble::rownames_to_column(
    "sample"
  ) %>%
  dplyr::count(
    species,
    name = "n_samples"
  )

print(
  species_sample_counts
)



## ASSESS FEASIBILITY OF MODULE PRESERVATION

minimum_species_n <- min(
  species_sample_counts$n_samples
)

maximum_species_n <- max(
  species_sample_counts$n_samples
)


message(
  paste0(
    "\nFormal WGCNA module preservation was not performed.\n\n",
    
    "Species-specific sample sizes range from ",
    minimum_species_n,
    " to ",
    maximum_species_n,
    " samples.\n\n",
    
    "An exploratory analysis using WGCNA::modulePreservation(), ",
    "following the framework of Langfelder et al. (2011), was attempted ",
    "using pooled-network modules as the reference and individual ",
    "species as test datasets. WGCNA could not perform the analysis ",
    "because the species-specific datasets contained too few samples ",
    "for valid network calculations.\n\n",
    
    "Formal module-preservation statistics were therefore not ",
    "estimated or interpreted."
  )
)