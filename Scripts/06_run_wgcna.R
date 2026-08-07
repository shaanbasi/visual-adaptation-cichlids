library(WGCNA)
library(dplyr)
library(tibble)
library(readr)

options(stringsAsFactors = FALSE)

## allow multithreading where supported
enableWGCNAThreads()

## INPUT AND OUTPUT SETTINGS

input_dir <- "Data/Processed"

output_dir <- "Results/WGCNA"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


## LOAD PREPARED WGCNA INPUT OBJECTS

datExpr <- readRDS(
  file.path(
    input_dir,
    "datExpr.rds"
  )
)

speciesTraits <- readRDS(
  file.path(
    input_dir,
    "speciesTraits.rds"
  )
)


##BASIC VALIDATION

datExpr <- as.data.frame(datExpr)
speciesTraits <- as.data.frame(speciesTraits)

stopifnot(
  identical(
    rownames(datExpr),
    rownames(speciesTraits)
  )
)

cat(
  "Samples:",
  nrow(datExpr),
  "\n"
)

cat(
  "Genes:",
  ncol(datExpr),
  "\n"
)

cat(
  "Species traits:",
  ncol(speciesTraits),
  "\n"
)



## CHECK GENES AND SAMPLES

goodSamplesGenesResult <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

if (!goodSamplesGenesResult$allOK) {
  
  if (sum(!goodSamplesGenesResult$goodGenes) > 0) {
    message(
      "Removing genes: ",
      paste(
        colnames(datExpr)[
          !goodSamplesGenesResult$goodGenes
        ],
        collapse = ", "
      )
    )
  }
  
  if (sum(!goodSamplesGenesResult$goodSamples) > 0) {
    message(
      "Removing samples: ",
      paste(
        rownames(datExpr)[
          !goodSamplesGenesResult$goodSamples
        ],
        collapse = ", "
      )
    )
  }
  
  datExpr <- datExpr[
    goodSamplesGenesResult$goodSamples,
    goodSamplesGenesResult$goodGenes
  ]
  
  speciesTraits <- speciesTraits[
    goodSamplesGenesResult$goodSamples,
    ,
    drop = FALSE
  ]
}

stopifnot(
  identical(
    rownames(datExpr),
    rownames(speciesTraits)
  )
)



## SAMPLE CLUSTERING

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

pdf(
  file.path(
    output_dir,
    "sample_clustering.pdf"
  ),
  width = 10,
  height = 7
)

plot(
  sampleTree,
  main = "Sample clustering",
  xlab = "",
  sub = "",
  cex = 0.8
)

dev.off()



##SOFT THRESHOLD ANALYSIS

## network settings
network_type <- "signed"

powers <- c(
  1:10,
  seq(
    12,
    30,
    by = 2
  )
)

softThresholdResults <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = network_type,
  corFnc = "bicor",
  corOptions = list(
    use = "pairwise.complete.obs",
    maxPOutliers = 0.10
  ),
  verbose = 5
)


softThresholdTable <- softThresholdResults$fitIndices

write_csv(
  softThresholdTable,
  file.path(
    output_dir,
    "soft_threshold_results.csv"
  )
)

pdf(
  file.path(
    output_dir,
    "soft_threshold_diagnostics.pdf"
  ),
  width = 12,
  height = 6
)

par(
  mfrow = c(1, 2)
)

plot(
  softThresholdTable$Power,
  -sign(
    softThresholdTable$slope
  ) *
    softThresholdTable$SFT.R.sq,
  xlab = "Soft-threshold power",
  ylab = "Signed R²",
  type = "n",
  main = "Scale-free topology fit"
)

text(
  softThresholdTable$Power,
  -sign(
    softThresholdTable$slope
  ) *
    softThresholdTable$SFT.R.sq,
  labels = softThresholdTable$Power,
  cex = 0.8
)

abline(
  h = 0.80,
  lty = 2
)

plot(
  softThresholdTable$Power,
  softThresholdTable$mean.k.,
  xlab = "Soft-threshold power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  softThresholdTable$Power,
  softThresholdTable$mean.k.,
  labels = softThresholdTable$Power,
  cex = 0.8
)

dev.off()


## SOFT-THRESHOLD SETTINGS

softPower <- 20


## unsigned tutorial, but other authors recommend signed
## WGCNA author Q&A: "We generally recommend signed 
## (or signed hybrid) networks because they produce modules 
## that are easier to interpret biologically."
## Ultimately, I want to see which genes are working together,
## It doesn't particularly matter if they behave in opposite directions
network_type <- "signed"
TOM_type <- "signed"


## CONSTRUCT THE NETWORK

set.seed(123)

net <- blockwiseModules(
  datExpr,
  
  power = softPower,
  
  networkType = network_type,
  TOMType = TOM_type,
  
  corType = "bicor",
  maxPOutliers = 0.10,
  
  minModuleSize = 30,
  deepSplit = 2,
  
  mergeCutHeight = 0.25,
  
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  
  maxBlockSize = ncol(datExpr),
  
  saveTOMs = FALSE,
  
  verbose = 5
)



## EXTRACT MODULE ASSIGNMENTS

moduleColors <- labels2colors(
  net$colors
)

geneModuleAssignments <- tibble(
  Gene = colnames(datExpr),
  Module = moduleColors
)

moduleSizes <- geneModuleAssignments %>%
  dplyr::count(
    Module,
    name = "module_size"
  ) %>%
  dplyr::arrange(
    dplyr::desc(module_size)
  )

print(
  moduleSizes,
  n = Inf
)

write_csv(
  geneModuleAssignments,
  file.path(
    output_dir,
    "gene_module_assignments.csv"
  )
)

write_csv(
  moduleSizes,
  file.path(
    output_dir,
    "module_sizes.csv"
  )
)



## CALCULATE MODULE EIGENGENES

MElist <- moduleEigengenes(
  datExpr,
  colors = moduleColors,
  impute = TRUE,
  nPC = 1
)

MEs <- orderMEs(
  MElist$eigengenes
)

stopifnot(
  identical(
    rownames(datExpr),
    rownames(MEs)
  )
)



##MODULE-SPECIES ASSOCIATIONS

moduleTraitStats <- bicorAndPvalue(
  x = MEs,
  y = speciesTraits,
  use = "pairwise.complete.obs",
  maxPOutliers = 0.10,
  robustX = TRUE,
  
  ## Species columns are binary dummy variables
  robustY = FALSE
)

moduleTraitCor <- moduleTraitStats$bicor
moduleTraitPvalue <- moduleTraitStats$p



##BH FDR CORRECTION

moduleTraitFDR <- matrix(
  p.adjust(
    as.vector(moduleTraitPvalue),
    method = "BH"
  ),
  nrow = nrow(moduleTraitPvalue),
  ncol = ncol(moduleTraitPvalue),
  dimnames = dimnames(moduleTraitPvalue)
)



## TIDY MODULE-SPECIES RESULTS

moduleSpeciesResults <- expand.grid(
  eigengene = rownames(moduleTraitCor),
  species = colnames(moduleTraitCor),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    module = sub(
      "^ME",
      "",
      eigengene
    ),
    
    correlation =
      moduleTraitCor[
        eigengene,
        species
      ],
    
    pvalue =
      moduleTraitPvalue[
        eigengene,
        species
      ],
    
    FDR =
      moduleTraitFDR[
        eigengene,
        species
      ],
    
    absolute_correlation =
      abs(correlation),
    
    direction = case_when(
      correlation > 0 ~ "Positive",
      correlation < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) %>%
  ungroup() %>%
  left_join(
    moduleSizes,
    by = c(
      "module" = "Module"
    )
  ) %>%
  arrange(
    species,
    FDR,
    desc(absolute_correlation)
  )

write_csv(
  moduleSpeciesResults,
  file.path(
    output_dir,
    "module_species_results.csv"
  )
)


## MODULE-TRAIT HEATMAP

heatmapText <- paste0(
  signif(
    moduleTraitCor,
    2
  ),
  "\n(",
  signif(
    moduleTraitFDR,
    2
  ),
  ")"
)

dim(heatmapText) <- dim(moduleTraitCor)

pdf(
  file.path(
    output_dir,
    "module_species_heatmap.pdf"
  ),
  width = 9,
  height = 11
)

par(
  mar = c(
    6,
    9,
    3,
    3
  )
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  
  xLabels = colnames(speciesTraits),
  yLabels = rownames(moduleTraitCor),
  
  ySymbols = rownames(moduleTraitCor),
  
  colorLabels = FALSE,
  
  colors = blueWhiteRed(50),
  
  textMatrix = heatmapText,
  
  setStdMargins = FALSE,
  
  cex.text = 0.6,
  
  zlim = c(
    -1,
    1
  ),
  
  main = paste(
    "Module-species relationships",
    "\ncorrelation (BH FDR)"
  )
)

dev.off()



## SAVE CORE NETWORK OBJECTS

save(
  datExpr,
  speciesTraits,
  softPower,
  network_type,
  TOM_type,
  net,
  moduleColors,
  geneModuleAssignments,
  moduleSizes,
  MEs,
  moduleTraitCor,
  moduleTraitPvalue,
  moduleTraitFDR,
  moduleSpeciesResults,
  
  file = file.path(
    output_dir,
    "WGCNA_network_objects.RData"
  )
)

saveRDS(
  net,
  file.path(
    output_dir,
    "WGCNA_network.rds"
  )
)



## FINAL SUMMARY

cat(
  "\n============================================\n"
)

cat(
  "WGCNA network construction complete\n"
)

cat(
  "Samples:",
  nrow(datExpr),
  "\n"
)

cat(
  "Genes:",
  ncol(datExpr),
  "\n"
)

cat(
  "Non-grey modules:",
  length(
    setdiff(
      unique(moduleColors),
      "grey"
    )
  ),
  "\n"
)

cat(
  "Grey genes:",
  sum(moduleColors == "grey"),
  "\n"
)

cat(
  "Soft-threshold power:",
  softPower,
  "\n"
)

cat(
  "Network type:",
  network_type,
  "\n"
)

cat(
  "Saved to:",
  output_dir,
  "\n"
)

cat(
  "============================================\n\n"
)
