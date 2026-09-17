library(WGCNA)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

## allow multithreading where supported
enableWGCNAThreads()

## INPUT AND OUTPUT SETTINGS

input_dir <- "Data/Processed"

output_dir <- "Results/WGCNA"

figure_dir <- "Figures/WGCNA"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

if (!dir.exists(figure_dir)) {
  dir.create(
    figure_dir,
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

wgcna_genes <- readRDS(
  file.path(
    input_dir,
    "wgcna_genes.rds"
  )
)

stopifnot(
  identical(
    colnames(datExpr),
    wgcna_genes
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


## CREATE SPECIES TRAIT MATRIX

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

stopifnot(
  identical(
    rownames(datExpr),
    rownames(speciesTraitMatrix)
  )
)



## SAMPLE CLUSTERING

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

pdf(
  file.path(
    figure_dir,
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
    figure_dir,
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
## selected after examining scale-free topology fit and mean connectivity diagnostics
softPower <- 20


## unsigned tutorial, but other authors recommend signed WGCNA author Q&A:
## "We generally recommend signed (or signed hybrid) networks because they produce
## modules that are easier to interpret biologically."


## Signed network used so that modules primarily contain genes with positively 
## correlated expression profiles.
## Genes with strong negative correlations are not treated as strongly connected 
## within the same module
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

##module dendogram and colour assignments
pdf(
  file.path(
    figure_dir,
    "module_dendrogram_colours.pdf"
  ),
  width = 14,
  height = 8
)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[
    net$blockGenes[[1]]
  ],
  groupLabels = "Module",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene clustering and WGCNA module assignments"
)

dev.off()


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


## Save
readr::write_csv(
  tibble::rownames_to_column(
    as.data.frame(MEs),
    var = "sample"
  ),
  file.path(
    output_dir,
    "module_eigengenes.csv"
  )
)

readr::write_csv(
  tibble::rownames_to_column(
    speciesTraitMatrix,
    var = "sample"
  ),
  file.path(
    output_dir,
    "species_trait_matrix.csv"
  )
)


##MODULE-SPECIES ASSOCIATIONS
moduleTraitStats <- bicorAndPvalue(
  x = MEs,
  y = speciesTraitMatrix,
  use = "pairwise.complete.obs",
  maxPOutliers = 0.10,
  robustX = TRUE,
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
    ),
    
    expression_pattern = case_when(
      correlation > 0 ~
        paste(
          species,
          "higher"
        ),
      
      correlation < 0 ~
        paste(
          species,
          "lower"
        ),
      
      TRUE ~
        "No direction"
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

## prepare nongrey module-species results
heatmap_df <- moduleSpeciesResults %>%
  
  dplyr::filter(
    module != "grey"
  ) %>%
  
  dplyr::mutate(
    
    species = factor(
      species,
      levels = c(
        "Ab",
        "Mz",
        "Nb",
        "On",
        "Pn"
      )
    ),
    
    module = factor(
      module,
      levels = rev(
        unique(
          sub(
            "^ME",
            "",
            rownames(moduleTraitCor)[
              rownames(moduleTraitCor) != "MEgrey"
            ]
          )
        )
      )
    ),
    
    label = sprintf(
      "%.2f",
      correlation
    ),
    
    significant =
      !is.na(FDR) &
      FDR < 0.05,
    
    label_sig = ifelse(
      significant,
      paste0(label, "*"),
      label
    )
  )


p_module_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = species,
    y = module,
    fill = correlation
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label = label_sig,
      fontface = ifelse(
        significant,
        "bold",
        "plain"
      )
    ),
    size = 3
  ) +
  
  scale_fill_gradient2(
    name = "Correlation",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -1,
      1
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 8
    ),
    
    axis.ticks = element_blank(),
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

## display
p_module_heatmap


## save
ggsave(
  filename = file.path(
    figure_dir,
    "Figure4A_module_species_heatmap.png"
  ),
  plot = p_module_heatmap,
  width = 8,
  height = 7,
  dpi = 600,
  bg = "white"
)



## SAVE CORE NETWORK OBJECTS

save(
  datExpr,
  speciesTraits,
  speciesTraitMatrix,
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

