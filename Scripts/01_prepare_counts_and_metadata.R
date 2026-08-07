library(dplyr)
library(tibble)

## ANALYSIS SETTINGS
## species levels
species_levels <- c("Ab", "Mz", "Nb", "On", "Pn")

## minimum raw count required
minimum_count <- 5

## minimum number of samples meeting the count threshold
minimum_samples <- 2

## output folder
prep_output_dir <- "Data/Processed"

dir.create(
  prep_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


##read metadata
coldata <- read.csv(
  "Data/Raw/metadata.csv",
  row.names = 1,
  check.names = FALSE
)

coldata

## check for species column 
if (!"species" %in% colnames(coldata)) {
  stop("metadata.csv must contain a column named 'species'.")
}

## ensure species is a factor
coldata$species <- factor(
  coldata$species,
  levels = species_levels
)

## check for unexpected or missing species
if (anyNA(coldata$species)) {
  stop(
    paste0(
      "Unexpected or missing species labels were found. Expected: ",
      paste(species_levels, collapse = ", ")
    )
  )
}


##READ HTSEQ COUNT FILES
ab <- read.delim("Data/Counts/Ab_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
mz <- read.delim("Data/Counts/Mz_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
nb <- read.delim("Data/Counts/Nb_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
on <- read.delim("Data/Counts/On_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
pn <- read.delim("Data/Counts/Pn_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)

dim(ab)
dim(mz)
dim(nb)
dim(on)
dim(pn)

##check sample names
colnames(ab)
colnames(mz)
colnames(nb)
colnames(on)
colnames(pn)

## check gene IDs
head(rownames(ab))
head(rownames(mz))
head(rownames(nb))
head(rownames(on))
head(rownames(pn))

## remove rows such as __no_feature, __ambiguous and __not_aligned
ab <- ab[
  !grepl("^__", rownames(ab)),
  ,
  drop = FALSE
]

mz <- mz[
  !grepl("^__", rownames(mz)),
  ,
  drop = FALSE
]

nb <- nb[
  !grepl("^__", rownames(nb)),
  ,
  drop = FALSE
]

on <- on[
  !grepl("^__", rownames(on)),
  ,
  drop = FALSE
]

pn <- pn[
  !grepl("^__", rownames(pn)),
  ,
  drop = FALSE
]



## STEP 1 - ortholog mapping

##clean up bed files
ab_bed <- read.delim("Data/Bed/Haplochromis_burtoni.AstBur1.0.101.bed", header = FALSE)
head(ab_bed)
dim(ab_bed)

mz_bed <- read.delim("Data/Bed/Maylandia_zebra.M_zebra_UMD2a.101.bed", header = FALSE)
nb_bed <- read.delim("Data/Bed/Neolamprologus_brichardi.NeoBri1.0.101.bed", header = FALSE)
on_bed <- read.delim("Data/Bed/Oreochromis_niloticus.O_niloticus_UMD_NMBU.101.bed", header = FALSE)
pn_bed <- read.delim("Data/Bed/Pundamilia_nyererei.PunNye1.0.101.bed", header = FALSE)

## assign BED column names
bed_column_names <- c(
  "chr",
  "start",
  "end",
  "dot",
  "score",
  "strand",
  "gene_id",
  "gene_name"
)

colnames(ab_bed) <- bed_column_names
colnames(mz_bed) <- bed_column_names
colnames(nb_bed) <- bed_column_names
colnames(on_bed) <- bed_column_names
colnames(pn_bed) <- bed_column_names


## ##merge bed annot into each count table
## add gene IDs from row names
ab$gene_id <- rownames(ab)
mz$gene_id <- rownames(mz)
nb$gene_id <- rownames(nb)
on$gene_id <- rownames(on)
pn$gene_id <- rownames(pn)

## merge gene names into count tables
ab_annot <- merge(
  ab,
  ab_bed[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

mz_annot <- merge(
  mz,
  mz_bed[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

nb_annot <- merge(
  nb,
  nb_bed[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

on_annot <- merge(
  on,
  on_bed[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

pn_annot <- merge(
  pn,
  pn_bed[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

## check missing
sum(is.na(ab_annot$gene_name) | ab_annot$gene_name == "")
sum(is.na(mz_annot$gene_name) | mz_annot$gene_name == "")
sum(is.na(nb_annot$gene_name) | nb_annot$gene_name == "")
sum(is.na(on_annot$gene_name) | on_annot$gene_name == "")
sum(is.na(pn_annot$gene_name) | pn_annot$gene_name == "")

##remove genes without gene name
ab_named <- ab_annot[!is.na(ab_annot$gene_name) & ab_annot$gene_name != "", ]
mz_named <- mz_annot[!is.na(mz_annot$gene_name) & mz_annot$gene_name != "", ]
nb_named <- nb_annot[!is.na(nb_annot$gene_name) & nb_annot$gene_name != "", ]
on_named <- on_annot[!is.na(on_annot$gene_name) & on_annot$gene_name != "", ]
pn_named <- pn_annot[!is.na(pn_annot$gene_name) & pn_annot$gene_name != "", ]

## check remaining
nrow(ab_named)
nrow(mz_named)
nrow(nb_named)
nrow(on_named)
nrow(pn_named)


## COMBINE DUPLICATED GENE NAMES
## find genes shared across all species
collapse_gene_counts <- function(
    annotated_counts,
    sample_columns
) {
  
  annotated_counts %>%
    select(
      gene_name,
      all_of(sample_columns)
    ) %>%
    group_by(gene_name) %>%
    summarise(
      across(
        all_of(sample_columns),
        ~ sum(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    as.data.frame() %>%
    column_to_rownames("gene_name")
}


## identify sample columns
ab_sample_columns <- setdiff(
  colnames(ab),
  "gene_id"
)

mz_sample_columns <- setdiff(
  colnames(mz),
  "gene_id"
)

nb_sample_columns <- setdiff(
  colnames(nb),
  "gene_id"
)

on_sample_columns <- setdiff(
  colnames(on),
  "gene_id"
)

pn_sample_columns <- setdiff(
  colnames(pn),
  "gene_id"
)


## collapse duplicate gene names
ab_named_counts <- collapse_gene_counts(
  ab_named,
  ab_sample_columns
)

mz_named_counts <- collapse_gene_counts(
  mz_named,
  mz_sample_columns
)

nb_named_counts <- collapse_gene_counts(
  nb_named,
  nb_sample_columns
)

on_named_counts <- collapse_gene_counts(
  on_named,
  on_sample_columns
)

pn_named_counts <- collapse_gene_counts(
  pn_named,
  pn_sample_columns
)

## check that duplicate gene names are gone
sum(duplicated(rownames(ab_named_counts)))
sum(duplicated(rownames(mz_named_counts)))
sum(duplicated(rownames(nb_named_counts)))
sum(duplicated(rownames(on_named_counts)))
sum(duplicated(rownames(pn_named_counts)))
## all zero


## FIND GENES SHARED ACROSS ALL SPECIES
##gene names into row names
shared_genes <- Reduce(
  intersect,
  list(
    rownames(ab_named_counts),
    rownames(mz_named_counts),
    rownames(nb_named_counts),
    rownames(on_named_counts),
    rownames(pn_named_counts)
  )
)

length(shared_genes)
head(shared_genes, 20)


##RETAIN SHARED GENES IN THE SAME ORDER
##reorder rows
ab_shared <- ab_named_counts[
  shared_genes,
  ,
  drop = FALSE
]

mz_shared <- mz_named_counts[
  shared_genes,
  ,
  drop = FALSE
]

nb_shared <- nb_named_counts[
  shared_genes,
  ,
  drop = FALSE
]

on_shared <- on_named_counts[
  shared_genes,
  ,
  drop = FALSE
]

pn_shared <- pn_named_counts[
  shared_genes,
  ,
  drop = FALSE
]


## check identical row order
stopifnot(
  identical(rownames(ab_shared), shared_genes),
  identical(rownames(mz_shared), shared_genes),
  identical(rownames(nb_shared), shared_genes),
  identical(rownames(on_shared), shared_genes),
  identical(rownames(pn_shared), shared_genes)
)
## all true


##COMBINE INTO ONE COUNT MATRIX
##check with metadata
count_mat <- cbind(
  ab_shared,
  mz_shared,
  nb_shared,
  on_shared,
  pn_shared
)

count_mat <- as.matrix(count_mat)

## DESeq2 requires raw integer counts
storage.mode(count_mat) <- "integer"

dim(count_mat)
head(count_mat)


## CHECK COUNTS AGAINST METADATA
## check that count matrix and metadata contain the same samples
if (!setequal(colnames(count_mat), rownames(coldata))) {
  
  samples_missing_from_metadata <- setdiff(
    colnames(count_mat),
    rownames(coldata)
  )
  
  samples_missing_from_counts <- setdiff(
    rownames(coldata),
    colnames(count_mat)
  )
  
  stop(
    paste0(
      "Sample names do not match between count_mat and coldata.\n",
      "Missing from metadata: ",
      paste(samples_missing_from_metadata, collapse = ", "),
      "\nMissing from count matrix: ",
      paste(samples_missing_from_counts, collapse = ", ")
    )
  )
}

## reorder metadata to match count matrix columns
coldata <- coldata[
  colnames(count_mat),
  ,
  drop = FALSE
]

## confirm identical sample order
stopifnot(
  identical(
    colnames(count_mat),
    rownames(coldata)
  )
)


## LOW-COUNT FILTERING
## keep genes with count =/>5 in at least 2 samples.

## check dimensions before filtering
dim(count_mat)

## keep genes with at least 5 counts in at least 2 samples
keep <- rowSums(count_mat >= 5) >= 2

table(keep)

count_mat_filt <- count_mat[keep, ]

## dimensions after filtering
dim(count_mat_filt)

##check
nrow(count_mat) - nrow(count_mat_filt)
nrow(count_mat_filt)

##Genes were filtered prior to differential expression analysis to retain genes with at least 5 raw counts in at least 2 samples.



## RECORD FILTERING SUMMARY
filter_summary <- tibble(
  filtering_stage = c(
    "Shared named genes before low-count filtering",
    "Shared named genes after low-count filtering"
  ),
  
  n_genes = c(
    nrow(count_mat),
    nrow(count_mat_filt)
  )
)

filter_summary

write.csv(
  filter_summary,
  file.path(
    prep_output_dir,
    "Gene_filter_summary.csv"
  ),
  row.names = FALSE
)

##record the common number of genes used for each species
## all species will be analysed using the same filtered gene set

species_filter_summary <- tibble(
  species = species_levels,
  
  n_samples = as.integer(
    table(coldata$species)[species_levels]
  ),
  
  n_genes_entering_DESeq2 = nrow(count_mat_filt)
)

species_filter_summary

write.csv(
  species_filter_summary,
  file.path(
    prep_output_dir,
    "Genes_retained_per_species.csv"
  ),
  row.names = FALSE
)


## FINAL CHECKS

## check for missing values
if (anyNA(count_mat_filt)) {
  stop("count_mat_filt contains missing values.")
}

## check for negative counts
if (any(count_mat_filt < 0)) {
  stop("count_mat_filt contains negative values.")
}

## check that counts are integers
if (any(count_mat_filt != round(count_mat_filt))) {
  stop("count_mat_filt must contain raw integer counts.")
}

storage.mode(count_mat_filt) <- "integer"


##SAVE PREPARED OBJECTS
## save count matrix and metadata as CSV files
write.csv(
  count_mat_filt,
  file.path(
    prep_output_dir,
    "count_matrix_filtered.csv"
  ),
  row.names = TRUE
)

write.csv(
  coldata,
  file.path(
    prep_output_dir,
    "metadata_prepared.csv"
  ),
  row.names = TRUE
)


## save R objects for direct use in the DESeq2 script
saveRDS(
  count_mat_filt,
  file.path(
    prep_output_dir,
    "count_mat_filt.rds"
  )
)

saveRDS(
  coldata,
  file.path(
    prep_output_dir,
    "coldata.rds"
  )
)










