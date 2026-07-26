coldata <- read.csv("metadata.csv", row.names = 1)
coldata

ab <- read.delim("Ab_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
mz <- read.delim("Mz_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
nb <- read.delim("Nb_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
on <- read.delim("On_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)
pn <- read.delim("Pn_Eye_htseq_counts.out", row.names = 1, check.names = FALSE)

dim(ab)
dim(mz)
dim(nb)
dim(on)
dim(pn)

head(ab)
head(on)

colnames(ab)
colnames(mz)
colnames(nb)
colnames(on)
colnames(pn)

tail(rownames(on), 10)

head(rownames(ab)); head(rownames(mz)); head(rownames(nb)); head(rownames(on)); head(rownames(pn))
##looks correct


## STEP 1 - ortholog mapping

##clean up bed files
ab_bed <- read.delim("Haplochromis_burtoni.AstBur1.0.101.bed", header = FALSE)
head(ab_bed)
dim(ab_bed)

mz_bed <- read.delim("Maylandia_zebra.M_zebra_UMD2a.101.bed", header = FALSE)
nb_bed <- read.delim("Neolamprologus_brichardi.NeoBri1.0.101.bed", header = FALSE)
on_bed <- read.delim("Oreochromis_niloticus.O_niloticus_UMD_NMBU.101.bed", header = FALSE)
pn_bed <- read.delim("Pundamilia_nyererei.PunNye1.0.101.bed", header = FALSE)

colnames(ab_bed) <- c("chr", "start", "end", "dot", "score", "strand", "gene_id", "gene_name")
colnames(mz_bed) <- c("chr", "start", "end", "dot", "score", "strand", "gene_id", "gene_name")
colnames(nb_bed) <- c("chr", "start", "end", "dot", "score", "strand", "gene_id", "gene_name")
colnames(on_bed) <- c("chr", "start", "end", "dot", "score", "strand", "gene_id", "gene_name")
colnames(pn_bed) <- c("chr", "start", "end", "dot", "score", "strand", "gene_id", "gene_name")

## ##merge bed annot into each count table
ab$gene_id <- rownames(ab)
ab_annot <- merge(ab, ab_bed[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
head(ab_annot)

mz$gene_id <- rownames(mz)
mz_annot <- merge(mz, mz_bed[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
head(mz_annot)

nb$gene_id <- rownames(nb)
nb_annot <- merge(nb, nb_bed[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
head(nb_annot)

on$gene_id <- rownames(on)
on_annot <- merge(on, on_bed[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
head(on_annot)

pn$gene_id <- rownames(pn)
pn_annot <- merge(pn, pn_bed[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
head(pn_annot)

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

##find genes shared across all species
shared_genes <- Reduce(intersect, list(
  ab_named$gene_name,
  mz_named$gene_name,
  nb_named$gene_name,
  on_named$gene_name,
  pn_named$gene_name
))

length(shared_genes)
head(shared_genes, 20)

## 12,603 shared genes to compare
ab_shared <- ab_named[ab_named$gene_name %in% shared_genes, ]
mz_shared <- mz_named[mz_named$gene_name %in% shared_genes, ]
nb_shared <- nb_named[nb_named$gene_name %in% shared_genes, ]
on_shared <- on_named[on_named$gene_name %in% shared_genes, ]
pn_shared <- pn_named[pn_named$gene_name %in% shared_genes, ]

##remove extra columns
ab_shared <- ab_shared[, c("gene_name", "Ab_5_Eye", "Ab_6_Eye")]
mz_shared <- mz_shared[, c("gene_name", "Mz_1_Eye", "Mz_2_Eye")]
nb_shared <- nb_shared[, c("gene_name", "Nb_4_Eye", "Nb_5_Eye")]
on_shared <- on_shared[, c("gene_name", "On_1_Eye", "On_2_Eye", "On_3_Eye")]
pn_shared <- pn_shared[, c("gene_name", "Pn_2_Eye", "Pn_3_Eye", "Pn_4_Eye")]

##remove duplicate gene names ready for deseq2
ab_shared <- ab_shared[!duplicated(ab_shared$gene_name), ]
mz_shared <- mz_shared[!duplicated(mz_shared$gene_name), ]
nb_shared <- nb_shared[!duplicated(nb_shared$gene_name), ]
on_shared <- on_shared[!duplicated(on_shared$gene_name), ]
pn_shared <- pn_shared[!duplicated(pn_shared$gene_name), ]

##check
sum(duplicated(ab_shared$gene_name))
sum(duplicated(mz_shared$gene_name))
sum(duplicated(nb_shared$gene_name))
sum(duplicated(on_shared$gene_name))
sum(duplicated(pn_shared$gene_name))
## all zero

##gene names into row names
rownames(ab_shared) <- ab_shared$gene_name
rownames(mz_shared) <- mz_shared$gene_name
rownames(nb_shared) <- nb_shared$gene_name
rownames(on_shared) <- on_shared$gene_name
rownames(pn_shared) <- pn_shared$gene_name

ab_shared <- ab_shared[, -1]
mz_shared <- mz_shared[, -1]
nb_shared <- nb_shared[, -1]
on_shared <- on_shared[, -1]
pn_shared <- pn_shared[, -1]

##reorder rows
ab_shared <- ab_shared[shared_genes, ]
mz_shared <- mz_shared[shared_genes, ]
nb_shared <- nb_shared[shared_genes, ]
on_shared <- on_shared[shared_genes, ]
pn_shared <- pn_shared[shared_genes, ]

all(rownames(ab_shared) == shared_genes)
all(rownames(mz_shared) == shared_genes)
all(rownames(nb_shared) == shared_genes)
all(rownames(on_shared) == shared_genes)
all(rownames(pn_shared) == shared_genes)
## all true

##combine into one matrix
count_mat <- cbind(ab_shared, mz_shared, nb_shared, on_shared, pn_shared)
dim(count_mat)
head(count_mat)

##check with metadata
##check with metadata
colnames(count_mat)
rownames(coldata)
coldata <- coldata[colnames(count_mat), ]

all(colnames(count_mat) == rownames(coldata))
## true

head(coldata)


## STEP 2 - LOW-COUNT FILTERING
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

