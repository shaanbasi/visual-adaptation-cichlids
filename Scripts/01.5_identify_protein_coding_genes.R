## Download Ensembl release 101 GTF annotations
dir.create(
  "Data/GTF_Ensembl101",
  recursive = TRUE,
  showWarnings = FALSE
)

gtf_urls <- c(
  Ab = paste0(
    "https://ftp.ensembl.org/pub/release-101/gtf/",
    "haplochromis_burtoni/",
    "Haplochromis_burtoni.AstBur1.0.101.gtf.gz"
  ),
  
  Mz = paste0(
    "https://ftp.ensembl.org/pub/release-101/gtf/",
    "maylandia_zebra/",
    "Maylandia_zebra.M_zebra_UMD2a.101.gtf.gz"
  ),
  
  Nb = paste0(
    "https://ftp.ensembl.org/pub/release-101/gtf/",
    "neolamprologus_brichardi/",
    "Neolamprologus_brichardi.NeoBri1.0.101.gtf.gz"
  ),
  
  On = paste0(
    "https://ftp.ensembl.org/pub/release-101/gtf/",
    "oreochromis_niloticus/",
    "Oreochromis_niloticus.O_niloticus_UMD_NMBU.101.gtf.gz"
  ),
  
  Pn = paste0(
    "https://ftp.ensembl.org/pub/release-101/gtf/",
    "pundamilia_nyererei/",
    "Pundamilia_nyererei.PunNye1.0.101.gtf.gz"
  )
)

## Keep the species names when creating local file paths
gtf_files <- setNames(
  file.path(
    "Data/GTF_Ensembl101",
    basename(gtf_urls)
  ),
  names(gtf_urls)
)

## Download
for (sp in names(gtf_urls)) {
  
  message("Downloading ", sp, "...")
  
  download.file(
    url = gtf_urls[[sp]],
    destfile = gtf_files[[sp]],
    mode = "wb",
    method = "libcurl"
  )
}

## Check names are present
print(gtf_files)

## Check that all five downloaded
data.frame(
  species = names(gtf_files),
  file = unname(gtf_files),
  exists = file.exists(gtf_files),
  size_MB = round(
    file.info(gtf_files)$size / 1024^2,
    2
  )
)




## EXTRACTION
library(dplyr)
library(stringr)
library(readr)
library(tibble)

extract_gtf_genes <- function(gtf_file) {
  
  message("Reading: ", basename(gtf_file))
  
  gtf <- read_tsv(
    gtf_file,
    comment = "#",
    col_names = c(
      "seqname",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attribute"
    ),
    col_types = cols(.default = "c"),
    progress = FALSE
  )
  
  genes <- gtf %>%
    filter(feature == "gene") %>%
    transmute(
      gene_id = str_match(
        attribute,
        'gene_id "([^"]+)"'
      )[, 2],
      
      gene_name = str_match(
        attribute,
        'gene_name "([^"]+)"'
      )[, 2],
      
      gene_biotype = coalesce(
        str_match(
          attribute,
          'gene_biotype "([^"]+)"'
        )[, 2],
        
        str_match(
          attribute,
          'gene_type "([^"]+)"'
        )[, 2]
      )
    ) %>%
    distinct()
  
  genes
}

gtf_annotations <- lapply(
  gtf_files,
  extract_gtf_genes
)

names(gtf_annotations) <- names(gtf_files)

## Check how many genes were extracted
sapply(gtf_annotations, nrow)

## Check biotype counts
for (sp in names(gtf_annotations)) {
  
  cat("\n\n", sp, "\n")
  
  print(
    gtf_annotations[[sp]] %>%
      dplyr::count(
        gene_biotype,
        sort = TRUE
      ),
    n = 20
  )
}



## MATCHING

bed_files <- c(
  Ab = "Data/Bed/Haplochromis_burtoni.AstBur1.0.101.bed",
  Mz = "Data/Bed/Maylandia_zebra.M_zebra_UMD2a.101.bed",
  Nb = "Data/Bed/Neolamprologus_brichardi.NeoBri1.0.101.bed",
  On = "Data/Bed/Oreochromis_niloticus.O_niloticus_UMD_NMBU.101.bed",
  Pn = "Data/Bed/Pundamilia_nyererei.PunNye1.0.101.bed"
)

count_mat_filt <- readRDS(
  "Data/Processed/count_mat_filt.rds"
)

final_gene_names <- rownames(count_mat_filt)

cat(
  "Genes in final matrix:",
  length(final_gene_names),
  "\n"
)

read_bed_gene_map <- function(bed_file) {
  
  bed <- read.delim(
    bed_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  bed <- bed[, c(7, 8)]
  
  colnames(bed) <- c(
    "gene_id",
    "gene_name"
  )
  
  bed %>%
    filter(
      !is.na(gene_id),
      gene_id != "",
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    distinct()
}

biotype_results <- lapply(
  names(bed_files),
  function(sp) {
    
    message("Matching ", sp, "...")
    
    bed_map <- read_bed_gene_map(
      bed_files[[sp]]
    )
    
    ann <- gtf_annotations[[sp]]
    
    matched <- bed_map %>%
      left_join(
        ann %>%
          select(
            gene_id,
            gene_biotype
          ),
        by = "gene_id"
      ) %>%
      filter(
        gene_name %in% final_gene_names
      )
    
    by_gene_name <- matched %>%
      group_by(gene_name) %>%
      summarise(
        has_annotation =
          any(!is.na(gene_biotype)),
        
        protein_coding =
          any(
            gene_biotype == "protein_coding",
            na.rm = TRUE
          ),
        
        biotypes = paste(
          sort(
            unique(
              gene_biotype[
                !is.na(gene_biotype)
              ]
            )
          ),
          collapse = "; "
        ),
        
        .groups = "drop"
      )
    
    by_gene_name <- tibble(
      gene_name = final_gene_names
    ) %>%
      left_join(
        by_gene_name,
        by = "gene_name"
      ) %>%
      mutate(
        has_annotation =
          tidyr::replace_na(
            has_annotation,
            FALSE
          ),
        
        protein_coding =
          tidyr::replace_na(
            protein_coding,
            FALSE
          )
      )
    
    list(
      matched = matched,
      by_gene_name = by_gene_name
    )
  }
)

names(biotype_results) <- names(bed_files)

protein_coding_summary <- bind_rows(
  lapply(
    names(biotype_results),
    function(sp) {
      
      x <- biotype_results[[sp]]$by_gene_name
      
      tibble(
        species = sp,
        final_genes = nrow(x),
        genes_with_biotype =
          sum(x$has_annotation),
        protein_coding =
          sum(x$protein_coding),
        annotated_non_protein_coding =
          sum(
            x$has_annotation &
              !x$protein_coding
          ),
        unmatched =
          sum(!x$has_annotation)
      )
    }
  )
)

print(
  protein_coding_summary,
  n = Inf
)


## Show non-protein coding genes for each species
noncoding_genes <- bind_rows(
  lapply(
    names(biotype_results),
    function(sp) {
      
      biotype_results[[sp]]$by_gene_name %>%
        filter(
          has_annotation,
          !protein_coding
        ) %>%
        mutate(
          species = sp
        ) %>%
        select(
          species,
          gene_name,
          biotypes
        )
    }
  )
)

print(
  noncoding_genes,
  n = Inf
)


## Confirm which non-protein-coding gene names occur in all species
noncoding_shared <- noncoding_genes %>%
  dplyr::count(gene_name) %>%
  dplyr::filter(
    n == length(unique(noncoding_genes$species))
  )

print(noncoding_shared)



## Remove non-coding from dataset ready for DESeq2
noncoding_to_remove <- c(
  "5S_rRNA",
  "Metazoa_SRP",
  "U6"
)

count_mat_pc <- count_mat_filt[
  !rownames(count_mat_filt) %in% noncoding_to_remove,
]

dim(count_mat_pc)


intersect(
  rownames(count_mat_pc),
  noncoding_to_remove
)

saveRDS(
  count_mat_pc,
  "Data/Processed/count_mat_filt.rds"
)
