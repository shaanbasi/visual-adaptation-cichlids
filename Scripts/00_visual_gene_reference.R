## VISUAL-SYSTEM CANDIDATE GENE REFERENCE

visual_genes <- c(
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
)

visual_genes <- unique(
  toupper(
    trimws(
      visual_genes
    )
  )
)

dir.create(
  "Data/Reference",
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  tibble::tibble(
    gene_name = visual_genes
  ),
  "Data/Reference/visual_candidate_genes.csv"
)

