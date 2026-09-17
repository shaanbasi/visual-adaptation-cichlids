# Retinal Transcriptomic Divergence in East African Cichlids

This repository contains the analysis associated with an MSc
Bioinformatics dissertation investigating retinal gene-expression
divergence across five East African cichlid species. The study combines
differential expression, Gene Ontology (GO) enrichment, weighted gene
co-expression network analysis (WGCNA), and targeted analysis of
visual-system genes to examine transcriptional differences within and
beyond canonical opsin pathways.

## Study design

Retinal RNA-seq data were analysed from 12 individuals representing five
cichlid species:

  Code   Species                        n
  ------ ---------------------------- ---
  Ab     *Astatotilapia burtoni*        2
  Mz     *Maylandia zebra*              2
  Nb     *Neolamprologus brichardi*     2
  On     *Oreochromis niloticus*        3
  Pn     *Pundamilia nyererei*          3

Species-specific read processing and quantification were performed prior
to the analyses contained in this repository. Reads were quality checked
using FastQC, trimmed with Trim Galore!, aligned to species-specific
reference genomes using HISAT2, and quantified using HTSeq.

The downstream analysis was restricted to 8,522 protein-coding genes
shared across all five species.

## Analysis workflow

1.  **Data preparation** - preparation and filtering of count matrices
    and sample metadata.
2.  **Protein-coding gene filtering** - identification and retention of
    shared protein-coding genes using Ensembl release 101 annotations.
3.  **Differential expression** - DESeq2 analysis, with the primary
    analysis comparing each species against the remaining four species
    using species-vs-rest contrasts.
4.  **Functional enrichment** - GO enrichment of genes showing
    significantly higher or lower expression, followed by semantic
    reduction of related GO terms.
5.  **Co-expression analysis** - WGCNA of the pooled retinal expression
    dataset to identify co-expression modules and module-species
    relationships.
6.  **Hub-gene analysis** - assessment of module membership and gene
    significance within supported module-species relationships.
7.  **Targeted visual-system analysis** - evaluation of 38 curated genes
    involved in opsins, phototransduction, photoreceptor recovery, the
    visual cycle, and retinal regulation.
8.  **Integration** - integration of differential-expression and WGCNA
    evidence to identify visual-system genes supported across multiple
    analytical layers.

## Repository structure

``` text
visual-adaptation-cichlids/
├── Data/
├── Scripts/
├── Results/
├── Figures/
├── Appendix/
├── README.md
└── visual-adaptation-cichlids.Rproj
```

-   `Data/` - input, reference and processed analysis data.
-   `Scripts/` - R scripts for the downstream analysis.
-   `Results/` - numerical outputs from DESeq2, GO and WGCNA analyses.
-   `Figures/` - main and supporting figures.
-   `Appendix/` - outputs used in the dissertation appendices.

### Data

The `Data/` directory contains HTSeq count data, sample metadata, BED
annotations, the curated visual-system gene reference, and processed
objects used by downstream analyses.

Ensembl release 101 GTF files are required for the protein-coding gene
filtering step but are excluded from the repository because they are
externally available reference annotations.

### Scripts

Scripts are numbered broadly in analysis order:

-   `00_visual_gene_reference.R` - visual-system gene reference setup.
-   `01_prepare_counts_and_metadata.R` - count and metadata preparation.
-   `01.5_identify_protein_coding_genes.R` - protein-coding gene
    identification and filtering.
-   `02_run_deseq2.R` - differential-expression analysis.
-   `03_run_go_enrichment.R` - GO enrichment analysis.
-   `04_GO_plots.R` - GO processing, semantic reduction and
    visualisation.
-   `05_wgcna_prep.R` - WGCNA input preparation.
-   `06_run_wgcna.R` - co-expression network construction and
    module-species analysis.
-   `07_identify_hub_genes.R` - module membership, gene significance and
    supported hub-gene analysis.
-   `08_WGCNA_module_GO_enrichment.R` - functional enrichment of
    selected WGCNA modules.
-   `09_WGCNA_module_expression.R` - module eigengene expression
    analysis.
-   `10_assess_WGCNA_module_preservation.R` - exploratory
    module-preservation assessment.
-   `11_targeted_visual_gene_analysis.R` - targeted analysis of curated
    visual-system genes.
-   `12_plot_targeted_visual_genes.R` - visual-system gene-expression
    visualisation.
-   `13_integrate_visual_genes_WGCNA.R` - integration of visual-gene DE
    and WGCNA evidence.
-   `14_plot_integrated_visual_candidates.R` - integrated
    visual-candidate visualisation.
-   `15_final_summaries_for_results.R` - final result summaries.

## Main outputs

-   `Results/DESeq2/species_vs_rest/` - species-vs-rest
    differential-expression results.
-   `Results/GO/species_vs_rest/` - GO enrichment and semantic-reduction
    results.
-   `Results/WGCNA/` - module assignments, module-species associations
    and hub-gene analyses.
-   `Results/targeted_visual_genes/` - targeted visual-system gene
    results and integrated DE/WGCNA evidence.
-   `Figures/` - final and supporting figures.
-   `Appendix/` - tables and diagnostic figures used as dissertation
    appendices.

## Summary of findings

Retinal gene expression showed substantial interspecific structure
across the five cichlid species. Differential expression within the
curated visual-system gene set extended beyond opsins to genes involved
in phototransduction, photoreceptor recovery, retinoid processing and
retinal regulation.

GO enrichment also identified differences in broader cellular processes,
including overlapping translation-related genes showing lower expression
in *A. burtoni* and higher expression in *O. niloticus*.

WGCNA identified two module-species relationships that remained
significant after multiple-testing correction, associated with *O.
niloticus* and *P. nyererei*. Integration of differential-expression and
co-expression evidence highlighted *rdh12* and *nr2e3* in *O. niloticus*
as candidates for further investigation.

These analyses identify patterns of retinal transcriptional divergence
among species but do not, by themselves, establish that observed
expression differences are adaptive.

## Software

The final downstream analyses were performed in R 4.5.3, principally
using DESeq2 1.50.2, WGCNA 1.74, gprofiler2 0.2.4, simplifyEnrichment
2.4.1, and ashr 2.2.63. Additional R packages used for data manipulation
and visualisation are specified within the individual scripts.

## Reproducibility and limitations

The study contains relatively small numbers of biological replicates per
species (n = 2-3). The WGCNA analysis therefore provides exploratory
evidence of co-expression structure rather than formal evidence of
network preservation across species.

Because reads were aligned to species-specific genome assemblies and
annotations, cross-species differences in genome and annotation quality
may also contribute to observed expression differences. The analysis
used whole-retina RNA-seq, so expression differences may reflect both
transcriptional regulation and variation in retinal cell composition.

## Project information

This analysis was conducted as part of an MSc Bioinformatics
dissertation at the University of Liverpool (2026).

**Supervisor:** Dr Tarang Mehta, Department of Biochemistry, Cell and
Systems Biology, University of Liverpool.
