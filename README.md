# Gene Expression Landscapes Underlying Visual Adaptation in East African Cichlids

## Overview
This repository contains a comparative transcriptomic analysis investigating how gene expression changes drive visual system adaptation across five East African cichlid species occupying distinct light environments.

## Project Aim
To characterise gene expression changes in eye tissue across cichlid species adapted to divergent ecological niches, uncovering the genetic basis for rapid visual system evolution and ecological divergence.

## Species Studied
- *Maylandia zebra* (M. zebra)
- *Oreochromis niloticus* (O. niloticus)
- *Pundamilia nyererei* (P. nyererei)
- *Astatotilapia burtoni* (A. burtoni)
- *Neolamprologus brichardi* (N. brichardi)

## Analysis Pipeline
1. **Quality Control** – FastQC & Trim Galore!
2. **Read Mapping** – HISAT2 to species-specific reference genomes
3. **Quantification** – HTSeq transcript abundance quantification
4. **Differential Expression** – DESeq2 (FDR < 0.05)
5. **Functional Enrichment** – Gene Ontology analysis via g:Profiler
6. **Co-expression Networks** – WGCNA to identify gene modules

## Key Findings


## Contents
- `data/` – Raw count tables and processed expression matrices
- `scripts/` – R scripts for DESeq2, GO enrichment, and WGCNA analyses
- `results/` – Figures, tables, and gene module assignments
- `docs/` – Methods flowcharts and supplementary information
- `references/` – Reference list
- `GAI declaration/` – Methods flowcharts and supplementary information

## Supervisor
Dr Tarang Mehta, Department of Biochemistry, Cell and Systems Biology, University of Liverpool
