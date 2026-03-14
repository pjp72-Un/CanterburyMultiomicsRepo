# CanterburyMultiomicsRepo

This repository contains scripts supporting the multiomics Canterbury, New Zealand analysis of the gut microbiome and metabolome in children with coeliac disease.

## About

### Publication

In this work we analysed both 16S rRNA sequencing data and targeted and discovery mode metabolomics in children in Canterbury, New Zealand. Analysis involved bacterial diveristy, differential abundance, gas chromatography mass spectrometry (GC-MS), nuclear magnetic resonance (NMR), and machine learning multiomics (MintTea & Glmnet). 

## Contents

### /DADA2

- `Dada2 pipline.R`- script containing the deonising pipeline taking raw reads to the raw phyloseq object.
- `Phyloseq processing.R`- script taking the raw phyloseq object and filtering and processing it into group based phyloseq objects for downstream analysis.

### /16S

- `Alpha diversity paired.R`- script that runs alpha diversity analysis for paired untreated-treated samples.
- `Alpha diversity.R`- script that runs alpha diversity analysis for non-paired comparisons.
- `Faith's PD.R`- script that runs phylogenetic analysis across all comparisons.
- `Beta diversity.R`- script that runs Aitchison distance beta diveristy across all comparisons.
- `Phylum abundance analysis.R`- script that runs phylum relative abundance analysis across the top 10 phyla for all comparisons.
- `Differential abundance.R`- script that runs ANCOMBC2 differential abundance analysis for non-paired comparisons.
- `Differential abundance paired.R`- script that runs ANCOMBC2 differential abundance analysis for paired (untreated-treated) comparisons.

### /Metabolomics

- `NMR lm.R`- script that runs LM analysis for the identified NMR metabolite bins across all comparisons.
- `GCMS lm.R`- script that runs LM analysis for all targeted GC-MS metabolites across all comparisons.

### /Machine Learning

- `Metabolome merge.R`- script that merges, filters, and median imputes log2 transformed metabolites using the GC-MS metabolites for any overlapping metabolites across methods.
- `16S prevalence filtering`- script that prevalence filters the 16S phyloseq object at both 15% & 30% for use in machine learning methods.
- `MintTea pipeline.R`- script that both runs a grid search using multiomics data to 1) find the optimial parameters 2) runs the  MintTea analysis pipline using those settings.
- `MintTea external validation.R`- script that runs external cross validation of the MintTea analysis.
- `Aligned omics generation.R`- script that produces the aligned omics blocks used in the Multinomial elastic net (glmnet) pipeline.
- `Glmnet pipeline.R`- script that runs the glmnet pipeline on each omic separately, then fuse predicted class probabilities.

## Requirements & Licenses

The code within this repository is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0) (see the `LICENSE` file).

The scripts rely on several external tools and libraries, each with its own license:

### R

[R version 4.5.1](https://www.r-project.org/) was used for many scripts in this repository. Packages utilised include:
