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
- `16S prevalence filtering.R`- script that prevalence filters the 16S phyloseq object at both 15% & 30% for use in machine learning methods.
- `MintTea pipeline.R`- script that both runs a grid search using multiomics data to 1) find the optimal parameters 2) runs the  MintTea analysis pipline using those settings.
- `MintTea external validation.R`- script that runs external cross validation of the MintTea analysis.
- `Out-of-fold MintTea boostrap.R`- script that applies bootstrapping validation to the external cross validation of the MintTea models.
- `Aligned blocks generation.R`- script that produces the aligned omics blocks used in the Multinomial elastic net (glmnet) pipeline.
- `Glmnet pipeline.R`- script that runs the glmnet pipeline on each omic separately, then fuse predicted class probabilities.

## Requirements & Licenses

The code within this repository is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0) (see the `LICENSE` file).

The scripts rely on several external tools and libraries, each with its own license:

### R

[R version 4.4.0](https://www.r-project.org/) was used for many scripts in this repository. Packages utilised include:

- `Phyloseq` 1.48.0 [AGPL-3](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
- `Microbiome` 1.26.0 [BSD_2](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html)
- `Picante` 1.8.2 [GPL-2](https://cran.r-project.org/web/packages/picante/index.html)
- `lme4` 1.1.35.5 [GPL-2/GPL-3](https://cran.r-project.org/web/packages/lme4/index.html)
- `lmerTest` 3.1.3 [GPL-2/GPL-3](https://cran.r-project.org/web/packages/lmerTest/index.html)
- `dplyr` 1.1.4 [MIT](https://cran.r-project.org/web/packages/dplyr/index.html)
- `vegan` 2.6.6.1 [GPL-2](https://cran.r-project.org/web/packages/vegan/index.html)
- `compositions` 2.0.8 [GPL-2/GPL-3](https://cran.r-project.org/web/packages/compositions/index.html)
- `pairwiseAdonis` 0.4.1 [GPL-3](https://github.com/pmartinezarbizu/pairwiseAdonis)
- `tidyr` 1.3.1 [MIT](https://cran.r-project.org/web/packages/tidyr/index.html)
- `broom` 1.0.8 [MIT](https://cran.r-project.org/web/packages/broom/index.html)
- `ANCOMBC` 2.6.0 [Artistic-2.0](https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html)
- `MintTea` 1.0.0 [MIT](https://github.com/efratmuller/MintTea)
- `glmnet` 4.1.10 [GPL-2](https://cran.r-project.org/web/packages/glmnet/index.html)

### FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used for checking for the presence of adapter sequences. (License: [GPL v3.0](https://github.com/s-andrews/FastQC))

### Trimmomatic

[Trimmomatic](https://github.com/usadellab/Trimmomatic/releases) was used for trimming adapter sequences from raw reads. (License: [GPL v3.0](https://github.com/usadellab/Trimmomatic/blob/main/LICENSE))

### DADA2

[DADA2](https://benjjneb.github.io/dada2/tutorial.html) 1.36.0 was used for processing raw 16S reads into ASVs and abundances. (License: [LGPL-3.0](https://github.com/benjjneb/dada2/blob/master/LICENSE))

## Authors

- **Peter Prendergast**: peter.prendergast@pg.canterbury.ac.nz
