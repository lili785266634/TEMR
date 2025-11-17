# TEMR Analysis Code for “Trans-ethnic Conditional-Likelihood Mendelian Randomization Reveals Novel Causal Links Between Depression and Allergic Rhinitis/Asthma Across Five Global Ancestries”

## Overview

This repository contains reproducible R code for the trans-ethnic conditional-likelihood Mendelian randomization (TEMR) analyses described in the article:

> “Trans-ethnic Conditional-Likelihood Mendelian Randomization Reveals Novel Causal Links Between Depression and Allergic Rhinitis/Asthma Across Five Global Ancestries”

The script implements the forward-direction TEMR analysis of **major depressive disorder (MDD) on asthma** as a worked example. All other analyses in the paper (e.g., MDD ↔ allergic rhinitis/asthma in different ancestries and directions) can be reproduced by analogy using the same workflow.

## Data Availability

All GWAS summary statistics used in this project are **publicly available** and can be downloaded according to the **Data Availability Statement** in the original article.  
Please refer to the paper for:

- The names of the GWAS consortia and datasets  
- Download links or accession IDs  
- Ancestry-specific summary statistics used in the TEMR framework  

This repository does **not** redistribute the original GWAS summary data; instead, it provides the analysis pipeline and example formatted input data.

## Repository Contents

- **`analysis_code.R`**  
  Main R script for reproducing the TEMR analysis.  
  Compared with the original TEMR code, this script:

  - Adds more detailed implementation steps  
  - Includes extensive in-line comments  
  - Explicitly shows data import and preprocessing  
  - Incorporates additional quality control procedures  

  The goal is to make it easier for users to understand and adapt the TEMR workflow for their own traits and datasets.

- **`MDD_asthma_156.csv`**  
  Example input file for the **MDD → asthma** TEMR analysis in the **EUR** ancestry.  
  This file represents the final set of sentinel SNPs used as instruments for MDD in a conventional two-sample MR framework after:

  - Applying the genome-wide significance threshold (P < 5 × 10⁻⁸)  
  - Linkage disequilibrium (LD) clumping  
  - Allele harmonization between exposure and outcome datasets  
  - Steiger filtering and other standard MR QC steps  and etc

  These 156 SNPs are then used as the input instruments for the conditional-likelihood TEMR analysis.

## How to Use

1. **Prepare the data**

   - Download the required GWAS summary statistics following the Data Availability Statement in the article.  
   - Format your exposure and outcome datasets to match the structure illustrated by `MDD_asthma_156.csv`.

2. **Run the example analysis**

   - Open `analysis_code.R` in R or RStudio.  
   - Set your working directory to the folder containing `analysis_code.R` and `MDD_asthma_156.csv`.  
   - Install any required R packages (see the top of the script) if needed.  
   - Run the script line by line or source it to reproduce the MDD → asthma TEMR results for the EUR ancestry.

3. **Extend to other analyses**
 
   - Follow the same workflow as in `analysis_code.R`, modifying file paths and parameter settings where appropriate.

## Citation

If you use this code in your work, please cite the original article:

> “TEMR: Trans-ethnic mendelian randomization method using large-scale GWAS summary datasets, The American Journal of Human Genetics 112 (2025) 28–43. https://doi.org/10.1016/j.ajhg.2024.11.006.” (full citation details to be added here once available).

