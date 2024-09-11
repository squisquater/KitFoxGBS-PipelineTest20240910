# Population Genomics and Landscape Genetics Pipeline

## Overview

This Snakemake workflow automates the process of aligning genomic reads, identifying SNPs, performing population structure analysis, and conducting landscape genetics analysis. It integrates tools like `bwa`, `stacks`, `bcftools`, `PLINK`, `ADMIXTURE`, and custom scripts for various stages of data processing and analysis. The workflow is modular and can be run in parts or as a complete pipeline.

## Key Features

- **Read Alignment and Quality Control**: Align paired-end reads to a reference genome, add read group information, and generate quality control metrics.
- **SNP Calling**: Use the Stacks `ref_map` and `populations` pipelines to call SNPs.
- **Population Genetics**: Perform PCA, MDS, and ADMIXTURE analyses to estimate population structure.
- **Landscape Genetics**: Generate pairwise FST and Nei's D matrices, calculate centroids, and test isolation-by-distance using Mantel tests.
- **Modular Workflow**: Individual snakefiles allow for running different sections independently, but the entire pipeline can also be run from the master snakefile.

## Directory Structure

```plaintext
├── 00.snakefile_master                 # Master Snakemake file that coordinates the entire workflow
├── accessory_scripts                   # Custom scripts for various steps in the pipeline
│   ├── GenDist.R                       # Script for generating distance matrices and heatmaps
│   ├── admixture_metrics.R             # Script for summarizing ADMIXTURE CV errors
│   ├── plot_admixture_barplots.R       # Script for plotting ADMIXTURE results
│   └── (other R/Python scripts)        # Additional scripts for specific tasks (e.g., isolation by distance, population centroids)
├── config                              # Configuration files for each module
│   ├── 01.snakefile_alignPE.yml
│   ├── 02.snakefile_stacks.yml
│   ├── 04.snakefile_admixture.yml
│   ├── 06.snakefile_landgen.yml
│   └── snakemake.yml                   # Snakemake environment configuration
├── envs                                # Conda environment files
│   ├── GBS.yml
│   ├── LandGen.yml
│   └── (other .yml files)
├── input_files                         # Input files for analysis
│   ├── SJKF_cleaned_database_n484.txt  # Master database containing sample information
│   ├── hexgrid5km.shp                  # Hexagonal grid shapefiles for landscape genetics (can be customized)
│   └── (additional input files)
├── modules                             # Subdirectories containing Snakemake workflows for each analysis module
│   ├── alignPE
│   ├── stacks
│   ├── admixture
│   └── landgen
