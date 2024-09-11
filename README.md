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
```
## Getting Started

1. Clone the Repository

First, clone the repository from GitHub to your local machine:

git clone https://github.com/<your-username>/<your-repository>.git
cd <your-repository>

Replace <your-username> and <your-repository> with your actual GitHub username and repository name.

2. Set Up Conda Environment

Once the repository is cloned, create and activate the Snakemake environment:

conda env create -f envs/snakemake.yml
conda activate snakemake

The pipeline will automatically create other necessary Conda environments as it runs, but they are defined in the envs/ directory (e.g., GBS.yml, LandGen.yml).

3. Customize Input Files

Place your input files in the input_files/ directory. Key files include:

- Master Database: This file (SJKF_cleaned_database_n484.txt) contains sample information like sample IDs and library information.
- Shapefiles: Hexagonal grid shapefiles for landscape genetics (hexgrid5km.shp, hexgrid10km.shp). If your project spans outside California or uses a different region, replace these shapefiles with your own.

4. Modify Configuration Files

Each module of the pipeline has a corresponding configuration file in the config/ directory. These files define paths to input files, output directories, and specific parameters for each analysis.

Steps to Modify:

Open the relevant .yml file in the config/ directory using a text editor.
Modify the paths to your input files, set parameters like depth thresholds, and update other analysis options as needed.

Key Files:

01.snakefile_alignPE.yml: Defines paths for read alignment, reference genome, and sample metadata.
02.snakefile_stacks.yml: Defines settings for SNP calling, population mapping, and filtering.
04.snakefile_admixture.yml: Defines the K-values range for ADMIXTURE analysis.
06.snakefile_landgen.yml: Defines settings for landscape genetics, including IBD analysis.

Example of modifying a configuration file:

# Example: config/01.snakefile_alignPE.yml
home_dir: "/path/to/your/project"
data_dir: "/path/to/your/data"
reference_file: "/path/to/reference_genome.fna"
masterDB: "SJKF_cleaned_database_n484.txt"
depth_threshold: 15

Make sure all file paths and parameters match your data and desired analysis.

5. Run the Pipeline

After customizing the input files and configuration files, you can run the entire workflow or individual modules.

Run the entire workflow:

snakemake --use-conda --cores <number_of_cores>

Run individual modules:

To run a specific analysis (e.g., read alignment, SNP calling), use the corresponding snakefile:

snakemake -s modules/alignPE/01.snakefile_alignPE --use-conda --cores <number_of_cores>

## Output Files

After running the workflow, the output will include:

- BAM files and alignment metrics.
- VCF files containing SNPs.
- PLINK files for population genetics analysis.
- ADMIXTURE results including Q-value summaries and bar plots.
- Distance matrices (FST and Nei’s D) and heatmaps for landscape genetics.
- Isolation-by-distance plots from the Mantel test.

Output directories and files are defined in the config files.

## Troubleshooting

1. Ensure that input file paths are correctly specified in the configuration files.
2. Activate the appropriate Conda environment before running the pipeline.
3. Check Snakemake logs for error messages if the workflow fails.

## Acknowledgements

This pipeline integrates multiple tools and custom scripts. If you use this workflow in your research, please cite the following:

- Snakemake: https://snakemake.readthedocs.io
- ADMIXTURE: https://github.com/chrchang/admixture
- Stacks: https://catchenlab.life.illinois.edu/stacks/
- PLINK: https://www.cog-genomics.org/plink/
- BCFtools: http://samtools.github.io/bcftools/bcftools.html
