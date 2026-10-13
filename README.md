# anastatica-population-transcriptomics
# Scripts for data analysis

This repository contains the scripts used to process and analyse the data presented in the following manuscript:

**Manuscript citation:**  
*[To be added once the citation is available.]*

## Repository overview

The repository is organised into six folders, corresponding to the main stages of the analytical workflow.

| Folder | Description |
|---|---|
| `01_climate/` | Extraction, correlation analysis, and principal component analysis (PCA) of bioclimatic variables. |
| `02_rnaseq_data_processing/` | Preparation of sample metadata and RNA-seq data processing. |
| `03_RDA/` | Redundancy analysis (RDA) to investigate isolation by distance (IBD), isolation by environment (IBE), and east–west differentiation, including functional enrichment analyses and visualisation. |
| `04_models/` | Single, combined, and interaction models, followed by post-processing, quality control, and diagnostic analyses. |
| `05_gene_classification/` | Classification of genes into four classes (climate-only, range-only, climate-and-range, and non-associated) and visualisation of gene classes. |
| `06_functional_enrichment/` | Functional overrepresentation analysis, Gene Ontology (GO) visualisation, and gene set enrichment analysis (GSEA) using `zenith`. |

## Data and metadata

The metadata required for the analyses are provided in the **Supplementary Materials of the associated manuscript**.

## Analysis workflow

The scripts are grouped by analytical objective. The folder numbering provides a general guide to the workflow, although some analyses may be run independently once their required input files are available.

### 1. Climate data analysis — `01_climate/`

Scripts for extracting bioclimatic variables from WorldClim, assessing correlations among variables, and performing PCA.

- `01_extract_worldclim_bioclim.py` — Extracts bioclimatic variables from WorldClim.
- `02_bioclimatic_correlations.r` — Assesses correlations among bioclimatic variables.
- `03_PCA_bioclimatic_variables.r` — Performs PCA of bioclimatic variables.
- `04_PCA_scores_loadings_plots.R` — Generates plots of PCA scores and loadings.

### 2. RNA-seq data processing — `02_rnaseq_data_processing/`

Scripts for preparing the metadata and processing RNA-seq data.

- `04_make_metadata_and_run_tximeta.r` — Prepares metadata and runs `tximeta` for transcript quantification data handling.

### 3. Redundancy analysis — `03_RDA/`

Scripts for investigating patterns of transcriptomic variation in relation to geographic and environmental factors, including IBD, IBE, and east–west differentiation.

- `01_IBD_IBE_transcriptome_RDA.R` — Performs transcriptome-level RDA analyses of IBD and IBE.
- `01_IBD_IBE_transcriptome_RDA_run.sh` — Shell script to run the corresponding R analysis.
- `02_IBD_IBE_transcriptome_RDA_East_West.r` — Performs east–west comparisons using RDA.
- `02_IBD_IBE_transcriptome_RDA_East_West_run.sh` — Shell script to run the corresponding R analysis.
- `03_RDA_functional_enrichment_East_West.R` — Performs functional enrichment analyses associated with east–west RDA results.
- `03_RDA_functional_enrichment_East_West_run.sh` — Shell script to run the corresponding R analysis.
- `04_RDA_plot.R` — Generates RDA visualisations.
- `04_RDA_plot_run.sh` — Shell script to run the plotting script.

### 4. Statistical models — `04_models/`

Scripts for fitting models for variance partitioning and differential expression analysis, and evaluating their results.

- `01_single_models.r` — Fits single-predictor models.
- `02_combined_models.r` — Fits combined-predictor models.
- `03_interaction_models.r` — Fits models incorporating interactions.
- `04_models_postprocessing_qc_diagnostics.r` — Performs model post-processing, quality control, and diagnostic analyses.
- `01_models_run.bash` and `04_models_postprocessing_qc_diagnostics_run.sh` — Shell scripts for running the corresponding analyses.

### 5. Gene classification — `05_gene_classification/`

Scripts for classifying genes into climate-only, range-only, climate-and-range, and non-associated, and visualising the resulting classes.

- `01_classify_genes_variance_by_class.R` — Classifies genes.
- `01_classify_genes_variance_by_class_run.sh` — Shell script to run the classification analysis.
- `02_gene_classes_plot.R` — Generates plots of gene classifications.
- `02_gene_classes_plot_run.sh` — Shell script to run the plotting script.

### 6. Functional enrichment — `06_functional_enrichment/`

Scripts for identifying and visualising functional categories associated with the analysed gene sets.

- `01_overrepresentation.R` — Performs functional overrepresentation analysis.
- `01_overrepresentation_run.sh` — Shell script to run the overrepresentation analysis.
- `02_GO_dotplot.R` — Generates Gene Ontology (GO) dot plots.
- `03_zenith_GSEA.R` — Performs gene set enrichment analysis using `zenith`.
- `03_zenith_GSEA_run.sh` — Shell script to run the GSEA analysis.

## Requirements

The analyses use a combination of Python, R, and shell scripts.

The required software, packages, and dependencies vary between analysis steps. These may include:

- **Python** for climate data extraction.
- **R** and the relevant R/Bioconductor packages for statistical analyses, RNA-seq data processing, visualisation, and functional enrichment.
- **Bash** for running the provided shell scripts.
- Access to the relevant input data files and reference resources required by each analysis.

Please consult the individual scripts for package dependencies, input file requirements, and configuration settings. Package versions and software environments are documented where necessary to facilitate reproducibility.

## Running the analyses

1. Clone or download this repository.
2. Obtain the required input data and metadata from the sources described in the manuscript and its Supplementary Materials.
3. Review the input and output paths specified in each script and adapt them to your local environment.
4. Install the required software and packages.
5. Run the scripts in the appropriate order, taking into account dependencies between analysis steps.

Shell scripts (`.sh` and `.bash`) are provided for selected analyses and can be used to launch the corresponding R scripts. Check their contents before execution to confirm the expected working directory, paths, and computational requirements.

## Contributors

- **Iris Sammarco** — Repository maintainer; development of the data analysis scripts. https://github.com/iris-sammarco/
- **Zofia Maria Szlachtowska** — Contribution of additional scripts for data analysis.

## Citation

If you use these scripts, please cite the associated manuscript:

*[Full manuscript citation to be added]*
