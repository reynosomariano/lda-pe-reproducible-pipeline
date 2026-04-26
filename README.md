
# Reproducible R pipeline

This repository contains the R scripts used to reproduce the statistical analyses of the study: NAME.

The pipeline is designed to run on a curated input dataset with the same variable structure as the original study dataset. Original clinical data are not publicly shared, but a data template and variable dictionary are provided so that users can adapt the scripts to their own data.

## Repository structure

- `00_setup.R`: loads packages, defines paths, and initializes the project environment.
- `01_data_preparation.R`: reads the curated input dataset and creates analytic variables, interaction terms, and profile groupings.
- `02_baseline_analyses.R`: reproduces baseline descriptive and inferential analyses.
- `03_discrete_outcomes.R`: analyzes discrete perinatal outcomes across maternal profiles.
- `04_birthweight_profiles.R`: analyzes birthweight differences across the 8 maternal profiles.
- `05_sem_female.R`: fits the moderated mediation SEM in female fetuses.
- `06_sem_male.R`: fits the moderated mediation SEM in male fetuses.
- `07_multigroup_sem.R`: fits the multi-group SEM to assess heterogeneity by fetal sex.
- `08_sensitivity_small_cells.R`: performs sensitivity analyses for small profile cells.
- `09_sensitivity_path_b.R`: performs sensitivity analyses focused on path b.

Folders: - `input/`: input dataset provided by the user. - `output/tables/`: exported result tables. - `output/rds/`: intermediate R objects and fitted model outputs. - `output/supplementary/`: supplementary result files. - `data_template/`: dataset template and data dictionary. - `R/`: optional helper functions.

## Requirements

This project was developed in R. Running the scripts is easier from RStudio, although this is not strictly required.

Main R packages used: - `readxl` - `dplyr` - `tidyr` - `tibble` - `purrr` - `openxlsx` - `DescTools` - `boot` - `coin` - `lavaan`

## Input data

The repository expects a curated Excel file placed at:

`input/dataset.xlsx`

The required worksheet name is:

`Datos Clínicos`

The scripts assume that the dataset already contains the curated study variables needed for analysis. The repository does not reconstruct the original data cleaning process performed outside R.

At minimum, the input dataset must include the variables required for: - pregnancy status, - low-dose aspirin exposure, - maternal obesity, - fetal sex, - gestational age, - birthweight, - preterm birth, - birthweight \< p10, - and all variables listed in the accompanying data dictionary.

Detailed variable definitions are provided in:

- `data_template/data_dictionary.csv`
- `data_template/template_dataset.csv`

## Execution order

Run the scripts in the following order:

1.  `00_setup.R`
2.  `01_data_preparation.R`
3.  `02_baseline_analyses.R`
4.  `03_discrete_outcomes.R`
5.  `04_birthweight_profiles.R`
6.  `05_sem_female.R`
7.  `06_sem_male.R`
8.  `07_multigroup_sem.R`
9.  `08_sensitivity_small_cells.R`
10. `09_sensitivity_path_b.R`

Each script assumes that the previous steps have already been run successfully.

### How to run the scripts

Open the project in RStudio, and run each script one by one using the "Source" button or by typing:

\`\`\`r source("script_name.R")

## Outputs

The scripts generate: - `.csv` tables for descriptive and inferential results, - `.xlsx` workbooks for SEM and sensitivity analyses, - `.rds` objects containing processed datasets, model fits, and intermediate outputs.

All outputs are written to the `output/` directory.

## Reproducibility notes

- The repository is intended to reproduce the analytical workflow, not the full original data curation process.
- The original clinical dataset is not publicly distributed.
- Users must provide their own dataset following the same variable structure described in `data_template/`.
- Figures produced outside R are not included in this repository.
- Random procedures were run with fixed seeds for reproducibility whenever applicable.

## Data availability

The original clinical dataset cannot be publicly shared due to privacy and confidentiality restrictions.

To facilitate reuse of the code, this repository includes: - a template dataset structure, - and a variable dictionary describing the expected input variables.
