# ============================================================
# 00_setup.R
# Project setup for reproducible analyses
# ============================================================

# -----------------------------
# 1. Clear environment
# -----------------------------
rm(list = ls())

# -----------------------------
# 2. Load required packages
# -----------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(openxlsx)
library(DescTools)
library(boot)
library(coin)
library(lavaan)

# -----------------------------
# 3. Reproducibility seed
# -----------------------------
set.seed(123)

# -----------------------------
# 4. Define project paths
# -----------------------------
path_input  <- "input"
path_output <- "output"
path_r      <- "R"

# Optional output subfolders
path_output_tables <- file.path(path_output, "tables")
path_output_rds    <- file.path(path_output, "rds")

# -----------------------------
# 5. Create output folders if needed
# -----------------------------
dir.create(path_output, showWarnings = FALSE)
dir.create(path_output_tables, showWarnings = FALSE)
dir.create(path_output_rds, showWarnings = FALSE)

# -----------------------------
# 6. Input file path
# -----------------------------
path_dataset <- file.path(path_input, "dataset.xlsx")

# -----------------------------
# 7. Check that input file exists
# -----------------------------
if (!file.exists(path_dataset)) {
  stop("Input file not found at 'input/dataset.xlsx'.")
}