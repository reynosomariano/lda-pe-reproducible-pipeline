# ============================================================
# 01_data_preparation.R
# Read curated dataset and create analytic variables
# ============================================================

# -----------------------------
# 1. Load project setup
# -----------------------------
source("00_setup.R")

# -----------------------------
# 2. Read Excel file
# -----------------------------
bd <- read_excel(path_dataset, sheet = "Datos Clínicos")

# -----------------------------
# 3. Check required variables
# -----------------------------
required_vars <- c(
  "n",
  "EMBARAZO",
  "TX_LDA",
  "INICIO_SEM_TX_LDA",
  "EDAD",
  "PESO_I",
  "PESO_F",
  "GANANCIA_PESO",
  "ALTURA",
  "BMI_NUM",
  "BMI_CAT",
  "OBESIDAD",
  "PAS",
  "PAD",
  "PROT_CAT",
  "PESO_RN",
  "SEXO_RN",
  "EG",
  "PARTO",
  "RCIU"
)

missing_vars <- setdiff(required_vars, names(bd))

if (length(missing_vars) > 0) {
  stop(
    paste(
      "The following required variables are missing from the dataset:",
      paste(missing_vars, collapse = ", ")
    )
  )
}

# -----------------------------
# 4. Create analytic variables
# -----------------------------
bd_clean <- bd %>%
  mutate(
    EMB = case_when(
      EMBARAZO == "NP" ~ 0,
      EMBARAZO == "PE" ~ 1,
      TRUE ~ NA_real_
    ),
    
    TX = case_when(
      TX_LDA == "no tratada" ~ 0,
      TX_LDA == "tratada"    ~ 1,
      TRUE ~ NA_real_
    ),
    
    OB = case_when(
      OBESIDAD == "no" ~ 0,
      OBESIDAD == "sí" ~ 1,
      TRUE ~ NA_real_
    ),
    
    FETAL_SEX = factor(
      SEXO_RN,
      levels = c("femenino", "masculino"),
      ordered = TRUE
    ),
    
    PTB = case_when(
      PARTO == "termino"    ~ 0,
      PARTO == "pretermino" ~ 1,
      TRUE ~ NA_real_
    ),
    
    BW_P10 = case_when(
      RCIU == "No" ~ 0,
      RCIU == "Sí" ~ 1,
      TRUE ~ NA_real_
    ),
    
    EMB_TX = EMB * TX,
    EMB_OB = EMB * OB,
    TX_OB = TX * OB,
    EMB_TX_OB = EMB * TX * OB
  )
# -----------------------------
# 5. Create 4-level pregnancy/treatment groups
# -----------------------------
bd_clean <- bd_clean %>%
  mutate(
    grupo4 = case_when(
      EMB == 0 & TX == 0 ~ "NP / LDA−",
      EMB == 0 & TX == 1 ~ "NP / LDA+",
      EMB == 1 & TX == 0 ~ "PE / LDA−",
      EMB == 1 & TX == 1 ~ "PE / LDA+",
      TRUE ~ NA_character_
    ),
    grupo4 = factor(
      grupo4,
      levels = c("NP / LDA−", "NP / LDA+", "PE / LDA−", "PE / LDA+")
    )
  )

# -----------------------------
# 6. Create 8-level maternal profile
# -----------------------------
bd_clean <- bd_clean %>%
  mutate(
    profile8 = case_when(
      EMB == 0 & TX == 0 & OB == 0 ~ "NP / LDA− / OB−",
      EMB == 0 & TX == 0 & OB == 1 ~ "NP / LDA− / OB+",
      EMB == 0 & TX == 1 & OB == 0 ~ "NP / LDA+ / OB−",
      EMB == 0 & TX == 1 & OB == 1 ~ "NP / LDA+ / OB+",
      EMB == 1 & TX == 0 & OB == 0 ~ "PE / LDA− / OB−",
      EMB == 1 & TX == 0 & OB == 1 ~ "PE / LDA− / OB+",
      EMB == 1 & TX == 1 & OB == 0 ~ "PE / LDA+ / OB−",
      EMB == 1 & TX == 1 & OB == 1 ~ "PE / LDA+ / OB+",
      TRUE ~ NA_character_
    ),
    profile8 = factor(
      profile8,
      levels = c(
        "NP / LDA− / OB−",
        "NP / LDA− / OB+",
        "NP / LDA+ / OB−",
        "NP / LDA+ / OB+",
        "PE / LDA− / OB−",
        "PE / LDA− / OB+",
        "PE / LDA+ / OB−",
        "PE / LDA+ / OB+"
      )
    )
  )

# -----------------------------
# 7. Split datasets by fetal sex
# -----------------------------
bd_female <- bd_clean %>%
  filter(FETAL_SEX == "femenino")

bd_male <- bd_clean %>%
  filter(FETAL_SEX == "masculino")

# -----------------------------
# 8. Quick checks
# -----------------------------
cat("\n============================================\n")
cat("SAMPLE SIZE\n")
cat("============================================\n")
cat("Total sample:", nrow(bd_clean), "\n")
cat("Female fetuses:", nrow(bd_female), "\n")
cat("Male fetuses:", nrow(bd_male), "\n")

cat("\n============================================\n")
cat("GROUP COUNTS: grupo4\n")
cat("============================================\n")
print(table(bd_clean$grupo4, useNA = "ifany"))

cat("\n============================================\n")
cat("GROUP COUNTS: profile8\n")
cat("============================================\n")
print(table(bd_clean$profile8, useNA = "ifany"))

# -----------------------------
# 9. Save cleaned datasets
# -----------------------------
saveRDS(bd_clean,  file = file.path(path_output_rds, "bd_clean.rds"))
saveRDS(bd_female, file = file.path(path_output_rds, "bd_female.rds"))
saveRDS(bd_male,   file = file.path(path_output_rds, "bd_male.rds"))