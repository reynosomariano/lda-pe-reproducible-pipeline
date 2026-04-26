# ============================================================
# 05_sem_female.R
# Moderated mediation SEM in female fetuses
# ============================================================

# -----------------------------
# 1. Load setup
# -----------------------------
source("00_setup.R")

# -----------------------------
# 2. Load female dataset
# -----------------------------
bd_female <- readRDS(file.path(path_output_rds, "bd_female.rds"))

# -----------------------------
# 3. Define SEM model
# -----------------------------
modelo_female <- '
# Mediator model
EG ~ a1*EMB + a2*TX + a3*OB +
     a12*EMB_TX + a13*EMB_OB + a23*TX_OB +
     a123*EMB_TX_OB

# Outcome model
PESO_RN ~ c1*EMB + c2*TX + c3*OB +
          c12*EMB_TX + c13*EMB_OB + c23*TX_OB +
          c123*EMB_TX_OB +
          b*EG

# EG contrasts anchored to reference profile
EG_PE_TX0_OB0 := a1
EG_PE_TX0_OB1 := a1 + a3 + a13
EG_PE_TX1_OB0 := a1 + a2 + a12
EG_PE_TX1_OB1 := a1 + a2 + a3 + a12 + a13 + a23 + a123
EG_NP_TX0_OB1 := a3
EG_NP_TX1_OB0 := a2
EG_NP_TX1_OB1 := a2 + a3 + a23

# Direct effects on birthweight
DIR_PE_TX0_OB0 := c1
DIR_PE_TX0_OB1 := c1 + c3 + c13
DIR_PE_TX1_OB0 := c1 + c2 + c12
DIR_PE_TX1_OB1 := c1 + c2 + c3 + c12 + c13 + c23 + c123
DIR_NP_TX0_OB1 := c3
DIR_NP_TX1_OB0 := c2
DIR_NP_TX1_OB1 := c2 + c3 + c23

# Indirect effects
IND_PE_TX0_OB0 := EG_PE_TX0_OB0 * b
IND_PE_TX0_OB1 := EG_PE_TX0_OB1 * b
IND_PE_TX1_OB0 := EG_PE_TX1_OB0 * b
IND_PE_TX1_OB1 := EG_PE_TX1_OB1 * b
IND_NP_TX0_OB1 := EG_NP_TX0_OB1 * b
IND_NP_TX1_OB0 := EG_NP_TX1_OB0 * b
IND_NP_TX1_OB1 := EG_NP_TX1_OB1 * b

# Total effects
TOT_PE_TX0_OB0 := DIR_PE_TX0_OB0 + IND_PE_TX0_OB0
TOT_PE_TX0_OB1 := DIR_PE_TX0_OB1 + IND_PE_TX0_OB1
TOT_PE_TX1_OB0 := DIR_PE_TX1_OB0 + IND_PE_TX1_OB0
TOT_PE_TX1_OB1 := DIR_PE_TX1_OB1 + IND_PE_TX1_OB1
TOT_NP_TX0_OB1 := DIR_NP_TX0_OB1 + IND_NP_TX0_OB1
TOT_NP_TX1_OB0 := DIR_NP_TX1_OB0 + IND_NP_TX1_OB0
TOT_NP_TX1_OB1 := DIR_NP_TX1_OB1 + IND_NP_TX1_OB1

# Proportion mediated (%)
PM_PE_TX0_OB0 := (IND_PE_TX0_OB0 / TOT_PE_TX0_OB0) * 100
PM_PE_TX0_OB1 := (IND_PE_TX0_OB1 / TOT_PE_TX0_OB1) * 100
PM_PE_TX1_OB0 := (IND_PE_TX1_OB0 / TOT_PE_TX1_OB0) * 100
PM_PE_TX1_OB1 := (IND_PE_TX1_OB1 / TOT_PE_TX1_OB1) * 100
PM_NP_TX0_OB1 := (IND_NP_TX0_OB1 / TOT_NP_TX0_OB1) * 100
PM_NP_TX1_OB0 := (IND_NP_TX1_OB0 / TOT_NP_TX1_OB0) * 100
PM_NP_TX1_OB1 := (IND_NP_TX1_OB1 / TOT_NP_TX1_OB1) * 100
'

# -----------------------------
# 4. Fit model with bootstrap
# -----------------------------
set.seed(123)

fit_female <- lavaan::sem(
  modelo_female,
  data = bd_female,
  meanstructure = TRUE,
  se = "bootstrap",
  bootstrap = 10000
)

# -----------------------------
# 5. Full parameter table
# -----------------------------
pe_female <- lavaan::parameterEstimates(
  fit_female,
  standardized = TRUE,
  ci = TRUE,
  boot.ci.type = "perc"
) %>%
  as.data.frame() %>%
  tibble::as_tibble()

# -----------------------------
# 6. Path a table
# -----------------------------
table_path_a_female <- pe_female %>%
  dplyr::filter(lhs == "EG", op == "~") %>%
  dplyr::select(lhs, op, rhs, label, est, se, z, pvalue, ci.lower, ci.upper, std.all) %>%
  dplyr::arrange(match(label, c("a1", "a2", "a3", "a12", "a13", "a23", "a123")), rhs)

# -----------------------------
# 7. Path b and c' table
# -----------------------------
table_path_bc_female <- pe_female %>%
  dplyr::filter(lhs == "PESO_RN", op == "~") %>%
  dplyr::select(lhs, op, rhs, label, est, se, z, pvalue, ci.lower, ci.upper, std.all) %>%
  dplyr::arrange(match(label, c("c1", "c2", "c3", "c12", "c13", "c23", "c123", "b")), rhs)

# -----------------------------
# 8. Defined effects table
# -----------------------------
table_effects_female <- pe_female %>%
  dplyr::filter(op == ":=") %>%
  dplyr::transmute(
    tipo = dplyr::case_when(
      grepl("^IND_", label) ~ "Indirecto",
      grepl("^DIR_", label) ~ "Directo",
      grepl("^TOT_", label) ~ "Total",
      grepl("^PM_", label)  ~ "Prop. mediada",
      TRUE ~ "Otros"
    ),
    efecto = label,
    definicion = rhs,
    est,
    se,
    z,
    pvalue,
    ci.lower,
    ci.upper,
    std.all
  ) %>%
  dplyr::arrange(
    factor(tipo, levels = c("Indirecto", "Directo", "Total", "Prop. mediada", "Otros")),
    efecto
  )

# -----------------------------
# 9. Technical tables
# -----------------------------
technical_female <- list(
  nobs = tryCatch(lavaan::lavInspect(fit_female, "nobs"), error = function(e) NA),
  estimator = tryCatch(fit_female@Options$estimator, error = function(e) NA),
  se = tryCatch(fit_female@Options$se, error = function(e) NA),
  bootstrap = tryCatch(fit_female@Options$bootstrap, error = function(e) NA)
) %>%
  tibble::as_tibble()

table_var_resid_female <- pe_female %>%
  dplyr::filter(op == "~~", lhs == rhs) %>%
  dplyr::transmute(
    parametro = paste0("Var(", lhs, ")"),
    est, se, z, pvalue, ci.lower, ci.upper, std.all
  ) %>%
  dplyr::arrange(parametro)

table_r2_female <- tryCatch(
  lavaan::lavInspect(fit_female, "r2") %>%
    tibble::enframe(name = "variable", value = "R2") %>%
    tibble::as_tibble(),
  error = function(e) tibble::tibble(variable = NA_character_, R2 = NA_real_)
)

table_fit_female <- tryCatch(
  lavaan::fitMeasures(
    fit_female,
    c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
  ) %>%
    tibble::enframe(name = "indice", value = "valor") %>%
    tibble::as_tibble(),
  error = function(e) tibble::tibble(indice = NA_character_, valor = NA_real_)
)

# -----------------------------
# 10. Print key results
# -----------------------------
cat("\n===== FEMALE: Path a (EG ~ ...) =====\n")
print(table_path_a_female, n = Inf)

cat("\n===== FEMALE: Path b and c' (PESO_RN ~ ... + EG) =====\n")
print(table_path_bc_female, n = Inf)

cat("\n===== FEMALE: Defined effects (:=) =====\n")
print(table_effects_female, n = Inf)

cat("\n===== FEMALE: Technical details =====\n")
print(technical_female)

cat("\n===== FEMALE: Residual variances =====\n")
print(table_var_resid_female, n = Inf)

cat("\n===== FEMALE: R-squared =====\n")
print(table_r2_female, n = Inf)

cat("\n===== FEMALE: Fit measures =====\n")
print(table_fit_female, n = Inf)

# -----------------------------
# 11. Save outputs
# -----------------------------
openxlsx::write.xlsx(
  list(
    Path_a = table_path_a_female,
    Path_bc = table_path_bc_female,
    Effects = table_effects_female,
    Residual_variances = table_var_resid_female,
    R2 = table_r2_female,
    Fit = table_fit_female
  ),
  file = file.path(path_output_tables, "SEM_results_female.xlsx"),
  rowNames = FALSE
)

saveRDS(fit_female, file = file.path(path_output_rds, "fit_sem_female.rds"))
saveRDS(pe_female, file = file.path(path_output_rds, "sem_parameters_female.rds"))
saveRDS(table_effects_female, file = file.path(path_output_rds, "sem_effects_female.rds"))