# ============================================================
# 07_multigroup_sem.R
# Multi-group SEM to assess heterogeneity by fetal sex
# ============================================================

# -----------------------------
# 1. Load setup
# -----------------------------
source("00_setup.R")

# -----------------------------
# 2. Load prepared dataset
# -----------------------------
bd_clean <- readRDS(file.path(path_output_rds, "bd_clean.rds"))

# -----------------------------
# 3. Prepare data for multi-group SEM
# -----------------------------
bd_mg <- bd_clean %>%
  filter(SEXO_RN %in% c("femenino", "masculino")) %>%
  mutate(
    SEXO_RN = factor(SEXO_RN, levels = c("femenino", "masculino"))
  ) %>%
  select(
    SEXO_RN,
    EMB, TX, OB,
    EMB_TX, EMB_OB, TX_OB, EMB_TX_OB,
    EG, PESO_RN
  ) %>%
  filter(complete.cases(.))

cat("\n============================================\n")
cat("SAMPLE SIZE BY FETAL SEX\n")
cat("============================================\n")
print(table(bd_mg$SEXO_RN))

# -----------------------------
# 4. Helper: build model syntax
# equal_block = "none", "b", "a", "c", "global"
# -----------------------------
build_mg_model <- function(equal_block = c("none", "b", "a", "c", "global")) {
  equal_block <- match.arg(equal_block)
  
  lab <- function(name, block) {
    if (equal_block %in% c(block, "global")) {
      paste0("c(", name, ",", name, ")")
    } else {
      paste0("c(", name, "_fem,", name, "_masc)")
    }
  }
  
  paste0(
    '
# =========================
# MEDIATOR EQUATION
# =========================
EG ~ ', lab("a1", "a"), '*EMB
   + ', lab("a2", "a"), '*TX
   + ', lab("a3", "a"), '*OB
   + ', lab("a12", "a"), '*EMB_TX
   + ', lab("a13", "a"), '*EMB_OB
   + ', lab("a23", "a"), '*TX_OB
   + ', lab("a123", "a"), '*EMB_TX_OB

# =========================
# OUTCOME EQUATION
# =========================
PESO_RN ~ ', lab("c1", "c"), '*EMB
        + ', lab("c2", "c"), '*TX
        + ', lab("c3", "c"), '*OB
        + ', lab("c12", "c"), '*EMB_TX
        + ', lab("c13", "c"), '*EMB_OB
        + ', lab("c23", "c"), '*TX_OB
        + ', lab("c123", "c"), '*EMB_TX_OB
        + ', lab("b", "b"), '*EG

# =========================
# INTERCEPTS AND RESIDUAL VARIANCES
# (free across groups)
# =========================
EG ~ c(int_EG_fem, int_EG_masc)*1
PESO_RN ~ c(int_BW_fem, int_BW_masc)*1

EG ~~ c(var_EG_fem, var_EG_masc)*EG
PESO_RN ~~ c(var_BW_fem, var_BW_masc)*PESO_RN
'
  )
}

# -----------------------------
# 5. Fit helper
# -----------------------------
fit_mg_model <- function(model_syntax, data) {
  lavaan::sem(
    model = model_syntax,
    data = data,
    group = "SEXO_RN",
    estimator = "MLR",
    meanstructure = TRUE,
    fixed.x = TRUE
  )
}

# -----------------------------
# 6. Fit free and constrained models
# -----------------------------
model_free   <- build_mg_model("none")
model_b_eq   <- build_mg_model("b")
model_a_eq   <- build_mg_model("a")
model_c_eq   <- build_mg_model("c")
model_glob_eq <- build_mg_model("global")

fit_free    <- fit_mg_model(model_free, bd_mg)
fit_b_eq    <- fit_mg_model(model_b_eq, bd_mg)
fit_a_eq    <- fit_mg_model(model_a_eq, bd_mg)
fit_c_eq    <- fit_mg_model(model_c_eq, bd_mg)
fit_glob_eq <- fit_mg_model(model_glob_eq, bd_mg)

# -----------------------------
# 7. Model comparisons
# -----------------------------
compare_models <- function(fit_restricted, fit_free, comparison_name) {
  lrt <- lavaan::lavTestLRT(fit_restricted, fit_free)
  
  tibble(
    comparison = comparison_name,
    chisq_diff = if ("Chisq diff" %in% names(lrt)) lrt[2, "Chisq diff"] else NA_real_,
    df_diff    = if ("Df diff" %in% names(lrt)) lrt[2, "Df diff"] else NA_real_,
    p_value    = if ("Pr(>Chisq)" %in% names(lrt)) lrt[2, "Pr(>Chisq)"] else NA_real_,
    aic_restricted = AIC(fit_restricted),
    bic_restricted = BIC(fit_restricted),
    aic_free       = AIC(fit_free),
    bic_free       = BIC(fit_free)
  )
}

table_model_comparisons <- bind_rows(
  compare_models(fit_b_eq, fit_free, "b_equal_vs_free"),
  compare_models(fit_a_eq, fit_free, "a_block_equal_vs_free"),
  compare_models(fit_c_eq, fit_free, "c_block_equal_vs_free"),
  compare_models(fit_glob_eq, fit_free, "global_equal_vs_free")
)

cat("\n============================================\n")
cat("MULTI-GROUP SEM MODEL COMPARISONS\n")
cat("============================================\n")
print(table_model_comparisons, n = Inf)

# -----------------------------
# 8. Parameter estimates from free model
# -----------------------------
pe_free <- lavaan::parameterEstimates(
  fit_free,
  standardized = TRUE,
  ci = TRUE
) %>%
  tibble::as_tibble()

cat("\n============================================\n")
cat("FREE MULTI-GROUP MODEL PARAMETERS\n")
cat("============================================\n")
print(pe_free, n = Inf)

# -----------------------------
# 9. Path a by sex
# -----------------------------
table_a_mg <- pe_free %>%
  filter(lhs == "EG", op == "~") %>%
  select(group, lhs, op, rhs, est, se, z, pvalue, ci.lower, ci.upper, std.all) %>%
  mutate(
    sexo = case_when(
      group == 1 ~ "femenino",
      group == 2 ~ "masculino",
      TRUE ~ as.character(group)
    )
  ) %>%
  select(sexo, everything(), -group)

cat("\n============================================\n")
cat("PATH a IN FREE MULTI-GROUP MODEL\n")
cat("============================================\n")
print(table_a_mg, n = Inf)

# -----------------------------
# 10. Path b and c by sex
# -----------------------------
table_bc_mg <- pe_free %>%
  filter(lhs == "PESO_RN", op == "~") %>%
  select(group, lhs, op, rhs, est, se, z, pvalue, ci.lower, ci.upper, std.all) %>%
  mutate(
    sexo = case_when(
      group == 1 ~ "femenino",
      group == 2 ~ "masculino",
      TRUE ~ as.character(group)
    )
  ) %>%
  select(sexo, everything(), -group)

cat("\n============================================\n")
cat("PATH b AND c IN FREE MULTI-GROUP MODEL\n")
cat("============================================\n")
print(table_bc_mg, n = Inf)

# -----------------------------
# 11. Fit measures for all models
# -----------------------------
extract_fit <- function(fit, model_name) {
  lavaan::fitMeasures(
    fit,
    c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
  ) %>%
    tibble::enframe(name = "indice", value = "valor") %>%
    mutate(modelo = model_name) %>%
    select(modelo, everything())
}

table_fit_mg <- bind_rows(
  extract_fit(fit_free, "free"),
  extract_fit(fit_b_eq, "b_equal"),
  extract_fit(fit_a_eq, "a_equal"),
  extract_fit(fit_c_eq, "c_equal"),
  extract_fit(fit_glob_eq, "global_equal")
)

cat("\n============================================\n")
cat("FIT INDICES FOR MULTI-GROUP MODELS\n")
cat("============================================\n")
print(table_fit_mg, n = Inf)

# -----------------------------
# 12. Save outputs
# -----------------------------
openxlsx::write.xlsx(
  list(
    Comparisons = table_model_comparisons,
    Free_model_parameters = pe_free,
    Path_a = table_a_mg,
    Path_bc = table_bc_mg,
    Fit = table_fit_mg
  ),
  file = file.path(path_output_tables, "MultiGroup_SEM_FetalSex.xlsx"),
  rowNames = FALSE
)

saveRDS(fit_free, file = file.path(path_output_rds, "fit_mg_free.rds"))
saveRDS(fit_b_eq, file = file.path(path_output_rds, "fit_mg_b_equal.rds"))
saveRDS(fit_a_eq, file = file.path(path_output_rds, "fit_mg_a_equal.rds"))
saveRDS(fit_c_eq, file = file.path(path_output_rds, "fit_mg_c_equal.rds"))
saveRDS(fit_glob_eq, file = file.path(path_output_rds, "fit_mg_global_equal.rds"))
saveRDS(table_model_comparisons, file = file.path(path_output_rds, "mg_model_comparisons.rds"))