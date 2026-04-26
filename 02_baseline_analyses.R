# ============================================================
# 02_baseline_analyses.R
# Baseline comparisons (Table 1)
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
# 3. Helpers
# -----------------------------
fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), NA_character_, format(round(x, digits), nsmall = digits))
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001) return("<0.0001")
  sprintf("%.4f", p)
}

med_iqr_string <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  med <- median(x)
  q1  <- quantile(x, 0.25, names = FALSE)
  q3  <- quantile(x, 0.75, names = FALSE)
  paste0(fmt_num(med, digits), " [", fmt_num(q1, digits), "-", fmt_num(q3, digits), "]")
}

# -----------------------------
# 4. Continuous variables summary
# -----------------------------
summarize_cont_4groups <- function(data, var, digits = 1) {
  out <- data %>%
    group_by(grupo4) %>%
    summarise(
      resumen = med_iqr_string(.data[[var]], digits = digits),
      .groups = "drop"
    )
  
  vals <- setNames(out$resumen, out$grupo4)
  
  tibble(
    variable = var,
    `NP / LDA−` = vals["NP / LDA−"],
    `NP / LDA+` = vals["NP / LDA+"],
    `PE / LDA−` = vals["PE / LDA−"],
    `PE / LDA+` = vals["PE / LDA+"]
  )
}

# -----------------------------
# 5. Kruskal-Wallis + post hoc
# -----------------------------
kw_and_conover_4groups <- function(data, var, alpha_posthoc = 0.05) {
  dd <- data %>%
    select(grupo4, all_of(var)) %>%
    filter(!is.na(grupo4), !is.na(.data[[var]]))
  
  if (nrow(dd) == 0) {
    return(list(
      p_kw = NA_real_,
      kw = NULL,
      posthoc = NULL
    ))
  }
  
  kw <- kruskal.test(as.formula(paste(var, "~ grupo4")), data = dd)
  p_kw <- kw$p.value
  
  posthoc <- NULL
  if (!is.na(p_kw) && p_kw < alpha_posthoc) {
    posthoc <- DescTools::ConoverTest(
      x = dd[[var]],
      g = dd$grupo4,
      method = "holm"
    )
  }
  
  list(
    p_kw = p_kw,
    kw = kw,
    posthoc = posthoc
  )
}

# -----------------------------
# 6. Categorical variables
# -----------------------------
summarize_binary_4groups <- function(data, var, positive_value, digits = 1) {
  out <- data %>%
    group_by(grupo4) %>%
    summarise(
      pct_positive = 100 * mean(.data[[var]] == positive_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  vals <- setNames(fmt_num(out$pct_positive, digits), out$grupo4)
  
  tibble(
    variable = var,
    `NP / LDA−` = vals["NP / LDA−"],
    `NP / LDA+` = vals["NP / LDA+"],
    `PE / LDA−` = vals["PE / LDA−"],
    `PE / LDA+` = vals["PE / LDA+"]
  )
}

fisher_4groups <- function(data, var) {
  dd <- data %>%
    select(grupo4, all_of(var)) %>%
    filter(!is.na(grupo4), !is.na(.data[[var]]))
  
  tab <- table(dd[[var]], dd$grupo4)
  fisher.test(tab)$p.value
}

# ============================================================
# ANALYSES
# ============================================================

# -----------------------------
# 7. Continuous variables
# -----------------------------
cont_vars <- c("EDAD", "BMI_NUM", "GANANCIA_PESO", "PAS", "PAD", "EG")

cat("\n================ CONTINUOUS VARIABLES =================\n")

for (v in cont_vars) {
  cat("\nVariable:", v, "\n")
  print(summarize_cont_4groups(bd_clean, v), n = Inf)
  
  res <- kw_and_conover_4groups(bd_clean, v)
  cat("Kruskal-Wallis p-value:", fmt_p(res$p_kw), "\n")
  
  if (v == "EG") {
    if (!is.null(res$posthoc)) {
      cat("\nConover-Iman post hoc for EG (Holm-adjusted p-values):\n")
      print(res$posthoc)
    } else {
      cat("No post hoc performed for EG (global test not significant).\n")
    }
  }
}

# -----------------------------
# 8. LDA initiation (treated only)
# -----------------------------
bd_lda <- bd_clean %>%
  filter(TX == 1) %>%
  mutate(
    grupo_lda = ifelse(EMB == 0, "NP / LDA+", "PE / LDA+")
  )

cat("\n================ LDA INITIATION =================\n")

print(
  bd_lda %>%
    group_by(grupo_lda) %>%
    summarise(
      resumen = med_iqr_string(INICIO_SEM_TX_LDA),
      .groups = "drop"
    ),
  n = Inf
)

kw_inicio <- kruskal.test(INICIO_SEM_TX_LDA ~ grupo_lda, data = bd_lda)
cat("Kruskal-Wallis p-value:", fmt_p(kw_inicio$p.value), "\n")

# -----------------------------
# 9. OBESIDAD
# -----------------------------
cat("\n================ OBESITY =================\n")

tabla_ob <- summarize_binary_4groups(bd_clean, "OBESIDAD", "sí")
print(tabla_ob, n = Inf)

p_ob <- fisher_4groups(bd_clean, "OBESIDAD")
cat("Fisher p-value:", fmt_p(p_ob), "\n")

# -----------------------------
# 10. SEXO FETAL
# -----------------------------
cat("\n================ FETAL SEX =================\n")

tabla_fem <- summarize_binary_4groups(bd_clean, "SEXO_RN", "femenino")
tabla_masc <- summarize_binary_4groups(bd_clean, "SEXO_RN", "masculino")

print(tabla_fem, n = Inf)
print(tabla_masc, n = Inf)

p_sexo <- fisher_4groups(bd_clean, "SEXO_RN")
cat("Fisher p-value:", fmt_p(p_sexo), "\n")

# -----------------------------
# 11. Birthweight by sex
# -----------------------------
analyze_bw <- function(data, sex_value) {
  dd <- data %>%
    filter(SEXO_RN == sex_value)
  
  cat("\nBirthweight -", sex_value, "\n")
  print(summarize_cont_4groups(dd, "PESO_RN"), n = Inf)
  
  res <- kw_and_conover_4groups(dd, "PESO_RN")
  cat("Kruskal-Wallis p-value:", fmt_p(res$p_kw), "\n")
  
  if (!is.null(res$posthoc)) {
    cat("\nConover-Iman post hoc for birthweight in", sex_value, "(Holm-adjusted p-values):\n")
    print(res$posthoc)
  } else {
    cat("No post hoc performed for birthweight in", sex_value, "(global test not significant).\n")
  }
}

cat("\n================ BIRTHWEIGHT =================\n")

analyze_bw(bd_clean, "femenino")
analyze_bw(bd_clean, "masculino")

# ============================================================
# SAVE KEY RESULTS
# ============================================================

# -----------------------------
# 1. P-values summary
# -----------------------------
pvals_table1 <- tibble(
  variable = c(
    "EDAD",
    "BMI_NUM",
    "GANANCIA_PESO",
    "PAS",
    "PAD",
    "EG"
  ),
  p_value = c(
    kw_and_conover_4groups(bd_clean, "EDAD")$p_kw,
    kw_and_conover_4groups(bd_clean, "BMI_NUM")$p_kw,
    kw_and_conover_4groups(bd_clean, "GANANCIA_PESO")$p_kw,
    kw_and_conover_4groups(bd_clean, "PAS")$p_kw,
    kw_and_conover_4groups(bd_clean, "PAD")$p_kw,
    kw_and_conover_4groups(bd_clean, "EG")$p_kw
  )
)

# -----------------------------
# 2. Post hoc EG and BW
# -----------------------------
res_eg <- kw_and_conover_4groups(bd_clean, "EG")

posthoc_eg <- if (!is.null(res_eg$posthoc)) {
  as.data.frame(res_eg$posthoc[[1]]) %>%
    tibble::rownames_to_column(var = "comparison") %>%
    dplyr::rename(
      mean_rank_diff = `mean rank diff`,
      p_value = pval
    )
} else {
  NULL
}

res_bw_f <- kw_and_conover_4groups(
  bd_clean %>% filter(SEXO_RN == "femenino"),
  "PESO_RN"
)

posthoc_bw_f <- if (!is.null(res_bw_f$posthoc)) {
  as.data.frame(res_bw_f$posthoc[[1]]) %>%
    tibble::rownames_to_column(var = "comparison") %>%
    dplyr::rename(
      mean_rank_diff = `mean rank diff`,
      p_value = pval
    )
} else {
  NULL
}

res_bw_m <- kw_and_conover_4groups(
  bd_clean %>% filter(SEXO_RN == "masculino"),
  "PESO_RN"
)

posthoc_bw_m <- if (!is.null(res_bw_m$posthoc)) {
  as.data.frame(res_bw_m$posthoc[[1]]) %>%
    tibble::rownames_to_column(var = "comparison") %>%
    dplyr::rename(
      mean_rank_diff = `mean rank diff`,
      p_value = pval
    )
} else {
  NULL
}
# -----------------------------
# 3. Save to files
# -----------------------------
write.csv(pvals_table1,
          file = file.path(path_output_tables, "table1_pvalues.csv"),
          row.names = FALSE)

if (!is.null(posthoc_eg)) {
  write.csv(posthoc_eg,
            file = file.path(path_output_tables, "posthoc_EG.csv"),
            row.names = FALSE)
}

if (!is.null(posthoc_bw_f)) {
  write.csv(posthoc_bw_f,
            file = file.path(path_output_tables, "posthoc_BW_female.csv"),
            row.names = FALSE)
}

if (!is.null(posthoc_bw_m)) {
  write.csv(posthoc_bw_m,
            file = file.path(path_output_tables, "posthoc_BW_male.csv"),
            row.names = FALSE)
}

# ============================================================
# SUPPLEMENTARY BASELINE TABLE BY 8 MATERNAL PROFILES
# Descriptive only (no outcomes)
# ============================================================

# -----------------------------
# 1. Helpers for 8-group summaries
# -----------------------------
summarize_cont_8groups <- function(data, var, digits = 1) {
  out <- data %>%
    group_by(profile8) %>%
    summarise(
      resumen = med_iqr_string(.data[[var]], digits = digits),
      .groups = "drop"
    )
  
  vals <- setNames(out$resumen, out$profile8)
  
  tibble(
    variable = var,
    `NP / LDA− / OB−` = vals["NP / LDA− / OB−"],
    `NP / LDA− / OB+` = vals["NP / LDA− / OB+"],
    `NP / LDA+ / OB−` = vals["NP / LDA+ / OB−"],
    `NP / LDA+ / OB+` = vals["NP / LDA+ / OB+"],
    `PE / LDA− / OB−` = vals["PE / LDA− / OB−"],
    `PE / LDA− / OB+` = vals["PE / LDA− / OB+"],
    `PE / LDA+ / OB−` = vals["PE / LDA+ / OB−"],
    `PE / LDA+ / OB+` = vals["PE / LDA+ / OB+"]
  )
}

# -----------------------------
# 2. Variables to include
# Baseline only (no outcomes)
# -----------------------------
cont_vars_8profiles <- c("EDAD", "BMI_NUM", "GANANCIA_PESO", "PAS", "PAD", "INICIO_SEM_TX_LDA")

# -----------------------------
# 3. Build supplementary table
# -----------------------------
tabla_supp_8profiles <- bind_rows(
  lapply(cont_vars_8profiles, function(v) {
    summarize_cont_8groups(
      bd_clean,
      v,
      digits = ifelse(v == "INICIO_SEM_TX_LDA", 0, 1)
    )
  })
)

# -----------------------------
# 4. Optional relabeling for readability
# -----------------------------
tabla_supp_8profiles <- tabla_supp_8profiles %>%
  mutate(
    variable = case_when(
      variable == "EDAD" ~ "Maternal age",
      variable == "BMI_NUM" ~ "Pregestational BMI",
      variable == "GANANCIA_PESO" ~ "Maternal weight gain",
      variable == "PAS" ~ "Systolic blood pressure",
      variable == "PAD" ~ "Diastolic blood pressure",
      variable == "INICIO_SEM_TX_LDA" ~ "Initiation of LDA exposure (gestational week)",
      TRUE ~ variable
    )
  )

# -----------------------------
# 5. Print
# -----------------------------
cat("\n============================================\n")
cat("SUPPLEMENTARY BASELINE TABLE BY 8 MATERNAL PROFILES\n")
cat("============================================\n")
print(tabla_supp_8profiles, n = Inf)

# -----------------------------
# 6. Save
# -----------------------------
write.csv(
  tabla_supp_8profiles,
  file = file.path(path_output_tables, "supplementary_baseline_8profiles.csv"),
  row.names = FALSE
)