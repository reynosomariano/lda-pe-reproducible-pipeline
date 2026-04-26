# ============================================================
# 08_sensitivity_small_cells.R
# Sensitivity analysis for small profile cells
# - Leave-one-out (LOO) for profiles with n <= 6
# - Sign stability by bootstrap
# ============================================================

# -----------------------------
# 1. Load setup
# -----------------------------
source("00_setup.R")

# -----------------------------
# 2. Load prepared datasets
# -----------------------------
bd_female <- readRDS(file.path(path_output_rds, "bd_female.rds"))
bd_male   <- readRDS(file.path(path_output_rds, "bd_male.rds"))

# -----------------------------
# 3. Recreate profile code
# -----------------------------
profile_levels <- c("P000", "P001", "P010", "P011", "P100", "P101", "P110", "P111")

bd_female <- bd_female %>%
  mutate(perfil = factor(paste0("P", EMB, TX, OB), levels = profile_levels))

bd_male <- bd_male %>%
  mutate(perfil = factor(paste0("P", EMB, TX, OB), levels = profile_levels))

# -----------------------------
# 4. Recreate SEM model syntax
# -----------------------------
modelo_female <- '
EG ~ a1*EMB + a2*TX + a3*OB +
     a12*EMB_TX + a13*EMB_OB + a23*TX_OB +
     a123*EMB_TX_OB

PESO_RN ~ c1*EMB + c2*TX + c3*OB +
          c12*EMB_TX + c13*EMB_OB + c23*TX_OB +
          c123*EMB_TX_OB +
          b*EG

EG_PE_TX0_OB0 := a1
EG_PE_TX0_OB1 := a1 + a3 + a13
EG_PE_TX1_OB0 := a1 + a2 + a12
EG_PE_TX1_OB1 := a1 + a2 + a3 + a12 + a13 + a23 + a123
EG_NP_TX0_OB1 := a3
EG_NP_TX1_OB0 := a2
EG_NP_TX1_OB1 := a2 + a3 + a23

DIR_PE_TX0_OB0 := c1
DIR_PE_TX0_OB1 := c1 + c3 + c13
DIR_PE_TX1_OB0 := c1 + c2 + c12
DIR_PE_TX1_OB1 := c1 + c2 + c3 + c12 + c13 + c23 + c123
DIR_NP_TX0_OB1 := c3
DIR_NP_TX1_OB0 := c2
DIR_NP_TX1_OB1 := c2 + c3 + c23

IND_PE_TX0_OB0 := EG_PE_TX0_OB0 * b
IND_PE_TX0_OB1 := EG_PE_TX0_OB1 * b
IND_PE_TX1_OB0 := EG_PE_TX1_OB0 * b
IND_PE_TX1_OB1 := EG_PE_TX1_OB1 * b
IND_NP_TX0_OB1 := EG_NP_TX0_OB1 * b
IND_NP_TX1_OB0 := EG_NP_TX1_OB0 * b
IND_NP_TX1_OB1 := EG_NP_TX1_OB1 * b

TOT_PE_TX0_OB0 := DIR_PE_TX0_OB0 + IND_PE_TX0_OB0
TOT_PE_TX0_OB1 := DIR_PE_TX0_OB1 + IND_PE_TX0_OB1
TOT_PE_TX1_OB0 := DIR_PE_TX1_OB0 + IND_PE_TX1_OB0
TOT_PE_TX1_OB1 := DIR_PE_TX1_OB1 + IND_PE_TX1_OB1
TOT_NP_TX0_OB1 := DIR_NP_TX0_OB1 + IND_NP_TX0_OB1
TOT_NP_TX1_OB0 := DIR_NP_TX1_OB0 + IND_NP_TX1_OB0
TOT_NP_TX1_OB1 := DIR_NP_TX1_OB1 + IND_NP_TX1_OB1

PM_PE_TX0_OB0 := (IND_PE_TX0_OB0 / TOT_PE_TX0_OB0) * 100
PM_PE_TX0_OB1 := (IND_PE_TX0_OB1 / TOT_PE_TX0_OB1) * 100
PM_PE_TX1_OB0 := (IND_PE_TX1_OB0 / TOT_PE_TX1_OB0) * 100
PM_PE_TX1_OB1 := (IND_PE_TX1_OB1 / TOT_PE_TX1_OB1) * 100
PM_NP_TX0_OB1 := (IND_NP_TX0_OB1 / TOT_NP_TX0_OB1) * 100
PM_NP_TX1_OB0 := (IND_NP_TX1_OB0 / TOT_NP_TX1_OB0) * 100
PM_NP_TX1_OB1 := (IND_NP_TX1_OB1 / TOT_NP_TX1_OB1) * 100
'

modelo_male <- '
EG ~ a1*EMB + a2*TX + a3*OB +
     a12*EMB_TX + a13*EMB_OB + a23*TX_OB +
     a123*EMB_TX_OB

PESO_RN ~ c1*EMB + c2*TX + c3*OB +
          c12*EMB_TX + c13*EMB_OB + c23*TX_OB +
          c123*EMB_TX_OB +
          b*EG

EG_PE_TX0_OB0 := a1
EG_PE_TX0_OB1 := a1 + a3 + a13
EG_PE_TX1_OB0 := a1 + a2 + a12
EG_PE_TX1_OB1 := a1 + a2 + a3 + a12 + a13 + a23 + a123
EG_NP_TX0_OB1 := a3
EG_NP_TX1_OB0 := a2
EG_NP_TX1_OB1 := a2 + a3 + a23

DIR_PE_TX0_OB0 := c1
DIR_PE_TX0_OB1 := c1 + c3 + c13
DIR_PE_TX1_OB0 := c1 + c2 + c12
DIR_PE_TX1_OB1 := c1 + c2 + c3 + c12 + c13 + c23 + c123
DIR_NP_TX0_OB1 := c3
DIR_NP_TX1_OB0 := c2
DIR_NP_TX1_OB1 := c2 + c3 + c23

IND_PE_TX0_OB0 := EG_PE_TX0_OB0 * b
IND_PE_TX0_OB1 := EG_PE_TX0_OB1 * b
IND_PE_TX1_OB0 := EG_PE_TX1_OB0 * b
IND_PE_TX1_OB1 := EG_PE_TX1_OB1 * b
IND_NP_TX0_OB1 := EG_NP_TX0_OB1 * b
IND_NP_TX1_OB0 := EG_NP_TX1_OB0 * b
IND_NP_TX1_OB1 := EG_NP_TX1_OB1 * b

TOT_PE_TX0_OB0 := DIR_PE_TX0_OB0 + IND_PE_TX0_OB0
TOT_PE_TX0_OB1 := DIR_PE_TX0_OB1 + IND_PE_TX0_OB1
TOT_PE_TX1_OB0 := DIR_PE_TX1_OB0 + IND_PE_TX1_OB0
TOT_PE_TX1_OB1 := DIR_PE_TX1_OB1 + IND_PE_TX1_OB1
TOT_NP_TX0_OB1 := DIR_NP_TX0_OB1 + IND_NP_TX0_OB1
TOT_NP_TX1_OB0 := DIR_NP_TX1_OB0 + IND_NP_TX1_OB0
TOT_NP_TX1_OB1 := DIR_NP_TX1_OB1 + IND_NP_TX1_OB1

PM_PE_TX0_OB0 := (IND_PE_TX0_OB0 / TOT_PE_TX0_OB0) * 100
PM_PE_TX0_OB1 := (IND_PE_TX0_OB1 / TOT_PE_TX0_OB1) * 100
PM_PE_TX1_OB0 := (IND_PE_TX1_OB0 / TOT_PE_TX1_OB0) * 100
PM_PE_TX1_OB1 := (IND_PE_TX1_OB1 / TOT_PE_TX1_OB1) * 100
PM_NP_TX0_OB1 := (IND_NP_TX0_OB1 / TOT_NP_TX0_OB1) * 100
PM_NP_TX1_OB0 := (IND_NP_TX1_OB0 / TOT_NP_TX1_OB0) * 100
PM_NP_TX1_OB1 := (IND_NP_TX1_OB1 / TOT_NP_TX1_OB1) * 100
'

# -----------------------------
# 5. User parameters
# -----------------------------
small_n_threshold <- 6
include_direct_effects <- FALSE
R_boot_sign <- 2000

set.seed(123)

# -----------------------------
# 6. Helpers
# -----------------------------
fmt2 <- function(x) round(as.numeric(x), 2)
fmt4 <- function(x) round(as.numeric(x), 4)

safe_sign <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x > 0 ~ "positive",
    x < 0 ~ "negative",
    TRUE ~ "zero"
  )
}

perfil_to_suffix <- function(perfil_code) {
  code <- gsub("^P", "", as.character(perfil_code))
  emb <- substr(code, 1, 1)
  tx  <- substr(code, 2, 2)
  ob  <- substr(code, 3, 3)
  
  emb_txt <- ifelse(emb == "1", "PE", "NP")
  paste0(emb_txt, "_TX", tx, "_OB", ob)
}

build_target_effect_labels <- function(profile_codes, include_direct = FALSE) {
  if (length(profile_codes) == 0) return(character(0))
  
  suffixes <- unique(vapply(profile_codes, perfil_to_suffix, character(1)))
  
  prefixes <- c("IND_", "TOT_")
  if (isTRUE(include_direct)) {
    prefixes <- c(prefixes, "DIR_")
  }
  
  as.vector(outer(prefixes, suffixes, paste0))
}

count_profiles <- function(df, sex_label) {
  df %>%
    filter(!is.na(perfil)) %>%
    count(perfil, name = "n") %>%
    arrange(n, perfil) %>%
    mutate(
      sexo = sex_label,
      suffix = vapply(perfil, perfil_to_suffix, character(1)),
      is_small = n <= small_n_threshold
    ) %>%
    select(sexo, perfil, suffix, n, is_small)
}

fit_sem_point <- function(model_syntax, data) {
  tryCatch(
    lavaan::sem(
      model = model_syntax,
      data = data,
      meanstructure = TRUE,
      se = "none"
    ),
    error = function(e) e
  )
}

extract_defined_effects <- function(fit_obj, labels_keep) {
  if (length(labels_keep) == 0) {
    return(tibble(effect = character(0), est = numeric(0)))
  }
  
  if (inherits(fit_obj, "error")) {
    return(tibble(effect = labels_keep, est = NA_real_))
  }
  
  pe <- lavaan::parameterEstimates(
    fit_obj,
    standardized = FALSE,
    ci = FALSE
  ) %>%
    as_tibble() %>%
    filter(op == ":=") %>%
    mutate(
      effect_name = case_when(
        !is.na(lhs) & lhs != "" ~ lhs,
        !is.na(label) & label != "" ~ label,
        TRUE ~ NA_character_
      )
    ) %>%
    select(effect = effect_name, est) %>%
    filter(!is.na(effect))
  
  tibble(effect = labels_keep) %>%
    left_join(pe, by = "effect")
}

# -----------------------------
# 7. Identify small profiles
# -----------------------------
small_profiles_female <- count_profiles(bd_female, "FEM")
small_profiles_male   <- count_profiles(bd_male, "MASC")

table_small_profiles <- bind_rows(small_profiles_female, small_profiles_male)

small_codes_female <- table_small_profiles %>%
  filter(sexo == "FEM", is_small) %>%
  pull(perfil) %>%
  as.character()

small_codes_male <- table_small_profiles %>%
  filter(sexo == "MASC", is_small) %>%
  pull(perfil) %>%
  as.character()

target_effects_female <- build_target_effect_labels(
  profile_codes = small_codes_female,
  include_direct = include_direct_effects
)

target_effects_male <- build_target_effect_labels(
  profile_codes = small_codes_male,
  include_direct = include_direct_effects
)

# -----------------------------
# 8. Leave-one-out
# -----------------------------
run_loo_small_cells <- function(df, model_syntax, sex_label, target_effects, small_profile_codes) {
  
  if (length(small_profile_codes) == 0 || length(target_effects) == 0) {
    return(list(
      full_effects = tibble(),
      loo_detail = tibble(),
      loo_summary = tibble()
    ))
  }
  
  dd <- df %>%
    mutate(.row_id_internal = row_number())
  
  fit_full <- fit_sem_point(model_syntax, dd %>% select(-.row_id_internal))
  full_effects <- extract_defined_effects(fit_full, target_effects) %>%
    rename(full_est = est) %>%
    mutate(
      full_sign = safe_sign(full_est),
      sexo = sex_label
    )
  
  loo_rows <- dd %>%
    filter(as.character(perfil) %in% small_profile_codes) %>%
    select(.row_id_internal, perfil)
  
  loo_detail <- purrr::map_dfr(seq_len(nrow(loo_rows)), function(i) {
    row_id_remove <- loo_rows$.row_id_internal[i]
    removed_profile <- as.character(loo_rows$perfil[i])
    
    dd_minus1 <- dd %>%
      filter(.row_id_internal != row_id_remove) %>%
      select(-.row_id_internal)
    
    fit_minus1 <- fit_sem_point(model_syntax, dd_minus1)
    est_minus1 <- extract_defined_effects(fit_minus1, target_effects) %>%
      rename(loo_est = est)
    
    full_effects %>%
      left_join(est_minus1, by = "effect") %>%
      mutate(
        sexo = sex_label,
        removed_profile = removed_profile,
        removed_row_id = row_id_remove,
        loo_sign = safe_sign(loo_est),
        sign_changed = if_else(full_sign != loo_sign, TRUE, FALSE, missing = NA),
        delta_est = loo_est - full_est
      )
  })
  
  loo_summary <- loo_detail %>%
    group_by(sexo, effect, removed_profile) %>%
    summarise(
      n_reestimates = n(),
      full_est = first(full_est),
      full_sign = first(full_sign),
      min_loo_est = ifelse(all(is.na(loo_est)), NA_real_, min(loo_est, na.rm = TRUE)),
      max_loo_est = ifelse(all(is.na(loo_est)), NA_real_, max(loo_est, na.rm = TRUE)),
      median_loo_est = ifelse(all(is.na(loo_est)), NA_real_, stats::median(loo_est, na.rm = TRUE)),
      pct_same_sign = ifelse(all(is.na(sign_changed)), NA_real_, mean(!sign_changed, na.rm = TRUE)),
      pct_sign_changed = ifelse(all(is.na(sign_changed)), NA_real_, mean(sign_changed, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      full_est = fmt2(full_est),
      min_loo_est = fmt2(min_loo_est),
      max_loo_est = fmt2(max_loo_est),
      median_loo_est = fmt2(median_loo_est),
      pct_same_sign = fmt4(pct_same_sign),
      pct_sign_changed = fmt4(pct_sign_changed)
    )
  
  list(
    full_effects = full_effects,
    loo_detail = loo_detail,
    loo_summary = loo_summary
  )
}

res_loo_female <- run_loo_small_cells(
  df = bd_female,
  model_syntax = modelo_female,
  sex_label = "FEM",
  target_effects = target_effects_female,
  small_profile_codes = small_codes_female
)

res_loo_male <- run_loo_small_cells(
  df = bd_male,
  model_syntax = modelo_male,
  sex_label = "MASC",
  target_effects = target_effects_male,
  small_profile_codes = small_codes_male
)

table_loo_summary <- bind_rows(
  res_loo_female$loo_summary,
  res_loo_male$loo_summary
)

table_loo_detail <- bind_rows(
  res_loo_female$loo_detail,
  res_loo_male$loo_detail
)

# -----------------------------
# 9. Sign stability by bootstrap
# -----------------------------
run_sign_stability <- function(df, model_syntax, sex_label, target_effects, R = 2000, seed = 123) {
  
  if (length(target_effects) == 0) {
    return(tibble())
  }
  
  fit0 <- fit_sem_point(model_syntax, df)
  
  if (inherits(fit0, "error")) {
    return(tibble(
      sexo = sex_label,
      effect = target_effects,
      full_est = NA_real_,
      full_sign = NA_character_,
      boot_median = NA_real_,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      prop_positive = NA_real_,
      prop_negative = NA_real_,
      prop_zero = NA_real_,
      prop_same_sign = NA_real_,
      R = R
    ))
  }
  
  full_effects <- extract_defined_effects(fit0, target_effects) %>%
    rename(full_est = est) %>%
    mutate(full_sign = safe_sign(full_est))
  
  fun_extract_defined <- function(fit) {
    pe <- lavaan::parameterEstimates(fit, standardized = FALSE, ci = FALSE) %>%
      as_tibble() %>%
      filter(op == ":=") %>%
      mutate(
        effect_name = case_when(
          !is.na(lhs) & lhs != "" ~ lhs,
          !is.na(label) & label != "" ~ label,
          TRUE ~ NA_character_
        )
      ) %>%
      select(effect = effect_name, est) %>%
      filter(!is.na(effect))
    
    vals <- pe$est[match(target_effects, pe$effect)]
    names(vals) <- target_effects
    vals
  }
  
  set.seed(seed)
  boot_out <- lavaan::bootstrapLavaan(
    object = fit0,
    R = R,
    FUN = fun_extract_defined,
    verbose = FALSE
  )
  
  boot_out <- as.matrix(boot_out)
  if (ncol(boot_out) == length(target_effects)) {
    colnames(boot_out) <- target_effects
  }
  
  sign_summary <- purrr::map_dfr(target_effects, function(lbl) {
    x <- boot_out[, lbl]
    
    tibble(
      effect = lbl,
      boot_median = stats::median(x, na.rm = TRUE),
      ci95_low = stats::quantile(x, probs = 0.025, na.rm = TRUE),
      ci95_high = stats::quantile(x, probs = 0.975, na.rm = TRUE),
      prop_positive = mean(x > 0, na.rm = TRUE),
      prop_negative = mean(x < 0, na.rm = TRUE),
      prop_zero = mean(x == 0, na.rm = TRUE)
    )
  })
  
  full_effects %>%
    left_join(sign_summary, by = "effect") %>%
    mutate(
      sexo = sex_label,
      prop_same_sign = case_when(
        full_sign == "positive" ~ prop_positive,
        full_sign == "negative" ~ prop_negative,
        full_sign == "zero" ~ prop_zero,
        TRUE ~ NA_real_
      ),
      R = R,
      full_est = fmt2(full_est),
      boot_median = fmt2(boot_median),
      ci95_low = fmt2(ci95_low),
      ci95_high = fmt2(ci95_high),
      prop_positive = fmt4(prop_positive),
      prop_negative = fmt4(prop_negative),
      prop_zero = fmt4(prop_zero),
      prop_same_sign = fmt4(prop_same_sign)
    ) %>%
    select(
      sexo, effect, full_est, full_sign,
      boot_median, ci95_low, ci95_high,
      prop_positive, prop_negative, prop_zero, prop_same_sign, R
    )
}

table_sign_female <- run_sign_stability(
  df = bd_female,
  model_syntax = modelo_female,
  sex_label = "FEM",
  target_effects = target_effects_female,
  R = R_boot_sign,
  seed = 123
)

table_sign_male <- run_sign_stability(
  df = bd_male,
  model_syntax = modelo_male,
  sex_label = "MASC",
  target_effects = target_effects_male,
  R = R_boot_sign,
  seed = 123
)

table_sign_summary <- bind_rows(table_sign_female, table_sign_male)

# -----------------------------
# 10. Print
# -----------------------------
cat("\n============================================\n")
cat("PROFILE COUNTS AND SMALL CELLS\n")
cat("============================================\n")
print(table_small_profiles, n = Inf)

cat("\n============================================\n")
cat("LEAVE-ONE-OUT SUMMARY\n")
cat("============================================\n")
print(table_loo_summary, n = Inf)

cat("\n============================================\n")
cat("SIGN STABILITY BY BOOTSTRAP\n")
cat("============================================\n")
print(table_sign_summary, n = Inf)

# -----------------------------
# 11. Save outputs
# -----------------------------
openxlsx::write.xlsx(
  list(
    Small_profiles = table_small_profiles,
    LOO_summary = table_loo_summary,
    LOO_detail = table_loo_detail,
    Sign_stability = table_sign_summary
  ),
  file = file.path(path_output_tables, "Sensitivity_small_cells_SEM.xlsx"),
  rowNames = FALSE
)

saveRDS(table_small_profiles, file = file.path(path_output_rds, "small_profiles_summary.rds"))
saveRDS(table_loo_summary, file = file.path(path_output_rds, "loo_summary_small_cells.rds"))
saveRDS(table_loo_detail, file = file.path(path_output_rds, "loo_detail_small_cells.rds"))
saveRDS(table_sign_summary, file = file.path(path_output_rds, "sign_stability_small_cells.rds"))