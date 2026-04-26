# ============================================================
# 03_discrete_outcomes.R
# Discrete perinatal outcomes by maternal profile and fetal sex
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
# 3. Create short profile code
#    P000 = EMB0 TX0 OB0, etc.
# -----------------------------
bd_female <- bd_female %>%
  mutate(
    perfil = factor(paste0("P", EMB, TX, OB))
  )

bd_male <- bd_male %>%
  mutate(
    perfil = factor(paste0("P", EMB, TX, OB))
  )

# -----------------------------
# 4. Core function:
#    compare one profile vs reference
# -----------------------------
compare_profile_vs_ref <- function(df,
                                   outcome_var,
                                   profile_level,
                                   ref_level = "P000",
                                   B = 10000,
                                   seed = 123) {
  set.seed(seed)
  
  dat <- df %>%
    filter(!is.na(perfil), !is.na(.data[[outcome_var]])) %>%
    mutate(
      perfil = as.character(perfil),
      .event = .data[[outcome_var]] == 1
    )
  
  d_ref <- dat %>% filter(perfil == ref_level)
  d_g   <- dat %>% filter(perfil == profile_level)
  
  n_ref <- nrow(d_ref)
  n_g   <- nrow(d_g)
  
  if (n_ref == 0 || n_g == 0) {
    return(NULL)
  }
  
  s_ref <- sum(d_ref$.event)
  s_g   <- sum(d_g$.event)
  
  p_ref <- s_ref / n_ref
  p_g   <- s_g / n_g
  
  # RR with Haldane–Anscombe correction
  rr_hat <- ((s_g + 0.5) / (n_g + 1)) / ((s_ref + 0.5) / (n_ref + 1))
  
  # Risk difference / attributable risk
  rd_hat <- p_g - p_ref
  ar_hat <- rd_hat
  
  # Fisher exact test
  mat <- matrix(
    c(s_g, n_g - s_g,
      s_ref, n_ref - s_ref),
    nrow = 2,
    byrow = TRUE
  )
  
  p_fisher <- tryCatch(
    fisher.test(mat)$p.value,
    error = function(e) NA_real_
  )
  
  # Stratified bootstrap within groups
  rr_boot <- numeric(B)
  rd_boot <- numeric(B)
  
  for (b in seq_len(B)) {
    idx_ref <- sample.int(n_ref, n_ref, replace = TRUE)
    idx_g   <- sample.int(n_g, n_g, replace = TRUE)
    
    s_ref_b <- sum(d_ref$.event[idx_ref])
    s_g_b   <- sum(d_g$.event[idx_g])
    
    p_ref_b <- s_ref_b / n_ref
    p_g_b   <- s_g_b / n_g
    
    rr_boot[b] <- ((s_g_b + 0.5) / (n_g + 1)) / ((s_ref_b + 0.5) / (n_ref + 1))
    rd_boot[b] <- p_g_b - p_ref_b
  }
  
  rr_ci <- quantile(rr_boot, c(0.025, 0.975), na.rm = TRUE)
  rd_ci <- quantile(rd_boot, c(0.025, 0.975), na.rm = TRUE)
  
  tibble(
    outcome     = outcome_var,
    perfil      = profile_level,
    n_perfil    = n_g,
    eventos     = s_g,
    prop_evento = p_g,
    ref         = ref_level,
    n_ref       = n_ref,
    eventos_ref = s_ref,
    prop_ref    = p_ref,
    RR          = rr_hat,
    RR_low      = as.numeric(rr_ci[1]),
    RR_high     = as.numeric(rr_ci[2]),
    RD          = rd_hat,
    RD_low      = as.numeric(rd_ci[1]),
    RD_high     = as.numeric(rd_ci[2]),
    AR          = ar_hat,
    AR_low      = as.numeric(rd_ci[1]),
    AR_high     = as.numeric(rd_ci[2]),
    p_raw       = p_fisher
  )
}

# -----------------------------
# 5. Run one outcome across profiles
#    and apply BH correction
# -----------------------------
analyze_outcome_profiles <- function(df,
                                     outcome_var,
                                     ref_level = "P000",
                                     B = 10000,
                                     seed = 123) {
  
  perfiles <- df %>%
    filter(!is.na(perfil)) %>%
    distinct(perfil) %>%
    mutate(perfil = as.character(perfil)) %>%
    pull(perfil) %>%
    unique()
  
  perfiles <- setdiff(perfiles, ref_level)
  
  res <- purrr::map_dfr(
    perfiles,
    ~ compare_profile_vs_ref(
      df = df,
      outcome_var = outcome_var,
      profile_level = .x,
      ref_level = ref_level,
      B = B,
      seed = seed
    )
  )
  
  if (nrow(res) == 0) {
    return(res)
  }
  
  res %>%
    mutate(
      p_BH = p.adjust(p_raw, method = "BH")
    )
}

# -----------------------------
# 6. Format final tables
# -----------------------------
format_results <- function(res, outcome_label) {
  res %>%
    mutate(
      Outcome = outcome_label,
      
      prop_evento_pp = round(100 * prop_evento, 2),
      prop_ref_pp    = round(100 * prop_ref, 2),
      
      RR      = round(RR, 2),
      RR_low  = round(RR_low, 2),
      RR_high = round(RR_high, 2),
      
      RD      = round(100 * RD, 2),
      RD_low  = round(100 * RD_low, 2),
      RD_high = round(100 * RD_high, 2),
      
      AR      = round(100 * AR, 2),
      AR_low  = round(100 * AR_low, 2),
      AR_high = round(100 * AR_high, 2),
      
      p_raw = round(p_raw, 4),
      p_BH  = round(p_BH, 4)
    ) %>%
    transmute(
      Outcome,
      Perfil = perfil,
      n = n_perfil,
      Eventos = eventos,
      `Proporción evento (%)` = prop_evento_pp,
      Ref = ref,
      n_ref,
      `Proporción ref (%)` = prop_ref_pp,
      `RR [IC95%]` = paste0(RR, " [", RR_low, "–", RR_high, "]"),
      `RD (pp) [IC95%]` = paste0(RD, " [", RD_low, "–", RD_high, "]"),
      `AR (pp) [IC95%]` = paste0(AR, " [", AR_low, "–", AR_high, "]"),
      `p (Fisher)` = sprintf("%.4f", p_raw),
      `p (BH)` = sprintf("%.4f", p_BH)
    ) %>%
    arrange(Outcome, Perfil)
}

# ============================================================
# 7. FEMALE FETUSES
# ============================================================

res_ptb_female <- analyze_outcome_profiles(
  df = bd_female,
  outcome_var = "PTB",
  ref_level = "P000",
  B = 10000,
  seed = 123
)

res_bwp10_female <- analyze_outcome_profiles(
  df = bd_female,
  outcome_var = "BW_P10",
  ref_level = "P000",
  B = 10000,
  seed = 123
)

table_discrete_female <- bind_rows(
  format_results(res_ptb_female, "Preterm birth"),
  format_results(res_bwp10_female, "Birthweight < p10")
)

cat("\n============================================\n")
cat("DISCRETE OUTCOMES - FEMALE FETUSES\n")
cat("============================================\n")
print(table_discrete_female, n = Inf)

# ============================================================
# 8. MALE FETUSES
# ============================================================

res_ptb_male <- analyze_outcome_profiles(
  df = bd_male,
  outcome_var = "PTB",
  ref_level = "P000",
  B = 10000,
  seed = 123
)

res_bwp10_male <- analyze_outcome_profiles(
  df = bd_male,
  outcome_var = "BW_P10",
  ref_level = "P000",
  B = 10000,
  seed = 123
)

table_discrete_male <- bind_rows(
  format_results(res_ptb_male, "Preterm birth"),
  format_results(res_bwp10_male, "Birthweight < p10")
)

cat("\n============================================\n")
cat("DISCRETE OUTCOMES - MALE FETUSES\n")
cat("============================================\n")
print(table_discrete_male, n = Inf)

# ============================================================
# 9. Save outputs
# ============================================================

write.csv(
  table_discrete_female,
  file = file.path(path_output_tables, "discrete_outcomes_female.csv"),
  row.names = FALSE
)

write.csv(
  table_discrete_male,
  file = file.path(path_output_tables, "discrete_outcomes_male.csv"),
  row.names = FALSE
)

saveRDS(
  res_ptb_female,
  file = file.path(path_output_rds, "res_ptb_female.rds")
)

saveRDS(
  res_bwp10_female,
  file = file.path(path_output_rds, "res_bwp10_female.rds")
)

saveRDS(
  res_ptb_male,
  file = file.path(path_output_rds, "res_ptb_male.rds")
)

saveRDS(
  res_bwp10_male,
  file = file.path(path_output_rds, "res_bwp10_male.rds")
)