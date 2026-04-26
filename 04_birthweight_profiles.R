# ============================================================
# 04_birthweight_profiles.R
# Birthweight comparisons across 8 maternal profiles
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
profile_levels <- c("P000", "P001", "P010", "P011", "P100", "P101", "P110", "P111")

bd_female <- bd_female %>%
  mutate(
    perfil = factor(paste0("P", EMB, TX, OB), levels = profile_levels)
  )

bd_male <- bd_male %>%
  mutate(
    perfil = factor(paste0("P", EMB, TX, OB), levels = profile_levels)
  )

# -----------------------------
# 4. Helper: bootstrap mean
# -----------------------------
mean_boot_fun <- function(data, indices) {
  mean(data[indices], na.rm = TRUE)
}

# -----------------------------
# 5. Helper: safe percentile CI from boot
# -----------------------------
safe_boot_ci <- function(boot_object) {
  ci <- tryCatch(
    boot::boot.ci(boot_object, type = "perc"),
    error = function(e) NULL
  )
  
  if (is.null(ci) || is.null(ci$percent)) {
    return(c(NA_real_, NA_real_))
  }
  
  c(ci$percent[4], ci$percent[5])
}

# -----------------------------
# 6. Descriptive summary by profile
# -----------------------------
summarize_profiles_boot <- function(df, outcome_var = "PESO_RN", R = 10000) {
  
  perfiles_presentes <- levels(df$perfil)
  
  res_list <- lapply(perfiles_presentes, function(p) {
    x <- df[[outcome_var]][df$perfil == p]
    x <- x[!is.na(x)]
    
    n <- length(x)
    
    if (n == 0) {
      return(
        tibble(
          perfil = p,
          n = 0,
          media = NA_real_,
          sd = NA_real_,
          sem = NA_real_,
          ic95_low = NA_real_,
          ic95_high = NA_real_
        )
      )
    }
    
    media <- mean(x)
    sd_x  <- if (n > 1) sd(x) else NA_real_
    sem_x <- if (n > 1) sd_x / sqrt(n) else NA_real_
    
    boot_out <- boot::boot(
      data = x,
      statistic = mean_boot_fun,
      R = R
    )
    
    ci_vals <- safe_boot_ci(boot_out)
    
    tibble(
      perfil = p,
      n = n,
      media = round(media, 2),
      sd = round(sd_x, 2),
      sem = round(sem_x, 2),
      ic95_low = round(ci_vals[1], 2),
      ic95_high = round(ci_vals[2], 2)
    )
  })
  
  bind_rows(res_list)
}

# -----------------------------
# 7. Bootstrap difference in means
# -----------------------------
diff_mean_boot_fun <- function(data, indices) {
  d <- data[indices, , drop = FALSE]
  
  g1 <- d$PESO_RN[d$perfil == levels(d$perfil)[1]]
  g2 <- d$PESO_RN[d$perfil == levels(d$perfil)[2]]
  
  mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE)
}

# -----------------------------
# 8. Pairwise comparisons
# -----------------------------
pairwise_profile_comparisons <- function(df, outcome_var = "PESO_RN", R = 10000) {
  
  niveles <- levels(df$perfil)
  pares <- combn(niveles, 2, simplify = FALSE)
  
  out_list <- lapply(pares, function(par) {
    p1 <- par[1]
    p2 <- par[2]
    
    datos_par <- df %>%
      filter(perfil %in% c(p1, p2), !is.na(.data[[outcome_var]])) %>%
      mutate(perfil = droplevels(perfil))
    
    n1 <- sum(datos_par$perfil == p1)
    n2 <- sum(datos_par$perfil == p2)
    
    if (nrow(datos_par) == 0 || n1 == 0 || n2 == 0) {
      return(NULL)
    }
    
    test_par <- coin::oneway_test(
      as.formula(paste(outcome_var, "~ perfil")),
      data = datos_par,
      distribution = coin::approximate(nresample = 10000)
    )
    
    x1 <- datos_par[[outcome_var]][datos_par$perfil == p1]
    x2 <- datos_par[[outcome_var]][datos_par$perfil == p2]
    
    diff_media <- mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE)
    
    boot_data <- datos_par %>%
      select(perfil, all_of(outcome_var)) %>%
      rename(PESO_RN = all_of(outcome_var))
    
    boot_out <- boot::boot(
      data = boot_data,
      statistic = diff_mean_boot_fun,
      R = R
    )
    
    ci_vals <- safe_boot_ci(boot_out)
    
    tibble(
      perfil1 = p1,
      perfil2 = p2,
      n1 = n1,
      n2 = n2,
      diff_media = diff_media,
      ic95_low_diff = ci_vals[1],
      ic95_high_diff = ci_vals[2],
      p_perm = as.numeric(coin::pvalue(test_par))
    )
  })
  
  res <- bind_rows(out_list)
  
  if (nrow(res) == 0) {
    return(res)
  }
  
  res %>%
    mutate(
      p_bh = p.adjust(p_perm, method = "BH")
    ) %>%
    mutate(
      across(c(diff_media, ic95_low_diff, ic95_high_diff), ~ round(., 2)),
      across(c(p_perm, p_bh), ~ signif(., 4))
    ) %>%
    arrange(p_perm)
}

# -----------------------------
# 9. Global permutation ANOVA
# -----------------------------
run_global_permutation_anova <- function(df, outcome_var = "PESO_RN") {
  
  dd <- df %>%
    filter(!is.na(perfil), !is.na(.data[[outcome_var]]))
  
  test <- coin::oneway_test(
    as.formula(paste(outcome_var, "~ perfil")),
    data = dd,
    distribution = coin::approximate(nresample = 10000)
  )
  
  tibble(
    outcome = outcome_var,
    p_perm_global = as.numeric(coin::pvalue(test))
  )
}

# ============================================================
# 10. FEMALE FETUSES
# ============================================================

summary_female <- summarize_profiles_boot(
  df = bd_female,
  outcome_var = "PESO_RN",
  R = 10000
)

anova_female <- run_global_permutation_anova(
  df = bd_female,
  outcome_var = "PESO_RN"
)

pairs_female <- pairwise_profile_comparisons(
  df = bd_female,
  outcome_var = "PESO_RN",
  R = 10000
)

cat("\n============================================\n")
cat("BIRTHWEIGHT BY PROFILE - FEMALE FETUSES\n")
cat("============================================\n")
print(summary_female, n = Inf)

cat("\nGlobal permutation ANOVA p-value (female):\n")
print(anova_female, n = Inf)

cat("\nPairwise profile comparisons (female):\n")
print(pairs_female, n = Inf)

# ============================================================
# 11. MALE FETUSES
# ============================================================

summary_male <- summarize_profiles_boot(
  df = bd_male,
  outcome_var = "PESO_RN",
  R = 10000
)

anova_male <- run_global_permutation_anova(
  df = bd_male,
  outcome_var = "PESO_RN"
)

pairs_male <- pairwise_profile_comparisons(
  df = bd_male,
  outcome_var = "PESO_RN",
  R = 10000
)

cat("\n============================================\n")
cat("BIRTHWEIGHT BY PROFILE - MALE FETUSES\n")
cat("============================================\n")
print(summary_male, n = Inf)

cat("\nGlobal permutation ANOVA p-value (male):\n")
print(anova_male, n = Inf)

cat("\nPairwise profile comparisons (male):\n")
print(pairs_male, n = Inf)

# ============================================================
# 12. Save outputs
# ============================================================

write.csv(
  summary_female,
  file = file.path(path_output_tables, "birthweight_profile_summary_female.csv"),
  row.names = FALSE
)

write.csv(
  summary_male,
  file = file.path(path_output_tables, "birthweight_profile_summary_male.csv"),
  row.names = FALSE
)

write.csv(
  anova_female,
  file = file.path(path_output_tables, "birthweight_global_permutation_anova_female.csv"),
  row.names = FALSE
)

write.csv(
  anova_male,
  file = file.path(path_output_tables, "birthweight_global_permutation_anova_male.csv"),
  row.names = FALSE
)

write.csv(
  pairs_female,
  file = file.path(path_output_tables, "birthweight_pairwise_profiles_female.csv"),
  row.names = FALSE
)

write.csv(
  pairs_male,
  file = file.path(path_output_tables, "birthweight_pairwise_profiles_male.csv"),
  row.names = FALSE
)

saveRDS(
  summary_female,
  file = file.path(path_output_rds, "birthweight_profile_summary_female.rds")
)

saveRDS(
  summary_male,
  file = file.path(path_output_rds, "birthweight_profile_summary_male.rds")
)

saveRDS(
  anova_female,
  file = file.path(path_output_rds, "birthweight_global_permutation_anova_female.rds")
)

saveRDS(
  anova_male,
  file = file.path(path_output_rds, "birthweight_global_permutation_anova_male.rds")
)

saveRDS(
  pairs_female,
  file = file.path(path_output_rds, "birthweight_pairwise_profiles_female.rds")
)

saveRDS(
  pairs_male,
  file = file.path(path_output_rds, "birthweight_pairwise_profiles_male.rds")
)
