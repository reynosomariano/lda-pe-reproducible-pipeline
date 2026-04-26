# ============================================================
# 09_sensitivity_path_b.R
# Sensitivity analyses for path b
# - Specific path b block test
# - Stepwise SEM sensitivity in male fetuses
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
# 3. Formatting helpers
# -----------------------------
fmt2_num <- function(x) round(as.numeric(x), 2)
fmt4_num <- function(x) round(as.numeric(x), 4)
fmt4_chr <- function(x) sprintf("%.4f", as.numeric(x))

# -----------------------------
# 4. Prepare data for path b testing
# -----------------------------
prep_path_b_data <- function(df) {
  df %>%
    filter(
      !is.na(PESO_RN),
      !is.na(EG),
      !is.na(EMB),
      !is.na(TX),
      !is.na(OB)
    ) %>%
    mutate(
      EMB = as.numeric(EMB),
      TX  = as.numeric(TX),
      OB  = as.numeric(OB),
      EG_c = EG - mean(EG, na.rm = TRUE)
    )
}

# -----------------------------
# 5. Freedman-Lane block permutation test
# -----------------------------
perm_block_test_FL <- function(dd, form0, form1, B = 10000, seed = 123) {
  set.seed(seed)
  
  m0 <- lm(form0, data = dd)
  m1 <- lm(form1, data = dd)
  
  rss0 <- sum(resid(m0)^2)
  rss1 <- sum(resid(m1)^2)
  stat_obs <- rss0 - rss1
  
  yhat0 <- fitted(m0)
  e0    <- resid(m0)
  yname <- all.vars(form0)[1]
  
  stat_perm <- numeric(B)
  
  for (b in seq_len(B)) {
    y_star <- yhat0 + sample(e0, replace = FALSE)
    dd_star <- dd
    dd_star[[yname]] <- y_star
    
    m0b <- lm(form0, data = dd_star)
    m1b <- lm(form1, data = dd_star)
    
    stat_perm[b] <- sum(resid(m0b)^2) - sum(resid(m1b)^2)
  }
  
  p_perm <- (1 + sum(stat_perm >= stat_obs, na.rm = TRUE)) / (B + 1)
  
  tibble(
    B = B,
    RSS_reducido = rss0,
    RSS_completo = rss1,
    Delta_RSS = stat_obs,
    prop_mejora = ifelse(rss0 > 0, stat_obs / rss0, NA_real_),
    p_perm_bloque = p_perm
  )
}

# -----------------------------
# 6. Bootstrap of path b coefficients
# -----------------------------
boot_path_b_coeffs <- function(df, R = 10000, seed = 123) {
  set.seed(seed)
  
  dd <- prep_path_b_data(df)
  
  stat_fun <- function(data, idx) {
    d <- data[idx, ]
    
    fit <- lm(PESO_RN ~ EG_c * EMB * TX * OB, data = d)
    cf  <- coef(fit)
    
    get_cf <- function(name) {
      if (name %in% names(cf)) cf[[name]] else NA_real_
    }
    
    c(
      b0          = get_cf("EG_c"),
      b_EMB       = get_cf("EG_c:EMB"),
      b_TX        = get_cf("EG_c:TX"),
      b_OB        = get_cf("EG_c:OB"),
      b_EMB_TX    = get_cf("EG_c:EMB:TX"),
      b_EMB_OB    = get_cf("EG_c:EMB:OB"),
      b_TX_OB     = get_cf("EG_c:TX:OB"),
      b_EMB_TX_OB = get_cf("EG_c:EMB:TX:OB")
    )
  }
  
  boot_obj <- boot::boot(data = dd, statistic = stat_fun, R = R)
  
  if (is.null(colnames(boot_obj$t))) {
    colnames(boot_obj$t) <- c(
      "b0", "b_EMB", "b_TX", "b_OB",
      "b_EMB_TX", "b_EMB_OB", "b_TX_OB", "b_EMB_TX_OB"
    )
  }
  
  boot_obj
}

# -----------------------------
# 7. Summarize bootstrap coefficients
# -----------------------------
summarize_boot_b_coeffs <- function(boot_obj, sexo_label = NULL) {
  Tmat <- boot_obj$t
  
  if (is.null(colnames(Tmat))) {
    colnames(Tmat) <- c(
      "b0", "b_EMB", "b_TX", "b_OB",
      "b_EMB_TX", "b_EMB_OB", "b_TX_OB", "b_EMB_TX_OB"
    )
  }
  
  out <- tibble(
    termino = colnames(Tmat),
    estimador = apply(Tmat, 2, mean, na.rm = TRUE),
    ic95_low = apply(Tmat, 2, quantile, probs = 0.025, na.rm = TRUE),
    ic95_high = apply(Tmat, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      estimador = fmt2_num(estimador),
      ic95_low  = fmt2_num(ic95_low),
      ic95_high = fmt2_num(ic95_high)
    )
  
  if (!is.null(sexo_label)) {
    out <- out %>%
      mutate(sexo = sexo_label) %>%
      select(sexo, everything())
  }
  
  out
}

# -----------------------------
# 8. Extract profile-specific slopes
# -----------------------------
extract_b_slopes_from_boot <- function(boot_obj, sexo_label = NULL) {
  Tmat <- boot_obj$t
  
  if (is.null(colnames(Tmat))) {
    colnames(Tmat) <- c(
      "b0", "b_EMB", "b_TX", "b_OB",
      "b_EMB_TX", "b_EMB_OB", "b_TX_OB", "b_EMB_TX_OB"
    )
  }
  
  slopes_mat <- cbind(
    `NP / LDA- / OB-` = Tmat[, "b0"],
    `PE / LDA- / OB-` = Tmat[, "b0"] + Tmat[, "b_EMB"],
    `NP / LDA+ / OB-` = Tmat[, "b0"] + Tmat[, "b_TX"],
    `NP / LDA- / OB+` = Tmat[, "b0"] + Tmat[, "b_OB"],
    `PE / LDA+ / OB-` = Tmat[, "b0"] + Tmat[, "b_EMB"] + Tmat[, "b_TX"] + Tmat[, "b_EMB_TX"],
    `PE / LDA- / OB+` = Tmat[, "b0"] + Tmat[, "b_EMB"] + Tmat[, "b_OB"] + Tmat[, "b_EMB_OB"],
    `NP / LDA+ / OB+` = Tmat[, "b0"] + Tmat[, "b_TX"] + Tmat[, "b_OB"] + Tmat[, "b_TX_OB"],
    `PE / LDA+ / OB+` = Tmat[, "b0"] + Tmat[, "b_EMB"] + Tmat[, "b_TX"] + Tmat[, "b_OB"] +
      Tmat[, "b_EMB_TX"] + Tmat[, "b_EMB_OB"] + Tmat[, "b_TX_OB"] + Tmat[, "b_EMB_TX_OB"]
  )
  
  out <- tibble(
    perfil = colnames(slopes_mat),
    pendiente_b = apply(slopes_mat, 2, mean, na.rm = TRUE),
    ic95_low = apply(slopes_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    ic95_high = apply(slopes_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      pendiente_b = fmt2_num(pendiente_b),
      ic95_low = fmt2_num(ic95_low),
      ic95_high = fmt2_num(ic95_high)
    )
  
  if (!is.null(sexo_label)) {
    out <- out %>%
      mutate(sexo = sexo_label) %>%
      select(sexo, everything())
  }
  
  out
}

# -----------------------------
# 9. Wrapper by sex
# -----------------------------
analyze_path_b_reviewer <- function(df, sexo_label = "FEM", B_perm = 10000, R_boot = 10000, seed = 123) {
  
  dd <- prep_path_b_data(df)
  
  form_reducido <- PESO_RN ~ EMB * TX * OB + EG_c
  form_completo <- PESO_RN ~ EMB * TX * OB + EG_c * EMB * TX * OB
  
  test_bloque <- perm_block_test_FL(
    dd = dd,
    form0 = form_reducido,
    form1 = form_completo,
    B = B_perm,
    seed = seed
  ) %>%
    mutate(sexo = sexo_label) %>%
    select(sexo, everything())
  
  mod_reducido <- lm(form_reducido, data = dd)
  mod_completo <- lm(form_completo, data = dd)
  
  comparacion_lm <- tibble(
    sexo = sexo_label,
    modelo = c("reducido_b_comun", "completo_b_moderado"),
    n = nrow(dd),
    df_residual = c(df.residual(mod_reducido), df.residual(mod_completo)),
    sigma = c(summary(mod_reducido)$sigma, summary(mod_completo)$sigma),
    r2 = c(summary(mod_reducido)$r.squared, summary(mod_completo)$r.squared),
    r2_adj = c(summary(mod_reducido)$adj.r.squared, summary(mod_completo)$adj.r.squared),
    AIC = c(AIC(mod_reducido), AIC(mod_completo)),
    BIC = c(BIC(mod_reducido), BIC(mod_completo))
  ) %>%
    mutate(
      sigma = fmt2_num(sigma),
      r2 = fmt4_num(r2),
      r2_adj = fmt4_num(r2_adj),
      AIC = fmt2_num(AIC),
      BIC = fmt2_num(BIC)
    )
  
  boot_obj <- boot_path_b_coeffs(df, R = R_boot, seed = seed)
  
  tabla_coef_b <- summarize_boot_b_coeffs(boot_obj, sexo_label = sexo_label)
  tabla_slopes_b <- extract_b_slopes_from_boot(boot_obj, sexo_label = sexo_label)
  
  list(
    test_bloque = test_bloque,
    comparacion_lm = comparacion_lm,
    tabla_coef_b = tabla_coef_b,
    tabla_slopes_b = tabla_slopes_b,
    modelo_reducido = mod_reducido,
    modelo_completo = mod_completo,
    boot_obj = boot_obj
  )
}

# -----------------------------
# 10. Run specific path b test by sex
# -----------------------------
set.seed(123)

res_b_female <- analyze_path_b_reviewer(
  df = bd_female,
  sexo_label = "FEM",
  B_perm = 10000,
  R_boot = 10000,
  seed = 123
)

res_b_male <- analyze_path_b_reviewer(
  df = bd_male,
  sexo_label = "MASC",
  B_perm = 10000,
  R_boot = 10000,
  seed = 123
)

table_test_bloque_b <- bind_rows(
  res_b_female$test_bloque,
  res_b_male$test_bloque
) %>%
  mutate(
    RSS_reducido = fmt2_num(RSS_reducido),
    RSS_completo = fmt2_num(RSS_completo),
    Delta_RSS = fmt2_num(Delta_RSS),
    prop_mejora = fmt4_num(prop_mejora),
    p_perm_bloque = fmt4_chr(p_perm_bloque)
  )

table_comparacion_lm_b <- bind_rows(
  res_b_female$comparacion_lm,
  res_b_male$comparacion_lm
)

table_coeficientes_b <- bind_rows(
  res_b_female$tabla_coef_b,
  res_b_male$tabla_coef_b
)

table_pendientes_b <- bind_rows(
  res_b_female$tabla_slopes_b,
  res_b_male$tabla_slopes_b
)

cat("\n==============================\n")
cat("PATH b BLOCK TEST\n")
cat("==============================\n")
print(table_test_bloque_b, n = Inf)

cat("\n====================================\n")
cat("OLS COMPARISON: REDUCED VS COMPLETE\n")
cat("====================================\n")
print(table_comparacion_lm_b, n = Inf)

cat("\n============================================\n")
cat("BOOTSTRAP COEFFICIENTS FOR EG INTERACTIONS\n")
cat("============================================\n")
print(table_coeficientes_b, n = Inf)

cat("\n=========================================\n")
cat("PATH b SLOPES BY CLINICAL PROFILE\n")
cat("=========================================\n")
print(table_pendientes_b, n = Inf)

# -----------------------------
# 11. Male SEM sensitivity: stepwise release of path b
# -----------------------------
bd_male_sens <- bd_male %>%
  filter(
    !is.na(EG),
    !is.na(PESO_RN),
    !is.na(EMB),
    !is.na(TX),
    !is.na(OB)
  ) %>%
  mutate(
    EG_EMB       = EG * EMB,
    EG_TX        = EG * TX,
    EG_OB        = EG * OB,
    EG_EMB_TX    = EG * EMB * TX,
    EG_EMB_OB    = EG * EMB * OB,
    EG_TX_OB     = EG * TX * OB,
    EG_EMB_TX_OB = EG * EMB * TX * OB
  )

perfil_info <- tibble(
  nombre = c(
    "PE_TX0_OB0", "PE_TX0_OB1", "PE_TX1_OB0", "PE_TX1_OB1",
    "NP_TX0_OB1", "NP_TX1_OB0", "NP_TX1_OB1"
  ),
  EMB = c(1, 1, 1, 1, 0, 0, 0),
  TX  = c(0, 0, 1, 1, 0, 1, 1),
  OB  = c(0, 1, 0, 1, 1, 0, 1),
  EG_expr = c(
    "a1",
    "a1 + a3 + a13",
    "a1 + a2 + a12",
    "a1 + a2 + a3 + a12 + a13 + a23 + a123",
    "a3",
    "a2",
    "a2 + a3 + a23"
  ),
  DIR_expr = c(
    "c1",
    "c1 + c3 + c13",
    "c1 + c2 + c12",
    "c1 + c2 + c3 + c12 + c13 + c23 + c123",
    "c3",
    "c2",
    "c2 + c3 + c23"
  )
)

build_b_expr <- function(tipo, EMB, TX, OB) {
  expr <- "b0"
  
  if (tipo %in% c("emb", "main", "full")) {
    if (EMB == 1) expr <- paste(expr, "+ bE")
  }
  
  if (tipo %in% c("main", "full")) {
    if (TX == 1) expr <- paste(expr, "+ bT")
    if (OB == 1) expr <- paste(expr, "+ bO")
  }
  
  if (tipo == "full") {
    if (EMB == 1 && TX == 1) expr <- paste(expr, "+ bET")
    if (EMB == 1 && OB == 1) expr <- paste(expr, "+ bEO")
    if (TX == 1 && OB == 1) expr <- paste(expr, "+ bTO")
    if (EMB == 1 && TX == 1 && OB == 1) expr <- paste(expr, "+ bETO")
  }
  
  expr
}

build_defined_block <- function(tipo = c("common", "emb", "main", "full")) {
  tipo <- match.arg(tipo)
  
  out <- c()
  
  for (i in seq_len(nrow(perfil_info))) {
    nm  <- perfil_info$nombre[i]
    EMB <- perfil_info$EMB[i]
    TX  <- perfil_info$TX[i]
    OB  <- perfil_info$OB[i]
    
    eg_expr  <- perfil_info$EG_expr[i]
    dir_expr <- perfil_info$DIR_expr[i]
    b_expr   <- build_b_expr(tipo, EMB, TX, OB)
    
    out <- c(
      out,
      paste0("B_", nm, " := ", b_expr),
      paste0("IND_", nm, " := (", eg_expr, ") * (", b_expr, ")"),
      paste0("TOT_", nm, " := (", dir_expr, ") + (", eg_expr, ") * (", b_expr, ")"),
      paste0("PM_", nm, " := (((", eg_expr, ") * (", b_expr, ")) / ((", dir_expr, ") + (", eg_expr, ") * (", b_expr, "))) * 100")
    )
  }
  
  paste(out, collapse = "\n")
}

build_model_male_b <- function(tipo = c("common", "emb", "main", "full")) {
  tipo <- match.arg(tipo)
  
  parte_a <- '
EG ~ a1*EMB + a2*TX + a3*OB +
     a12*EMB_TX + a13*EMB_OB + a23*TX_OB +
     a123*EMB_TX_OB
'
  
  parte_c <- '
PESO_RN ~ c1*EMB + c2*TX + c3*OB +
          c12*EMB_TX + c13*EMB_OB + c23*TX_OB +
          c123*EMB_TX_OB
'
  
  parte_b <- switch(
    tipo,
    common = '
PESO_RN ~ b0*EG
',
    emb = '
PESO_RN ~ b0*EG + bE*EG_EMB
',
    main = '
PESO_RN ~ b0*EG + bE*EG_EMB + bT*EG_TX + bO*EG_OB
',
    full = '
PESO_RN ~ b0*EG + bE*EG_EMB + bT*EG_TX + bO*EG_OB +
          bET*EG_EMB_TX + bEO*EG_EMB_OB + bTO*EG_TX_OB +
          bETO*EG_EMB_TX_OB
'
  )
  
  contrasts_block <- '
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
'
  
  defined_block <- build_defined_block(tipo)
  
  paste(parte_a, parte_c, parte_b, contrasts_block, defined_block, sep = "\n")
}

fit_sem_safe <- function(model_syntax, data, model_name, bootstrap = 2000, seed = 123) {
  set.seed(seed)
  
  fit <- tryCatch(
    lavaan::sem(
      model = model_syntax,
      data = data,
      meanstructure = TRUE,
      se = "bootstrap",
      bootstrap = bootstrap
    ),
    error = function(e) e
  )
  
  list(
    name = model_name,
    fit = fit,
    is_error = inherits(fit, "error")
  )
}

extract_fit_info <- function(fit_obj) {
  if (fit_obj$is_error) {
    return(tibble(
      modelo = fit_obj$name,
      converged = FALSE,
      npar = NA_integer_,
      chisq = NA_real_,
      df = NA_real_,
      pvalue = NA_real_,
      cfi = NA_real_,
      tli = NA_real_,
      rmsea = NA_real_,
      srmr = NA_real_,
      aic = NA_real_,
      bic = NA_real_,
      error = conditionMessage(fit_obj$fit)
    ))
  }
  
  fit <- fit_obj$fit
  
  tibble(
    modelo = fit_obj$name,
    converged = tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE),
    npar = tryCatch(lavaan::lavInspect(fit, "npar"), error = function(e) NA_integer_),
    chisq = tryCatch(lavaan::fitMeasures(fit, "chisq"), error = function(e) NA_real_),
    df = tryCatch(lavaan::fitMeasures(fit, "df"), error = function(e) NA_real_),
    pvalue = tryCatch(lavaan::fitMeasures(fit, "pvalue"), error = function(e) NA_real_),
    cfi = tryCatch(lavaan::fitMeasures(fit, "cfi"), error = function(e) NA_real_),
    tli = tryCatch(lavaan::fitMeasures(fit, "tli"), error = function(e) NA_real_),
    rmsea = tryCatch(lavaan::fitMeasures(fit, "rmsea"), error = function(e) NA_real_),
    srmr = tryCatch(lavaan::fitMeasures(fit, "srmr"), error = function(e) NA_real_),
    aic = tryCatch(lavaan::fitMeasures(fit, "aic"), error = function(e) NA_real_),
    bic = tryCatch(lavaan::fitMeasures(fit, "bic"), error = function(e) NA_real_),
    error = NA_character_
  )
}

extract_b_table <- function(fit_obj) {
  if (fit_obj$is_error) {
    return(tibble(
      modelo = fit_obj$name,
      lhs = NA_character_,
      op = NA_character_,
      rhs = NA_character_,
      label = NA_character_,
      est = NA_real_,
      se = NA_real_,
      pvalue = NA_real_,
      ci.lower = NA_real_,
      ci.upper = NA_real_,
      std.all = NA_real_
    ))
  }
  
  lavaan::parameterEstimates(
    fit_obj$fit,
    standardized = TRUE,
    ci = TRUE,
    boot.ci.type = "perc"
  ) %>%
    as_tibble() %>%
    filter(
      lhs == "PESO_RN",
      op == "~",
      rhs %in% c("EG", "EG_EMB", "EG_TX", "EG_OB", "EG_EMB_TX", "EG_EMB_OB", "EG_TX_OB", "EG_EMB_TX_OB")
    ) %>%
    mutate(modelo = fit_obj$name) %>%
    select(modelo, lhs, op, rhs, label, est, se, pvalue, ci.lower, ci.upper, std.all)
}

extract_defined_effects_sem <- function(fit_obj) {
  if (fit_obj$is_error) {
    return(tibble(
      modelo = fit_obj$name,
      tipo = NA_character_,
      efecto = NA_character_,
      est = NA_real_,
      se = NA_real_,
      pvalue = NA_real_,
      ci.lower = NA_real_,
      ci.upper = NA_real_
    ))
  }
  
  lavaan::parameterEstimates(
    fit_obj$fit,
    standardized = TRUE,
    ci = TRUE,
    boot.ci.type = "perc"
  ) %>%
    as_tibble() %>%
    filter(op == ":=") %>%
    mutate(
      tipo = case_when(
        grepl("^IND_", label) ~ "Indirecto",
        grepl("^DIR_", label) ~ "Directo",
        grepl("^TOT_", label) ~ "Total",
        grepl("^PM_", label)  ~ "Prop_mediada",
        grepl("^B_", label)   ~ "Pendiente_b",
        TRUE ~ "Otros"
      ),
      modelo = fit_obj$name
    ) %>%
    select(modelo, tipo, efecto = label, est, se, pvalue, ci.lower, ci.upper)
}

model_male_common <- build_model_male_b("common")
model_male_emb    <- build_model_male_b("emb")
model_male_main   <- build_model_male_b("main")
model_male_full   <- build_model_male_b("full")

set.seed(123)

fit_common <- fit_sem_safe(model_male_common, bd_male_sens, "M0_b_comun", bootstrap = 2000, seed = 123)
fit_emb    <- fit_sem_safe(model_male_emb,    bd_male_sens, "M1_b_libre_EMB", bootstrap = 2000, seed = 123)
fit_main   <- fit_sem_safe(model_male_main,   bd_male_sens, "M2_b_libre_EMB_TX_OB", bootstrap = 2000, seed = 123)
fit_full   <- fit_sem_safe(model_male_full,   bd_male_sens, "M3_b_libre_completo", bootstrap = 2000, seed = 123)

fits_list <- list(fit_common, fit_emb, fit_main, fit_full)

table_sens_fit <- purrr::map_dfr(fits_list, extract_fit_info)
table_sens_b <- purrr::map_dfr(fits_list, extract_b_table)
table_sens_defined <- purrr::map_dfr(fits_list, extract_defined_effects_sem)

table_sens_key_effects <- table_sens_defined %>%
  filter(tipo %in% c("Indirecto", "Total", "Pendiente_b"))

cat("\n============================================\n")
cat("MALE PATH b SEM SENSITIVITY\n")
cat("============================================\n")

cat("\n--- Model fit diagnostics ---\n")
print(table_sens_fit, n = Inf)

cat("\n--- Path b coefficients ---\n")
print(table_sens_b, n = Inf)

cat("\n--- Key effects (Indirect, Total, Path b slopes) ---\n")
print(table_sens_key_effects, n = Inf)

# -----------------------------
# 12. Save outputs
# -----------------------------
openxlsx::write.xlsx(
  list(
    Test_bloque_b = table_test_bloque_b,
    Comparacion_OLS = table_comparacion_lm_b,
    Coeficientes_b = table_coeficientes_b,
    Pendientes_b_por_perfil = table_pendientes_b,
    Male_SEM_fit = table_sens_fit,
    Male_SEM_coef_b = table_sens_b,
    Male_SEM_key_effects = table_sens_key_effects,
    Male_SEM_all_defined = table_sens_defined
  ),
  file = file.path(path_output_tables, "Sensitivity_path_b.xlsx"),
  rowNames = FALSE
)

saveRDS(table_test_bloque_b, file = file.path(path_output_rds, "path_b_block_test.rds"))
saveRDS(table_comparacion_lm_b, file = file.path(path_output_rds, "path_b_lm_comparison.rds"))
saveRDS(table_coeficientes_b, file = file.path(path_output_rds, "path_b_boot_coefficients.rds"))
saveRDS(table_pendientes_b, file = file.path(path_output_rds, "path_b_profile_slopes.rds"))
saveRDS(table_sens_fit, file = file.path(path_output_rds, "male_path_b_sem_sensitivity_fit.rds"))
saveRDS(table_sens_b, file = file.path(path_output_rds, "male_path_b_sem_sensitivity_coefficients.rds"))
saveRDS(table_sens_key_effects, file = file.path(path_output_rds, "male_path_b_sem_sensitivity_key_effects.rds"))