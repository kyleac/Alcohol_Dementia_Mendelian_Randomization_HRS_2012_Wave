#Set linear regression tidy model
lm_mod <-
  linear_reg()

#Set logistic regression tidy model
lr_mod <- 
  logistic_reg() %>% 
  set_engine("glm")

### Naive 1SMR Model function ###
#Exposure ~ exposure
#Outcome ~ outcome
#Data ~ dat
run_snp_mr_naive_model <- function(exposure, outcome, dat) {
  
  # Fit the linear model across all instruments, adjusted for ancestry covars
  fit <-
    map(
      .x = snps,
      .f = ~ lm_mod %>%
        fit(formula(paste0(exposure, " ~ ", .x)),
            data = dat)
    )
  
  fit_res <- map2(fit, snps, ~ tidy(.x) %>% filter(term == .y))
  names(fit_res) <- snps
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- snps
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") %>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all.stage1 <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  fit_all.stage1
  
  bx <- fit_all.stage1$estimate
  bxse <- fit_all.stage1$std.error
  
  
  fit <-
    map(
      .x = snps,
      .f = ~ lr_mod %>%
        fit(formula(paste0(outcome, " ~ ", .x, " + ", paste(combined.pgs.covariates, collapse = " + "))),
            data = dat)
    )
  
  fit_res <- map2(fit, snps, ~ tidy(.x, conf.int = TRUE, exponentiate = FALSE) %>% filter(term == .y))
  names(fit_res) <- snps
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- snps
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") #%>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all.stage2 <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  fit_all.stage2
  
  by <- fit_all.stage2$estimate
  byse <- fit_all.stage2$std.error
  
  MRInputObject <- 
    mr_input(bx = bx,
             bxse = bxse,
             by = by,
             byse = byse,
             snps = snps)
  
  MRAllObject_all <- mr_allmethods(MRInputObject, method = "main")
  
  MRAllObject_all
  
  # Interactive plot
  #mr_plot(MRInputObject,
  #error = TRUE, orientate = FALSE, line = "ivw")
  
  print(mr_plot(MRInputObject,
                error = TRUE, orientate = FALSE, line = "ivw", interactive = F))
  return(list(MRAllObject = MRAllObject_all, MRInputObject = MRInputObject, fit.stage1 = fit_all.stage1, fit.stage2 = fit_all.stage2))
}

### Covariate-adjusted 1SMR Model ###
#Exposure ~ exposure
#Outcome ~ outcome
#Data ~ dat
#Covariates ~ covariates
run_snp_mr <- function(exposure, outcome, covariates, dat) {
  
  # Fit the linear model across all instruments, adjusted for ancestry covars
  fit <-
    map(
      .x = snps,
      .f = ~ lm_mod %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(covariates, collapse = " + "))),
            data = dat)
    )
  
  fit_res <- map2(fit, snps, ~ tidy(.x) %>% filter(term == .y))
  names(fit_res) <- snps
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- snps
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") %>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all.stage1 <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  fit_all.stage1
  
  bx <- fit_all.stage1$estimate
  bxse <- fit_all.stage1$std.error
  
  
  fit <-
    map(
      .x = snps,
      .f = ~ lr_mod %>%
        fit(formula(paste0(outcome, " ~ ", .x, " + ", paste(combined.pgs.covariates, collapse = " + "))),
            data = dat)
    )
  
  fit_res <- map2(fit, snps, ~ tidy(.x, conf.int = TRUE, exponentiate = FALSE) %>% filter(term == .y))
  names(fit_res) <- snps
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- snps
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") #%>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all.stage2 <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  fit_all.stage2
  
  by <- fit_all.stage2$estimate
  byse <- fit_all.stage2$std.error
  
  MRInputObject <- 
    mr_input(bx = bx,
             bxse = bxse,
             by = by,
             byse = byse,
             snps = snps)
  
  MRAllObject_all <- mr_allmethods(MRInputObject, method = "main")
  
  MRAllObject_all
  
  # Interactive plot
  #mr_plot(MRInputObject,
  #error = TRUE, orientate = FALSE, line = "ivw")
  
  print(mr_plot(MRInputObject,
                error = TRUE, orientate = FALSE, line = "ivw", interactive = F))
  return(list(MRAllObject = MRAllObject_all, MRInputObject = MRInputObject, fit.stage1 = fit_all.stage1, fit.stage2 = fit_all.stage2))
}

### Covariate-adjusted 1SMR Model w/ SNP option argument ###
# This doesn't work w/ allele scores or only individual instruments b/c requires >2 variants
#Exposure ~ exposure
#Outcome ~ outcome
#Data ~ dat
#Covariates ~ covariates
#Genetic variants of interest ~ snps
run_snp_mr_pgs <- function(exposure, outcome, covariates, dat, snps) {
  
  # Fit the linear model across all instruments, adjusted for ancestry covars
  fit <-
    map(
      .x = snps,
      .f = ~ lm_mod %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(covariates, collapse = " + "))),
            data = dat)
    )
  
  fit_res <- map2(fit, snps, ~ tidy(.x) %>% filter(term == .y))
  names(fit_res) <- snps
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- snps
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") %>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all.stage1 <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  fit_all.stage1
  
  bx <- fit_all.stage1$estimate
  bxse <- fit_all.stage1$std.error
  
  
  fit <-
    map(
      .x = snps,
      .f = ~ lr_mod %>%
        fit(formula(paste0(outcome, " ~ ", .x, " + ", paste(combined.pgs.covariates, collapse = " + "))),
            data = dat)
    )
  
  fit_res <- map2(fit, snps, ~ tidy(.x, conf.int = TRUE, exponentiate = FALSE) %>% filter(term == .y))
  names(fit_res) <- snps
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- snps
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") #%>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all.stage2 <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  fit_all.stage2
  
  by <- fit_all.stage2$estimate
  byse <- fit_all.stage2$std.error
  
  MRInputObject <- 
    mr_input(bx = bx,
             bxse = bxse,
             by = by,
             byse = byse,
             snps = snps)
  
  MRAllObject_all <- mr_allmethods(MRInputObject, method = "main")
  
  MRAllObject_all
  
  # Interactive plot
  #mr_plot(MRInputObject,
  #error = TRUE, orientate = FALSE, line = "ivw")
  
  print(mr_plot(MRInputObject,
                error = TRUE, orientate = FALSE, line = "ivw", interactive = F))
  return(list(MRAllObject = MRAllObject_all, MRInputObject = MRInputObject, fit.stage1 = fit_all.stage1, fit.stage2 = fit_all.stage2))
}

#IVGLM model for both ts and g-est with continuous exposure, continuous IV, and binary outcome, implemented with ivtools package
ivglm_model_exp_num <- function(data, ancestry, variant, exposure, outcome, covariates) {
  
  ## Two-stage estimation
  fitX.LZ <-
    glm(formula = as.formula(
      paste0(
        exposure, " ~ ",
        add.backtick(variant),
        " + ", paste(covariates, collapse = " + ")
      )
    ),
    data = data)
  
  fitY.LX <-
    glm(formula = as.formula(
      paste0(
        outcome,
        " ~ ", exposure, " + ", paste(covariates, collapse = " + ")
      )
    ),
    family = "binomial",
    data = data)
  
  fitIV_ts  <-
    ivglm(
      estmethod = "ts",
      X = exposure,
      Y = outcome,
      fitX.LZ = fitX.LZ,
      fitY.LX = fitY.LX,
      data = data,
      ctrl = F
    )
  
  # summarize ts model
  summary.ts <- summary(fitIV_ts)
  
  # export to results table
  res.ts <- data.frame(
    variant = variant,
    model = "ts",
    ancestry = ancestry,
    term = exposure,
    estimate = round(exp(fitIV_ts$est[exposure]), 2),
    std.error = round(summary.ts$coefficients[2,2], 2),
    statistic = round(summary.ts$coefficients[2,3], 2),
    p.value =  round(summary.ts$coefficients[2,4], 2),
    conf.low = round(exp(confint(fitIV_ts)[2,1]), 2),
    conf.high = round(exp(confint(fitIV_ts)[2,2]), 2),
    covariates = paste(covariates, collapse = " + ")
  )
  
  ## G-estimation
  fitZ.L <-
    glm(formula = as.formula(
      paste0(
        add.backtick(variant),
        " ~ ", paste(covariates, collapse = " + ")
      )
    ),
    data = data)
  
  fitY.LZX <-
    glm(formula = as.formula(
      paste0(
        outcome,
        " ~ ",
        add.backtick(variant),
        " + ", exposure, " + ", paste(covariates, collapse = " + ")
      )
    ),
    family = "binomial",
    data = data)
  
  fitIV_g  <-
    ivglm(
      estmethod = "g",
      X = exposure,
      Y = outcome,
      fitZ.L = fitZ.L,
      fitY.LZX = fitY.LZX,
      data = data,
      link = "logit"
    )
  
  # Summarize g-est model
  summary.g <- summary(fitIV_g)
  
  # Export g-est model results
  res.g <- data.frame(
    variant = variant,
    model = "g",
    ancestry = ancestry,
    term = exposure,
    estimate = round(exp(fitIV_g$est[exposure]), 2),
    std.error = round(summary.g$coefficients[2], 2),
    statistic = round(summary.g$coefficients[3], 2),
    p.value =  round(summary.g$coefficients[4], 2),
    conf.low = round(exp(confint(fitIV_g)[1]), 2),
    conf.high = round(exp(confint(fitIV_g)[2]), 2),
    covariates = paste(covariates, collapse = " + ")
  )
  
  # Concatenate model results
  res <- rbind(res.ts, res.g)
  rownames(res) <- NULL
  
  # Return concatenated results
  return(res)
}

# IVGLM model for both ts and g-est with continuous exposure, continuous IV, and binary outcome, implemented with ivtools package
ivglm_model_exp_num_no_covars <- function(data, ancestry, variant, exposure, outcome) {
  
  ## Two-stage estimation
  fitX.LZ <-
    glm(formula = as.formula(
      paste0(
        exposure, " ~ ",
        add.backtick(variant)
      )
    ),
    data = data)
  
  fitY.LX <-
    glm(formula = as.formula(
      paste0(
        outcome,
        " ~ ", exposure
      )
    ),
    family = "binomial",
    data = data)
  
  fitIV_ts  <-
    ivglm(
      estmethod = "ts",
      X = exposure,
      Y = outcome,
      fitX.LZ = fitX.LZ,
      fitY.LX = fitY.LX,
      data = data,
      ctrl = F
    )
  
  # summarize ts model
  summary.ts <- summary(fitIV_ts)
  
  # export to results table
  res.ts <- data.frame(
    variant = variant,
    model = "ts",
    ancestry = ancestry,
    term = exposure,
    estimate = round(exp(fitIV_ts$est[exposure]), 2),
    std.error = round(summary.ts$coefficients[2,2], 2),
    statistic = round(summary.ts$coefficients[2,3], 2),
    p.value =  round(summary.ts$coefficients[2,4], 2),
    conf.low = round(exp(confint(fitIV_ts)[2,1]), 2),
    conf.high = round(exp(confint(fitIV_ts)[2,2]), 2),
    covariates = "none"
  )
  
  ## G-estimation
  fitZ.L <-
    glm(formula = as.formula(
      paste0(
        add.backtick(variant),
        " ~ ", paste(covariates, collapse = " + ")
      )
    ),
    data = data)
  
  fitY.LZX <-
    glm(formula = as.formula(
      paste0(
        outcome,
        " ~ ",
        add.backtick(variant),
        " + ", exposure, " + ", paste(covariates, collapse = " + ")
      )
    ),
    family = "binomial",
    data = data)
  
  fitIV_g  <-
    ivglm(
      estmethod = "g",
      X = exposure,
      Y = outcome,
      fitZ.L = fitZ.L,
      fitY.LZX = fitY.LZX,
      data = data,
      link = "logit"
    )
  
  # Summarize g-est model
  summary.g <- summary(fitIV_g)
  
  # Export g-est model results
  res.g <- data.frame(
    variant = variant,
    model = "g",
    ancestry = ancestry,
    term = exposure,
    estimate = round(exp(fitIV_g$est[exposure]), 2),
    std.error = round(summary.g$coefficients[2], 2),
    statistic = round(summary.g$coefficients[3], 2),
    p.value =  round(summary.g$coefficients[4], 2),
    conf.low = round(exp(confint(fitIV_g)[1]), 2),
    conf.high = round(exp(confint(fitIV_g)[2]), 2),
    covariates = "none"
  )
  
  # Concatenate model results
  res <- rbind(res.ts, res.g)
  rownames(res) <- NULL
  
  # Return concatenated results
  return(res)
}
