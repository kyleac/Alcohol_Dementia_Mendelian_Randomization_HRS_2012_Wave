#Custom function to tidy parsnip model mapping results apply to all genetic instrument tests (EA_PGS, pgs.norm, and 99 SNPs)
tidy_parsnip_map <- function(fit, map.x) {
  
  fit_res <- map2(fit, genetics, ~ tidy(.x) %>% filter(term == .y))
  names(fit_res) <- genetics
  fit_res_df <- bind_rows(fit_res, .id = "Instrument")
  
  fit_glance <- map(fit, ~ glance(.x))
  names(fit_glance) <- genetics
  fit_glance_df <- bind_rows(fit_glance, .id = "Instrument") %>% rename("F.stat" = statistic, "F.p.value" = p.value)
  
  fit_all <- left_join(fit_res_df, fit_glance_df, by = "Instrument")
  
  return(fit_all)
}

fit_models <- function(exposure, outcome, tidy_model, data) {
  
  # Set the covariate vectors, one model for naive and each of the three combinations below
  combined.pgs.covariates <- c("AncestryPC_1_5A", "AncestryPC_1_5B", "AncestryPC_1_5C", "AncestryPC_1_5D", "AncestryPC_1_5E")
  demographics <- c("NAGE", "GENDER", "DEGREE")
  additional_covars <- c("R11SMOKEV", "NMARST", "R11CESD", "R11CONDE", "R11SAYRET")
  
  # Fit the linear model across all instruments
  fit <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x)),
            data = data)
    )
  
  # Tidy output
  relevance_naive <- tidy_parsnip_map(fit = fit, map.x = instruments)
  
  #Additionally adjusted for ancestry PCs
  fit <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(combined.pgs.covariates, collapse = " + "))),
            data = data)
    )
  
  relevance_ancestryPC <- tidy_parsnip_map(fit = fit, map.x = instruments)
  
  #Additionally adjusted for demographics: age, sex, degree (on top of ancestry PCs)
  fit <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(c(demographics, combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  relevance_demographics <- tidy_parsnip_map(fit = fit, map.x = instruments)
  
  # Additionally adjusted for APOE012 (on top of ancestry PCs and demographics)
  fit <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(c(demographics, "APOE012", combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  relevance_demographics_apoe <- tidy_parsnip_map(fit = fit, map.x = instruments)
  
  #Additionally adjusted for all potential confounders identified in cross-sectional analysis lit. review
  # Fit the linear model across all instruments
  fit <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(c(demographics, "APOE012", additional_covars, combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  relevance_full <- tidy_parsnip_map(fit = fit, map.x = instruments)
  
  #Helpful webpage in doing ANOVA in tidymodels: https://www.thomasvanhoey.com/post/2021-10-12-tidymodels-interactions/
  #anova(reduced, full)
  fit_ancestry_reduced <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", paste(combined.pgs.covariates, collapse = " + "))),
            data = data)
    )
  
  fit_ancestry_full <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(combined.pgs.covariates, collapse = " + "))),
            data = data)
    )
  
  fit_ancestry_anova <- map2(.x = fit_ancestry_reduced,
                             .y = fit_ancestry_full,
                             ~ anova(.x %>% extract_fit_engine, .y %>% extract_fit_engine))
  names(fit_ancestry_anova) <- instruments
  
  fit_demo_reduced <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", paste(c(demographics, combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  fit_demo_full <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(c(demographics, combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  fit_demo_anova <- map2(.x = fit_demo_reduced,
                         .y = fit_demo_full,
                         ~ anova(.x %>% extract_fit_engine, .y %>% extract_fit_engine))
  names(fit_demo_anova) <- instruments
  
  fit_demo_apoe_reduced <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", paste(c(demographics, "APOE012", combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  fit_demo_apoe_full <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(c(demographics, "APOE012", combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  fit_demo_apoe_anova <- map2(.x = fit_demo_apoe_reduced,
                              .y = fit_demo_apoe_full,
                              ~ anova(.x %>% extract_fit_engine, .y %>% extract_fit_engine))
  names(fit_demo_anova) <- instruments
  
  fit_full_reduced <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", paste(c(demographics, "APOE012", additional_covars, combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  fit_full_full <-
    map(
      .x = instruments,
      .f = ~ tidy_model %>%
        fit(formula(paste0(exposure, " ~ ", .x, " + ", paste(c(demographics, "APOE012", additional_covars, combined.pgs.covariates), collapse = " + "))),
            data = data)
    )
  
  fit_full_anova <- map2(.x = fit_full_reduced,
                         .y = fit_full_full,
                         ~ anova(.x %>% extract_fit_engine, .y %>% extract_fit_engine))
  names(fit_full_anova) <- instruments
  
  #Retrieve partial F statistics
  x <- lapply(fit_ancestry_anova, function(x) {x$`F`[2]}) %>% unlist
  
  anova_ancestry_df <- data.frame(
    Instrument = instruments,
    Model = "AncestryPCs",
    `Partial F` = x
  )
  
  x <- lapply(fit_demo_anova, function(x) {x$`F`[2]}) %>% unlist
  
  anova_demo_df <- data.frame(
    Instrument = instruments,
    Model = "Demographics",
    `Partial F` = x
  )
  
  x <- lapply(fit_demo_apoe_anova, function(x) {x$`F`[2]}) %>% unlist
  
  anova_demo_apoe_df <- data.frame(
    Instrument = instruments,
    Model = "Demographics, APOE",
    `Partial F` = x
  )  
  
  x <- lapply(fit_full_anova, function(x) {x$`F`[2]}) %>% unlist
  
  anova_full_df <- data.frame(
    Instrument = instruments,
    Model = "Full",
    `Partial F` = x
  )
  
  partialf <-
    rbind(
      anova_ancestry_df,
      anova_demo_df,
      anova_demo_apoe_df,
      anova_full_df
    )
  
  #Collate instrument relevance results from EA
  relevance <- 
    rbind(
      relevance_naive %>% dplyr::select(Instrument, estimate, std.error, F.stat, F.p.value) %>% mutate(Model = "Unadjusted"),
      relevance_ancestryPC %>% dplyr::select(Instrument, estimate, std.error, F.stat, F.p.value) %>% mutate(Model = "AncestryPCs"),
      relevance_demographics %>% dplyr::select(Instrument, estimate, std.error, F.stat, F.p.value) %>% mutate(Model = "Demographics"),
      relevance_demographics_apoe %>% dplyr::select(Instrument, estimate, std.error, F.stat, F.p.value) %>% mutate(Model = "Demographics, APOE"),
      relevance_full %>% dplyr::select(Instrument, estimate, std.error, F.stat, F.p.value) %>% mutate(Model = "Full")
    )
  
  return <- left_join(relevance, partialf) %>%
    mutate(graph_F = ifelse(
      test = is.na(Partial.F),
      F.stat,
      Partial.F
    )) %>%
    mutate(Model = fct_relevel(Model, "Unadjusted")) %>%
    arrange(desc(Instrument), Model)
  
  return(return)
}