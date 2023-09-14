library(tidyverse)

# Race/ethnicity
factor_race = function(var) {
  var = factor(var, levels = c(0, 1, 2, 7), labels =
                 c("Not obtained", "Non-Hispanic White, European ancestry",
                   "Non-Hispanic Black, African ancestry", "Other"))
  var = relevel(var, ref = "Non-Hispanic White, European ancestry")
  var = droplevels(var)
  return(var)
}

# Langa-Weir dementia classification
factor_cog_lw = function(var) {
  var = factor(var, levels = c(1, 2, 3), labels = c("Dementia", "CIND", "Normal"))
  var = fct_infreq(var)
  return(var)
}

# APOE_012 dose variable as factor
factor_APOE_012 = function(var) {
  var = factor(var, levels = c(0, 1, 2), labels = c("Zero", "One", "Two"))
  var = relevel(var, ref = "Zero")
  return(var)
}

# R11DRINK, Currently drinks
factor_ever_drink = function(var) {
  var = factor(var, levels = c(0, 1), labels = c("Doesn't drink", "Drinks"))
  var = relevel(var, ref = "Doesn't drink")
  return(var)
}

# GENDER, sex
factor_gender = function(var) {
  var = factor(var, levels = c(1, 2), labels = c("Male", "Female"))
  var = relevel(var, ref = "Male")
  return(var)
}

# DEGREE, educational attainment
factor_degree = function(var) {
  var = factor(var, levels = c(0, 1, 2, 3, 4, 5, 6, 9), labels =
                 c("No degree", "GED", "High School Diploma", "Two year college degree",
                   "Four year college degree", "Master degree",
                   "Professional Degree", "Degree Unknown/Some college"))
  var = relevel(var, ref = "No degree")
  return(var)
}

# NMARST, marital status at wave N (2012)
factor_marital_status = function(var) {
  var = factor(var, levels = c(1, 2, 3, 4, 5), labels =
                 c("Married", "Separated/Divorced", "Widowed", "Never married",
                   "Unknown"))
  var = relevel(var, ref = "Married")
  var = droplevels(var)
  return(var)
}

# R11SAYRET
factor_says_retired = function(var) {
  var = factor(var, levels = c(0, 1, 2, 3), labels = 
                 c("Not retired", "Completely retired", "Partly retired",
                   "Question irrelevant"))
  var = relevel(var, ref = "Not retired")
  var = droplevels(var)
  return(var)
}

# R11SMOKEV, Ever smoke
factor_ever_smoke = function(var) {
  var = factor(var, levels = c(0, 1), labels = c("Never smoker", "Ever smoker"))
  var = relevel(var, ref = "Never smoker")
  return(var)
}

