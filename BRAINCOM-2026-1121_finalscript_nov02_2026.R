library(dplyr)
library(ggplot2)
library(car)
library(effectsize)
library(WRS2)
library(broom)
library(tidyr)
library(nnet)
library(tidyverse)
library(glmmTMB)
library(patchwork)
library(ggtext)
library(ggnewscale)
library(purrr)


# Load data_final_nov04_2026
data1 <- data_final_nov04_2026

# ---- Construct table 1 ----

## ---- Look at distribution of continuous variables to see if they need to be tranformed ----
continuous <- c(
  "time_since_injury_months", "neuro_yrs_education_2",
  "demo_1_age_3", "mmq_total", "rpq_total",
  "fmdi_total", "sf36_total", "whodas_total_rasch"
)

# histograms
par(mfrow = c(2, 4))

for (v in continuous) {
  hist(data1[[v]],
       main = v,
       xlab = "Score",
       col = "steelblue",
       breaks = 15)
}

par(mfrow = c(1, 1))

# logarithmic transformation for whodas_total_rasch and demo_1_age_3 (log1p(x) calculates log(x+1) so it can handle 0s)
data1$whodas_log <- log1p(data1$whodas_total_rasch)
data1$age_log <- log1p(data1$demo_1_age_3)
#  histograms
par(mfrow = c(1, 2))
hist(data1$whodas_log, main = "WHODAS (Log)", col = "lightgreen")
hist(data1$age_log, main = "Age (Log)", col = "lightgreen")

# make sure nothing is missing
sum(!is.na(data1$whodas_total_rasch) & is.na(data1$whodas_log))
sum(!is.na(data1$demo_1_age_3) & is.na(data1$age_log))
# verify transformations are correct
all.equal(data1$whodas_log, log1p(data1$whodas_total_rasch))
all.equal(data1$age_log, log1p(data1$demo_1_age_3))

## ---- Now create table 1 & 2 ----
continuous <- c(
  "time_since_injury_months", "neuro_yrs_education_2",
  "demo_1_age_3", "mmq_total", "rpq_total", "fmdi_total",
  "sf36_total", "whodas_total_rasch"
)

categorical <- c(
  "moi", "elig_screen_15_conciousness_eligiblity_screen_arm_1_M4",
  "elig_screen_18_gap_eligiblity_screen_arm_1_M4",
  "benefits_recoded", "demo_2_sex_3", "demo_4_gender_3",
  "ethnicity", "work_status_num", "demo_9_preconcussion_3",
  "demo_depression_3", "demo_anxiety_3", "demo_somatic_3"
)

transformed <- c(
  whodas_total_rasch = "whodas_log",
  demo_1_age_3 = "age_log"
)

families <- list(
  c("mmq_total", "rpq_total", "fmdi_total",
    "sf36_total", "whodas_total_rasch"),
  c(categorical[1:4], "time_since_injury_months"),
  c(categorical[5:12], "neuro_yrs_education_2", "demo_1_age_3")
)

data1 <- data1 %>%
  mutate(
    group = relevel(factor(group_assignment_1), ref = "2"),
    across(all_of(continuous), ~as.numeric(as.character(.x))),
    across(all_of(categorical), ~factor(as.character(.x)))
  )

# categories retained in regression model
keep <- list(
  demo_4_gender_3 = c("1", "2"),
  ethnicity = c("1", "3", "10"),
  work_status_num = c("1", "3", "4"),
  moi = c("1", "2", "4", "5", "6")
)

analyze <- function(v) {
  
  model_var <- if (v %in% names(transformed)) transformed[[v]] else v
  is_cat <- v %in% categorical
  
  d <- data1 %>%
    transmute(group, raw = .data[[v]], y = .data[[model_var]])
  
  d_desc <- d %>% drop_na(group, raw)
  desc <- d_desc %>%
    group_by(group) %>%
    summarise(
      summary = if (is_cat) {
        t <- table(raw)
        paste(names(t), sprintf("%d (%.1f%%)",
                                as.integer(t), 100 * prop.table(t)), collapse = "; ")
      } else {
        sprintf("%.2f (%.2f)", mean(raw), sd(raw))
      },
      .groups = "drop"
    )
  
  if (v == categorical[2])
    d$y <- dplyr::recode(as.character(d$y), "2" = "1", "3" = "0")
  
  if (v %in% names(keep))
    d <- d %>% filter(as.character(y) %in% keep[[v]])
  
  d <- d %>% drop_na(group, raw, y)
  d$group <- droplevels(d$group)
  if (is_cat) d$y <- droplevels(factor(d$y))
  
  if (is_cat && nlevels(d$y) > 2) {
    model <- "multinomial"
    full <- multinom(y ~ group, data = d, trace = FALSE)
    null <- multinom(y ~ 1, data = d, trace = FALSE)
  } else {
    model <- if (is_cat) "binomial" else "gaussian"
    fam <- if (is_cat) binomial() else gaussian()
    full <- glm(y ~ group, data = d, family = fam)
    null <- glm(y ~ 1, data = d, family = fam)
  }
  
  LR <- 2 * as.numeric(logLik(full) - logLik(null))
  df <- attr(logLik(full), "df") - attr(logLik(null), "df")
  
  if (!is_cat) {
    x <- d$y[d$group == "1"]
    y <- d$y[d$group == "2"]
    nx <- length(x); ny <- length(y)
    
    sp <- sqrt(((nx-1)*var(x) + (ny-1)*var(y))/(nx+ny-2))
    smd <- (mean(x)-mean(y))/sp
    se <- sqrt((nx+ny)/(nx*ny) + smd^2/(2*(nx+ny)))
    
    effect <- sprintf("SMD %.3f [%.3f, %.3f]",
                      smd, smd-1.96*se, smd+1.96*se)
  } else {
    b <- broom::tidy(full) %>% filter(grepl("group", term))
    
    if (model == "multinomial")
      b$term <- paste(b$y.level, b$term, sep = ": ")
    
    effect <- paste(sprintf("%s: OR %.2f [%.2f, %.2f]",
                            b$term, exp(b$estimate),
                            exp(b$estimate-1.96*b$std.error),
                            exp(b$estimate+1.96*b$std.error)), collapse = "; ")
  }
  
  tibble(
    variable = v,
    model = model,
    FCD = desc$summary[match("1", desc$group)],
    Control = desc$summary[match("2", desc$group)],
    LR_X2 = LR,
    df = df,
    p = pchisq(LR, df, lower.tail = FALSE),
    effect = effect,
    n_descriptive = nrow(d_desc),
    n_model = nrow(d),
    n_missing = sum(is.na(data1[[v]])),
    n_excluded = nrow(data1) - nrow(d) -
      sum(is.na(data1[[v]]) | is.na(data1$group))
  )
}

#run analyses and apply Holm corrections
results <- map_dfr(c(continuous, categorical), analyze) %>%
  mutate(
    family = map_int(variable,
                     ~which(map_lgl(families, \(f) .x %in% f)))
  ) %>%
  group_by(family) %>%
  mutate(holm_p = p.adjust(p, method = "holm")) %>%
  ungroup()

results <- results %>%
  mutate(LR_X2 = sprintf("%.2f", LR_X2))

print(results, n = Inf, width = Inf)

# ---- Neuropsychological tests ----
graphics.off()

neuropsych <- c(
  "topf_score_2", "neuro_ravlt_dr_cr_1_2",
  "neuro_ravlt_dr_cr_2_2", "neuro_ravlt_recog_score_2",
  "tmt_time_a_2", "tmt_time_b_2",
  "wais_ds_score_2", "wais_ss_score_2",
  "wais_c_score_2", "vf_total_correct_2"
)

## ---- print histograms to see if some values need to be transformed ----
data1 |>
  select(all_of(neuropsych)) |>
  pivot_longer(everything()) |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 15, fill = "grey70", color = "white") +
  facet_wrap(~name, scales = "free") +
  theme_minimal()

# best transformation approach for neuro_ravlt_dr_cr_1_2
x <- data1$neuro_ravlt_dr_cr_1_2
par(mfrow = c(2, 2))
hist(x, main = "Raw")
hist(log1p(x), main = "Log")
hist(x^(1/3), main = "Cube root")
hist(x^2, main = "Squared")
par(mfrow = c(1, 1))

# best transformation approach for neuro_ravlt_dr_cr_2_2
x <- data1$neuro_ravlt_dr_cr_2_2
par(mfrow = c(2, 2))
hist(x, main = "Raw")
hist(log1p(x), main = "Log")
hist(x^(1/3), main = "Cube root")
hist(x^2, main = "Squared")
par(mfrow = c(1, 1))

# squared for both
data1 <- data1 %>%
  mutate(
    ravlt_ir_squared = neuro_ravlt_dr_cr_1_2^2,
    ravlt_dr_squared = neuro_ravlt_dr_cr_2_2^2
  )
all.equal(data1$ravlt_ir_squared, data1$neuro_ravlt_dr_cr_1_2^2)
all.equal(data1$ravlt_dr_squared, data1$neuro_ravlt_dr_cr_2_2^2)
colSums(is.na(data1[c(
  "neuro_ravlt_dr_cr_1_2", "ravlt_ir_squared",
  "neuro_ravlt_dr_cr_2_2", "ravlt_dr_squared" )]))

## ---- Look at group comparisons ---- 

transformed <- c(
  neuro_ravlt_dr_cr_1_2 = "ravlt_ir_squared",
  neuro_ravlt_dr_cr_2_2 = "ravlt_dr_squared"
)

inverse_gaussian <- c(
  "tmt_time_a_2", "tmt_time_b_2",
  "wais_ds_score_2", "wais_ss_score_2", "wais_c_score_2"
)

covariates <- c("demo_1_age_3", "neuro_yrs_education_2", "demo_2_sex_3")

data1 <- data1 %>%
  mutate(
    group = factor(group_assignment_1, levels = c("2", "1")),
    across(all_of(c(neuropsych, covariates[1:2])),
           ~as.numeric(as.character(.x))),
    demo_2_sex_3 = factor(demo_2_sex_3)
  )

analyze <- function(v) {
  
  model_var <- if (v %in% names(transformed)) transformed[[v]] else v
  is_ig <- v %in% inverse_gaussian
  
  d <- data1 %>%
    transmute(group, raw = .data[[v]], y = .data[[model_var]],
              across(all_of(covariates))) %>%
    drop_na() %>% droplevels()
  
  fam <- if (is_ig) inverse.gaussian(link = "log") else gaussian()
  full <- glm(reformulate(c("group", covariates), "y"), d, family = fam)
  null <- glm(reformulate(covariates, "y"), d, family = fam)
  
  LR <- 2 * as.numeric(logLik(full) - logLik(null))
  df <- attr(logLik(full), "df") - attr(logLik(null), "df")
  
  x <- d$y[d$group == "1"]
  y <- d$y[d$group == "2"]
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx-1)*var(x) + (ny-1)*var(y))/(nx+ny-2))
  
  if (!is_ig) {
    est <- coef(full)["group1"]
    se <- sqrt(vcov(full)["group1", "group1"])
  } else {
    b <- coef(full)
    ctrl <- d
    ctrl$group <- factor("2", levels = c("2", "1"))
    X <- model.matrix(delete.response(terms(full)), ctrl)
    X <- X[, names(b), drop = FALSE]
    
    mu0 <- as.numeric(exp(X %*% b))
    mu1 <- mu0 * exp(b["group1"])
    diff <- mu1 - mu0
    est <- mean(diff)
    
    e <- as.numeric(names(b) == "group1")
    grad <- colMeans(sweep(X, 1, diff, "*") + outer(mu1, e))
    se <- sqrt(as.numeric(t(grad) %*% vcov(full) %*% grad))
  }
  
  tibble(
    variable = v,
    model_variable = model_var,
    family = if (is_ig) "inverse gaussian" else "gaussian",
    FCD = sprintf("%.2f (%.2f)",
                  mean(d$raw[d$group == "1"]), sd(d$raw[d$group == "1"])),
    Control = sprintf("%.2f (%.2f)",
                      mean(d$raw[d$group == "2"]), sd(d$raw[d$group == "2"])),
    LR_X2 = LR,
    df = df,
    p = pchisq(LR, df, lower.tail = FALSE),
    SMD = est/sp,
    CI_lower = (est-1.96*se)/sp,
    CI_upper = (est+1.96*se)/sp,
    n_model = nrow(d),
    n_missing = sum(is.na(data1[[v]]))
  )
}

neuro_results <- map_dfr(neuropsych, analyze) %>%
  mutate(holm_p = p.adjust(p, method = "holm"))

neuro_table <- neuro_results %>%
  mutate(
    LR_X2 = sprintf("%.3f", LR_X2),
    across(c(p, holm_p), ~sprintf("%.12f", .x)),
    across(c(SMD, CI_lower, CI_upper), ~sprintf("%.3f", .x))
  )

print(neuro_table, n = Inf, width = Inf)

## ---- Number of participants with "cognitive impairment" - per https://pmc.ncbi.nlm.nih.gov/articles/PMC8967334/ ----
vars <- c("ravlt_ir_percentile_2","ravlt_dr_percentile_2",
          "tmt_a_percentile_2","tmt_b_percentile_2",
          "wais_ps_percentile_2")

d <- data1 %>%
  mutate(
    group = factor(group_assignment_1, levels=c("2","1")),
    across(all_of(vars), ~.x <= 9),
    ravlt_impaired = ravlt_ir_percentile_2 | ravlt_dr_percentile_2,
    tmt_impaired = tmt_a_percentile_2 | tmt_b_percentile_2,
    wais_impaired = wais_ps_percentile_2,
    domain_count = rowSums(across(c(ravlt_impaired,tmt_impaired,wais_impaired)),
                           na.rm=TRUE),
    any_impaired = case_when(
      domain_count >= 1 ~ TRUE,
      if_all(all_of(vars), ~!is.na(.x)) ~ FALSE,
      TRUE ~ NA
    ),
    two_domains = case_when(
      domain_count >= 2 ~ TRUE,
      rowSums(!is.na(across(c(ravlt_impaired,tmt_impaired,wais_impaired)))) == 3 ~ FALSE,
      TRUE ~ NA
    )
  )

for(v in vars) {
  cat("\n",v,"\n")
  print(table(d$group,d[[v]],useNA="always"))
}

for(v in c("any_impaired","two_domains","domain_count")) {
  cat("\n",v,"\n")
  print(table(d$group,d[[v]],useNA="always"))
}

# ---- Behav obs ----
obs_vars <- c("neuro_affect_2", "neuro_response_2",
              "neuro_aware_2", "neuro_verb_2")

covariates <- c("demo_1_age_3", "neuro_yrs_education_2", "demo_2_sex_3")

keep <- list(
  neuro_affect_2 = c("1", "2"),
  neuro_aware_2 = c("2", "3")
)

refs <- c(neuro_aware_2 = "2", neuro_verb_2 = "0")

data1 <- data1 %>%
  mutate(
    group = factor(group_assignment_1, levels = c("2", "1")),
    across(all_of(obs_vars), ~factor(as.character(.x))),
    demo_2_sex_3 = factor(demo_2_sex_3)
  )

analyze_obs <- function(v) {
  
  d <- data1 %>%
    transmute(group, y = .data[[v]], across(all_of(covariates)))
  
  desc <- d %>%
    drop_na(group, y) %>%
    group_by(group) %>%
    summarise(
      summary = {
        t <- table(y)
        paste(names(t), sprintf("%d (%.1f%%)",
                                as.integer(t), 100 * prop.table(t)), collapse = "; ")
      },
      .groups = "drop"
    )
  
  if (v %in% names(keep))
    d <- d %>% filter(as.character(y) %in% keep[[v]])
  
  d <- d %>% drop_na() %>% droplevels()
  
  if (v %in% names(refs))
    d$y <- relevel(d$y, ref = refs[[v]])
  
  is_multi <- nlevels(d$y) > 2
  
  f1 <- reformulate(c("group", covariates), "y")
  f0 <- reformulate(covariates, "y")
  
  if (is_multi) {
    full <- multinom(f1, data = d, trace = FALSE)
    null <- multinom(f0, data = d, trace = FALSE)
  } else {
    full <- glm(f1, data = d, family = binomial())
    null <- glm(f0, data = d, family = binomial())
  }
  
  LR <- 2 * as.numeric(logLik(full) - logLik(null))
  df <- attr(logLik(full), "df") - attr(logLik(null), "df")
  
  b <- broom::tidy(full) %>%
    filter(term == "group1")
  
  if (is_multi)
    b$term <- paste(b$y.level, "vs", levels(d$y)[1])
  
  OR <- paste(sprintf("%s: %.2f [%.2f, %.2f]",
                      b$term, exp(b$estimate),
                      exp(b$estimate - 1.96*b$std.error),
                      exp(b$estimate + 1.96*b$std.error)),
              collapse = "; ")
  
  tibble(
    variable = v,
    family = if (is_multi) "multinomial" else "binomial",
    FCD = desc$summary[match("1", desc$group)],
    Control = desc$summary[match("2", desc$group)],
    LR_X2 = LR,
    df = df,
    p = pchisq(LR, df, lower.tail = FALSE),
    OR_95CI = OR,
    n_model = nrow(d),
    n_missing = sum(is.na(data1[[v]]))
  )
}

# analyses and Holm correction
obs_results <- map_dfr(obs_vars, analyze_obs) %>%
  mutate(holm_p = p.adjust(p, method = "holm"))

obs_table <- obs_results %>%
  mutate(
    LR_X2 = sprintf("%.3f", LR_X2),
    across(c(p, holm_p), ~sprintf("%.12f", .x))
  )

print(obs_table, n = Inf, width = Inf)

# ---- Questionnaires ----

## ---- Between group comparisons ----

### ---- First, look at histograms ----
vars <- c("fam_total","mcq_total","pain_total","gad7_total","phq9_total",
          "mia_total","phq15_total","isi_total","des_total","tas_total",
          "soms7_total","cdrisc_total","ssd_total","rts_rnt_total","lines_pediatric")

data1 %>%
  pivot_longer(all_of(vars), names_to = "Questionnaire", values_to = "Score") %>%
  ggplot(aes(Score)) +
  geom_histogram(bins = 15, fill = "steelblue", colour = "white") +
  facet_wrap(~Questionnaire, scales = "free", ncol = 3) +
  theme_bw()

# Based on visual inspection, variables to transform include soms7_total, pain_total, and ssd_total
# best transformation approach for soms7_total
x <- data1$soms7_total
par(mfrow = c(2, 2))
hist(x, main = "Raw")
hist(log1p(x), main = "Log")
hist(x^(1/3), main = "Cube root")
hist(x^2, main = "Squared")
par(mfrow = c(1, 1))
#best transformation approach for catastrophizing
x <- data1$pain_total
par(mfrow = c(2, 2))
hist(x, main = "Raw")
hist(log1p(x), main = "Log")
hist(x^(1/3), main = "Cube root")
hist(x^2, main = "Squared")
par(mfrow = c(1, 1))
# best transformation approach for ssd_total
x <- data1$ssd_total
par(mfrow = c(2, 2))
hist(x, main = "Raw")
hist(log1p(x), main = "Log")
hist(x^(1/3), main = "Cube root")
hist(x^2, main = "Squared")
par(mfrow = c(1, 1))
# Choose Log for all

# logarithmic transformations
data1$soms7_total_log <- log1p(data1$soms7_total)
data1$pain_total_log <- log1p(data1$pain_total)
data1$ssd_total_log <- log1p(data1$ssd_total)


par(mfrow = c(1, 3))
hist(data1$soms7_total_log, main = "SOMS7 (Log)", col = "lightgreen")
hist(data1$pain_total_log, main = "Catastrophizing (Log)", col = "lightgreen")
hist(data1$ssd_total_log, main = "SSD (Log)", col = "lightgreen")
# make sure nothing is missing
sum(!is.na(data1$soms7_total) & is.na(data1$soms7_total_log))
sum(!is.na(data1$pain_total) & is.na(data1$pain_total_log))
sum(!is.na(data1$ssd_total) & is.na(data1$ssd_total_log))
# verify transformations are correct
all.equal(data1$soms7_total_log, log1p(data1$soms7_total))
all.equal(data1$pain_total_log, log1p(data1$pain_total))
all.equal(data1$ssd_total_log, log1p(data1$ssd_total))

### ---- run GLMs 1 ----
vars <- c("fam_total","mcq_total","pain_total","gad7_total","phq9_total",
          "mia_total","phq15_total","isi_total","des_total","tas_total",
          "soms7_total","cdrisc_total","ssd_total","rts_rnt_total","lines_pediatric")
logs <- c("soms7_total","pain_total","ssd_total")
covars <- c("demo_1_age_3","neuro_yrs_education_2","demo_2_sex_3")
holm <- list(
  c("cdrisc_total","lines_pediatric","mia_total","tas_total"),
  c("des_total","fam_total","mcq_total","pain_total","rts_rnt_total"),
  c("gad7_total","isi_total","phq9_total","phq15_total","soms7_total","ssd_total")
)

data1 <- data1 %>%
  mutate(group = factor(group_assignment_1, levels=c("2","1")),
         demo_2_sex_3 = factor(demo_2_sex_3),
         across(all_of(logs), log1p, .names="{.col}_log"))

analyze <- function(v) {
  outcome <- if (v %in% logs) paste0(v,"_log") else v
  d <- data1 %>% select(group, all_of(c(v,outcome,covars))) %>% drop_na()
  d$y <- d[[outcome]]
  
  families <- list(gaussian=gaussian())
  if (all(d$y > 0))
    families <- c(families, list(
      gamma=Gamma(link="log"),
      inverse_gaussian=inverse.gaussian(link="log")))
  if (all(d$y >= 0))
    families$tweedie <- tweedie(link="log")
  
  f <- reformulate(c("group",covars),"y")
  f0 <- reformulate(covars,"y")
  
  fit <- function(formula, fam, name) {
    tryCatch(suppressWarnings(
      if (name=="tweedie")
        glmmTMB(formula, family=fam, data=d)
      else
        glm(formula, family=fam, data=d)
    ), error=function(e) NULL)
  }
  
  valid <- function(m) {
    if (is.null(m)) return(FALSE)
    ok <- if (inherits(m,"glmmTMB"))
      isTRUE(m$fit$convergence==0) && isTRUE(m$sdr$pdHess)
    else isTRUE(m$converged)
    ok && is.finite(AIC(m))
  }
  
  fits <- imap(families, ~fit(f,.x,.y)) %>% keep(valid)
  if (!length(fits)) stop(paste("No valid model:",v))
  
  best <- names(which.min(map_dbl(fits,AIC)))
  full <- fits[[best]]
  null <- fit(f0,families[[best]],best)
  if (!valid(null)) stop(paste("Reduced model failed:",v))
  
  LR <- max(0,2*as.numeric(logLik(full)-logLik(null)))
  df <- attr(logLik(full),"df")-attr(logLik(null),"df")
  
  d0 <- d1 <- d
  d0$group <- factor("2",levels=c("2","1"))
  d1$group <- factor("1",levels=c("2","1"))
  
  b <- if (best=="tweedie") fixef(full)$cond else coef(full)
  V <- if (best=="tweedie") vcov(full)$cond else vcov(full)
  V <- as.matrix(V)
  
  X0 <- model.matrix(f,d0)[,names(b),drop=FALSE]
  X1 <- model.matrix(f,d1)[,names(b),drop=FALSE]
  eta0 <- drop(X0%*%b)
  eta1 <- drop(X1%*%b)
  
  if (best=="gaussian") {
    effect <- mean(eta1-eta0)
    grad <- colMeans(X1-X0)
  } else {
    mu0 <- exp(eta0); mu1 <- exp(eta1)
    effect <- mean(mu1-mu0)
    grad <- colMeans(sweep(X1,1,mu1,"*")-
                       sweep(X0,1,mu0,"*"))
  }
  
  se <- sqrt(drop(t(grad)%*%V%*%grad))
  x <- d$y[d$group=="1"]
  y <- d$y[d$group=="2"]
  sp <- sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y)) /
               (length(x)+length(y)-2))
  fmt <- function(x) sprintf("%.2f (%.2f)",mean(x),sd(x))
  
  tibble(
    variable=v,
    transformation=ifelse(v %in% logs,"Log1p","Raw"),
    family=best,
    FCD=fmt(d[[v]][d$group=="1"]),
    Control=fmt(d[[v]][d$group=="2"]),
    LR_X2=LR, df=df,
    p=pchisq(LR,df,lower.tail=FALSE),
    SMD=effect/sp,
    CI_lower=(effect-1.96*se)/sp,
    CI_upper=(effect+1.96*se)/sp,
    n_model=nrow(d),
    n_missing=sum(is.na(data1[[v]]))
  )
}

results <- map_dfr(vars,analyze) %>%
  mutate(holm_family=map_int(variable,
                             ~which(map_lgl(holm,\(f) .x %in% f)))) %>%
  group_by(holm_family) %>%
  mutate(holm_p=p.adjust(p,"holm")) %>%
  ungroup()

table <- results %>%
  mutate(
    LR_X2=sprintf("%.3f",LR_X2),
    across(c(p,holm_p),~sprintf("%.12f",.x)),
    across(c(SMD,CI_lower,CI_upper),~sprintf("%.3f",.x))
  )

print(table,n=Inf,width=Inf)

### ---- Supplementary Figure 1 ----

#Qs
vars <- c("cdrisc_total","lines_pediatric","mia_total","tas_total",
          "des_total","fam_total","mcq_total","pain_total","rts_rnt_total",
          "gad7_total","isi_total","phq9_total","phq15_total","soms7_total","ssd_total")

labels <- c(
  "CD-RISC Total","LiNES Developmental","MIA-A Total","TAS-20 Total",
  "DES-II Total","FAM Total","MCQ Total","SCS Total","RTQ-RNT Total",
  "GAD-7 Total","ISI Total","PHQ-9 Total","PHQ-15 Total",
  "SOMS-7 Total","SSD-12 Total"
)
names(labels) <- vars

# prepare data
d <- data1 %>%
  mutate(group_label = factor(group_assignment_1,
                              levels=c("1","2"),
                              labels=c("FCD","Controls")))

plots <- lapply(vars, function(v) {
  x <- d[[v]]
  bw <- diff(range(x, na.rm = TRUE))/20
  if (bw == 0) bw <- 1
  breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE) + bw, by = bw)
  
  h <- hist(d[[v]][d$group_label == "FCD"],
            breaks = breaks, plot = FALSE)
  
  fcd <- tibble(
    x = h$breaks,
    y = c(h$counts, 0)
  )
  
  ggplot() +
    geom_histogram(
      data = filter(d, group_label == "Controls"),
      aes(x = .data[[v]]),
      breaks = breaks, fill = "#a8a8a8", alpha = .7
    ) +
    geom_step(
      data = fcd, aes(x, y),
      direction = "hv", colour = "black", linewidth = .7
    ) +
    labs(x = labels[[v]], y = "Frequency") +
    theme_classic(base_size = 9) +
    theme(
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
})

final <- wrap_plots(plots, ncol = 3) +
  plot_annotation(caption = "Black line = FCD     |     Grey bars = Controls")

print(final)

## --- Correlations within FCD group ----
### ---- Also manuscript Figure 2 ----

# Variables and labels
targets <- c(MMQ_S = "mmq_total", FMDI = "fmdi_total")

sections <- list(
  "Symptom and functional measures" = c("RPQ","SF-36","WHODAS"),
  "Candidate predisposing factors" = c("CD-RISC","LiNES-Developmental","MIA-A","TAS-20"),
  "Candidate perpetuating factors" = c("DES-II","FAM","MCQ","SCS","RTS-RNT"),
  "Comorbid factors" = c("GAD-7","ISI","PHQ-9","PHQ-15","SOMS-7","SSD-12")
)

labels <- c(
  rpq_total="RPQ", sf36_total="SF-36", whodas_total_rasch="WHODAS",
  cdrisc_total="CD-RISC", lines_pediatric="LiNES-Developmental",
  mia_total="MIA-A", tas_total="TAS-20", des_total="DES-II",
  fam_total="FAM", mcq_total="MCQ", pain_total="SCS",
  rts_rnt_total="RTS-RNT", gad7_total="GAD-7", isi_total="ISI",
  phq9_total="PHQ-9", phq15_total="PHQ-15",
  soms7_total="SOMS-7", ssd_total="SSD-12"
)

constructs <- c(
  RPQ="Post-concussion symptoms", `SF-36`="Health-related quality of life",
  WHODAS="Disability", `CD-RISC`="Resilience",
  `LiNES-Developmental`="Attachment insecurity",
  `MIA-A`="Memory perfectionism", `TAS-20`="Alexithymia",
  `DES-II`="Dissociation", FAM="Avoidance behaviour",
  MCQ="Use of compensatory strategies", SCS="Symptom catastrophizing",
  `RTS-RNT`="Rumination", `GAD-7`="Anxiety symptoms",
  ISI="Insomnia symptoms", `PHQ-9`="Depression symptoms",
  `PHQ-15`="Somatic symptoms",
  `SOMS-7`="Functional sensory-motor symptoms",
  `SSD-12`="SSD DSM-5 B-criteria"
)

# spearman
set.seed(12345)
d <- data1 %>% filter(as.character(group_assignment_1) == "1")

cors <- map_dfr(names(targets), function(t) {
  map_dfr(names(labels), function(v) {
    x <- d %>% select(all_of(c(targets[[t]], v))) %>% drop_na()
    a <- x[[1]]; b <- x[[2]]
    test <- suppressWarnings(cor.test(a, b, method="spearman", exact=FALSE))
    
    boot <- replicate(5000, {
      i <- sample.int(length(a), replace=TRUE)
      suppressWarnings(cor(a[i], b[i], method="spearman"))
    })
    
    ci <- quantile(boot, c(.025,.975), na.rm=TRUE)
    tibble(target=t, variable=labels[[v]], n=length(a),
           r=unname(test$estimate), low=ci[1], high=ci[2],
           p=test$p.value)
  }) %>% mutate(p_holm=p.adjust(p,"holm"))
})

order <- unlist(sections, use.names=FALSE)
lookup <- setNames(rep(names(sections), lengths(sections)), order)

hm <- cors %>%
  mutate(
    section=lookup[variable],
    target=factor(target, levels=c("MMQ_S","FMDI"),
                  labels=c("MMQ-S","FMDI")),
    label=sprintf("%.2f\n[%.2f, %.2f]",r,low,high),
    label=gsub("(-?)0\\.","\\1.",label)
  )

# plot each section
make_section <- function(title, vars) {
  x <- hm %>%
    filter(variable %in% vars) %>%
    mutate(variable=factor(variable,levels=rev(vars)))
  
  p <- ggplot()
  
  for (t in c("MMQ-S","FMDI")) {
    if (t=="FMDI") p <- p + ggnewscale::new_scale_fill()
    
    cols <- if (t=="MMQ-S") c("#7851A9","#D9A400") else
      c("#D9A400","#7851A9")
    
    p <- p +
      geom_tile(data=filter(x,target==t),
                aes(target,variable,fill=r),
                width=.95,height=.95,colour="white",linewidth=.7) +
      geom_text(data=filter(x,target==t),
                aes(target,variable,label=label),
                size=3.5,lineheight=.85) +
      scale_fill_gradient2(
        low=cols[1],mid="white",high=cols[2],
        midpoint=0,limits=c(-1,1),
        name=paste0("<b>",t,"</b><br>Spearman's ρ [95% CI]"),
        guide=guide_colourbar(
          order=ifelse(t=="MMQ-S",1,2),
          barheight=grid::unit(3.2,"cm"),
          barwidth=grid::unit(.45,"cm"))
      )
  }
  
  p +
    scale_x_discrete(
      limits=c("MMQ-S","FMDI"),drop=FALSE,
      labels=function(z) paste0(
        c("MMQ-S"="Memory concern","FMDI"="Memory lapse\nfrequency")[z],
        "\n",z)
    ) +
    scale_y_discrete(
      labels=function(z) paste0(constructs[z],"\n",z)
    ) +
    coord_fixed(ratio=.38) +
    labs(title=title,x=NULL,y=NULL) +
    theme_minimal(base_size=10) +
    theme(
      panel.grid=element_blank(),
      plot.title = element_text(face="bold", size=13, hjust=0,
                                margin=margin(b=5, l=20)),
      plot.title.position = "plot",
      axis.text.x=element_text(size=10,lineheight=.95),
      axis.text.y=element_text(size=10,lineheight=.9),
      axis.ticks=element_blank(),
      legend.title=ggtext::element_markdown(size=10),
      legend.text=element_text(size=9.5)
    )
}

# Combine sections
p <- wrap_plots(
  Map(make_section,names(sections),sections),
  ncol=1,heights=lengths(sections),guides="collect"
) & theme(legend.position="right",
          legend.box="vertical",
          legend.title=ggtext::element_markdown(size=10))

print(p)

# Optional: inspect correlation results
print(cors,n=Inf,width=Inf)

# Save file
 ggsave(
  "correlation_heatmap_Nov02.tiff",
  p,
  width = 8,
  height = 13.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ---- Structural MRI data ----
## ---- clinical data ----
d <- data1 %>%
  mutate(
    group = factor(group_assignment_1, levels=c("2","1")),
    demo_2_sex_3 = factor(demo_2_sex_3),
    WMH_binary = case_when(
      WMH_normative == "Abnormal" ~ 1,
      WMH_normative == "Normal" ~ 0,
      TRUE ~ NA_real_
    )
  )

outcomes <- c("WMH_binary","incidentals_dummy")
covars <- c("demo_1_age_3","demo_2_sex_3","neuro_yrs_education_2")

analyze <- function(v) {
  x <- d %>% select(group,all_of(c(v,covars))) %>% tidyr::drop_na()
  f <- reformulate(c("group",covars),v)
  f0 <- reformulate(covars,v)
  
  m <- glm(f,data=x,family=binomial())
  m0 <- glm(f0,data=x,family=binomial())
  
  lr <- 2*as.numeric(logLik(m)-logLik(m0))
  b <- coef(m)["group1"]
  se <- sqrt(vcov(m)["group1","group1"])
  
  tibble(
    outcome=v,
    n_model=nrow(x),
    n_missing=nrow(d)-nrow(x),
    LR_X2=round(lr,3),
    p=signif(pchisq(lr,1,lower.tail=FALSE),4),
    OR=round(exp(b),3),
    CI_lower=round(exp(b-1.96*se),3),
    CI_upper=round(exp(b+1.96*se),3)
  )
}

results <- map_dfr(outcomes,analyze)
print(results)

# Counts by group, including missing
for(v in outcomes) {
  cat("\n",v,"\n")
  print(table(d$group,d[[v]],useNA="always"))
}

## ---- metrics ---- 
regions <- c("ROI3","ROI4","ROI5","Whole_brain")
metrics <- c("FA","MD","RD","AD","NDI","ODI")
vars <- unlist(lapply(metrics, \(m) paste0(regions,"_",m)))
covars <- c("demo_1_age_3","neuro_yrs_education_2",
            "demo_2_sex_3","rms_rel_avg_mm")

d <- data1 %>%
  mutate(group=factor(group_assignment_1,levels=c("2","1")),
         demo_2_sex_3=factor(demo_2_sex_3),
         across(all_of(grep("_(MD|RD|AD)$",vars,value=TRUE)),~.x*1000))

fmt <- function(x) sprintf("%.3f (%.3f)",mean(x),sd(x))

analyze <- function(v) {
  x <- d %>% select(group,all_of(c(v,covars))) %>%
    drop_na() %>% mutate(y=.data[[v]])
  
  f <- reformulate(c("group",covars),"y")
  f0 <- reformulate(covars,"y")
  fam <- list(gaussian=gaussian())
  if(all(x$y>0))
    fam <- c(fam,list(gamma=Gamma(link="log"),
                      inverse_gaussian=inverse.gaussian(link="log")))
  
  fit <- function(fm,ff) tryCatch(
    suppressWarnings(glm(fm,data=x,family=ff)),
    error=function(e) NULL)
  
  fits <- map(fam,~fit(f,.x)) %>%
    keep(~!is.null(.x) && isTRUE(.x$converged) && is.finite(AIC(.x)))
  if(!length(fits)) stop(paste("No valid model:",v))
  
  best <- names(which.min(map_dbl(fits,AIC)))
  m <- fits[[best]]
  m0 <- fit(f0,fam[[best]])
  if(is.null(m0) || !isTRUE(m0$converged))
    stop(paste("Reduced model failed:",v))
  
  lr <- max(0,2*as.numeric(logLik(m)-logLik(m0)))
  df <- attr(logLik(m),"df")-attr(logLik(m0),"df")
  
  x0 <- x1 <- x
  x0$group <- factor("2",levels=c("2","1"))
  x1$group <- factor("1",levels=c("2","1"))
  
  b <- coef(m); V <- vcov(m)
  X0 <- model.matrix(f,x0)[,names(b),drop=FALSE]
  X1 <- model.matrix(f,x1)[,names(b),drop=FALSE]
  u0 <- drop(X0%*%b); u1 <- drop(X1%*%b)
  
  if(best=="gaussian") {
    effect <- mean(u1-u0)
    grad <- colMeans(X1-X0)
  } else {
    u0 <- exp(u0); u1 <- exp(u1)
    effect <- mean(u1-u0)
    grad <- colMeans(sweep(X1,1,u1,"*")-sweep(X0,1,u0,"*"))
  }
  
  se <- sqrt(drop(t(grad)%*%V%*%grad))
  a <- x$y[x$group=="1"]; b0 <- x$y[x$group=="2"]
  sp <- sqrt(((length(a)-1)*var(a)+(length(b0)-1)*var(b0))/
               (length(a)+length(b0)-2))
  
  tibble(variable=v,family=best,FCD=fmt(a),Control=fmt(b0),
         LR_X2=lr,df=df,p=pchisq(lr,df,lower.tail=FALSE),
         SMD=effect/sp,CI_lower=(effect-1.96*se)/sp,
         CI_upper=(effect+1.96*se)/sp,
         n_model=nrow(x),n_missing=nrow(d)-nrow(x))
}

results <- map_dfr(vars,analyze) %>%
  mutate(holm_p=p.adjust(p,"holm"))

table <- results %>%
  mutate(across(c(LR_X2,SMD,CI_lower,CI_upper),~sprintf("%.3f",.x)),
         across(c(p,holm_p),~sprintf("%.12f",.x)))

print(table,n=Inf,width=Inf)
write.csv(table,"DTI_NODDI_results.csv",row.names=FALSE)

#head motion descriptives
d %>% group_by(group) %>%
  summarise(n=sum(!is.na(rms_rel_avg_mm)),
            mean_sd=fmt(na.omit(rms_rel_avg_mm)),.groups="drop") %>%
  print()

# ---- Neuropsychiatric interview findings ----
#limit to FCD group
d <- data1 %>%
  filter(as.character(group_assignment_1) == "1")

neuropsychiatric <- c(
  "dx_functioning_8", "dx_clin_response_8", "dx_clin_memofailure_8",
  "dx_clin_lapse_8", "dx_clin_neuro_8", "dx_clin_onset_8",
  "dx_clin_trauma_8", "dx_clin_decline_8", "dx_clin_fluc_8",
  "dx_clin_cog_8", "dx_clin_stress_8", "ds_clin_fam_8",
  "dx_clin_memory_8", "dx_somatic_8", "dx_somatic_comorbid_8",
  "dx_mdd_8", "dx_gad_8", "dx_panic_8", "dx_ptsd_8", "dx_sud_8",
  "dx_patient_rx_8", "dx_clin_appoint_8", "dx_clin_answer_8",
  "dx_clin_alexi_8", "dx_clin_disso_8", "consensus_1_8",
  "consensus_1a_8", "consensus_1b_8", "consensus_1c_8",
  "consensus_1d_8", "consensus_2_8", "consensus_3_8",
  "consensus_4_8", "consensus_5_8"
)

freq <- d %>%
  select(all_of(neuropsychiatric)) %>%
  pivot_longer(everything(), names_to="variable", values_to="response") %>%
  count(variable, response, name="n") %>%
  group_by(variable) %>%
  mutate(
    pct=round(100*n/sum(n),1),
    response=replace_na(as.character(response),"Missing/NA")
  ) %>%
  ungroup()

print(freq, n=Inf)

#psych dx
psych <- d %>%
  transmute(
    MDD = case_when(is.na(dx_mdd_8) ~ NA, TRUE ~ dx_mdd_8 %in% c(1,2)),
    GAD = dx_gad_8 == 1,
    Panic = dx_panic_8 == 1,
    PTSD = dx_ptsd_8 == 1,
    SUD = dx_sud_8 == 1,
    Somatic = dx_somatic_8 == 1
  )

psych_min <- rowSums(psych, na.rm=TRUE)
psych_complete <- complete.cases(psych)
psych_count <- ifelse(psych_complete, psych_min, NA_integer_)

cat("\nTotal FCD participants:", nrow(d),
    "\nAt least 1 confirmed psychiatric diagnosis:", sum(psych_min >= 1),
    "\nAt least 2 confirmed psychiatric diagnoses:", sum(psych_min >= 2),
    "\nConfirmed zero psychiatric diagnoses:", sum(psych_complete & psych_min == 0),
    "\nParticipants with incomplete psychiatric data:", sum(!psych_complete),
    "\nDistribution of exact counts (NA = incomplete):\n")

print(table(psych_count, useNA="always"))

#internal inconsistencies
inconsistency_vars <- c(
  "consensus_1a_8", "consensus_1b_8", "consensus_1c_8",
  "consensus_1d_8", "consensus_2_8", "consensus_3_8"
)

inc <- d[inconsistency_vars] == 1

inc_min <- rowSums(inc, na.rm=TRUE)
inc_complete <- complete.cases(inc)
inc_count <- ifelse(inc_complete, inc_min, NA_integer_)

cat("\nAt least 1 confirmed internal inconsistency:", sum(inc_min >= 1),
    "\nConfirmed zero internal inconsistencies:", sum(inc_complete & inc_min == 0),
    "\nParticipants with incomplete inconsistency data:", sum(!inc_complete),
    "\nDistribution of exact counts (NA = incomplete):\n")

print(table(inc_count, useNA="always"))

