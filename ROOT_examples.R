library(gamlss)

source("analyze_ROOT.R")

# Create example dataset
set.seed(937)
n <- 100
trt_grps <- c("placebo", "treatment1", "treatment2")

# baseline pain
baseline_pain_placebo <- round(runif(n, 4, 9))
baseline_pain_treat1  <- round(runif(n, 4, 9))
baseline_pain_treat2  <- round(runif(n, 4, 9))

# age
age_placebo <- rnorm(100, 45, 8)
age_treat1  <- rnorm(100, 45, 8)
age_treat2  <- rnorm(100, 45, 8)

# non-responders (structural zeros)
zero_placebo <- rbinom(n, 1, ifelse(baseline_pain_placebo > 7, 0.5, 0.2))
zero_treat1  <- rbinom(n, 1, ifelse(baseline_pain_treat1 > 7, 0.35, 0.15))
zero_treat2  <- rbinom(n, 1, ifelse(baseline_pain_treat2 > 7, 0.35, 0.15))

# ROOT
placebo <- ifelse(zero_placebo == 1, 0,
                  ifelse(baseline_pain_placebo > 7, rbeta(n, 0.8, 2), rbeta(n, 1, 2)))
treat1 <- ifelse(zero_treat1 == 1, 0,
                 ifelse(baseline_pain_treat1 > 7, rbeta(n, 2, 2), rbeta(n, 3, 2)))
treat2 <- ifelse(zero_treat2 == 1, 0,
                 ifelse(baseline_pain_treat2 > 7, rbeta(n, 2, 2), rbeta(n, 3, 2)))

# 3-arm example
example_3arms <- data.frame(
  ROOT = c(placebo, treat1, treat2),
  treatment = factor(rep(trt_grps, each = n)),
  blpain = c(baseline_pain_placebo, baseline_pain_treat1, baseline_pain_treat2),
  blage = c(age_placebo, age_treat1, age_treat2)
)

# 2-arm example
example_2arms <- data.frame(
  ROOT = c(placebo, treat1),
  treatment = factor(rep(trt_grps[1:2], each = n)),
  blpain = c(baseline_pain_placebo, baseline_pain_treat1),
  blage = c(age_placebo, age_treat1)
)

# Observed Means
mean(placebo)
mean(treat1)
mean(treat2)

## Two-arm comparison

# Analyze ROOT using zero-inflated beta regression model, no covariate
analysis_2arms <- analyzeROOT(
  data = example_2arms,
  outcome = "ROOT",
  treatment_var = "treatment",
  baseline_vars = NULL
)

analysis_2arms$marginal_means
analysis_2arms$pairwise_contrasts

# Analyze ROOT using zero-inflated beta regression model, with baseline pain
analysis_2arms <- analyzeROOT(
  data = example_2arms,
  outcome = "ROOT",
  treatment_var = "treatment",
  baseline_vars = "blpain"
)

analysis_2arms$marginal_means
analysis_2arms$pairwise_contrasts


## Three-arm comparison

# Analyze ROOT using zero-inflated beta regression model, with both covariates
analysis_3arms <- analyzeROOT(
  data = example_3arms,
  outcome = "ROOT",
  treatment_var = "treatment",
  baseline_vars = c("blpain", "blage")
)

analysis_3arms$marginal_means
analysis_3arms$pairwise_contrasts

