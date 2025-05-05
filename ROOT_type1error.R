library(gamlss)

# refer to analyzeROOT function
source("analyze_ROOT.R")

#################################################################################################
# simulate type-I error under various structural zeros and shape parameters without covariates  #
#################################################################################################

set.seed(94769)
sim_type1_nocovar <- function(nsim, pi0_vals, n_per_group, a, b) {
  results <- expand.grid(sim = 1:nsim, pi0 = pi0_vals)
  results$p_value <- NA_real_
  results$converged <- FALSE
  
  for (i in seq_len(nrow(results))) {
    pi0 <- results$pi0[i]
    
    # Simulate zero-inflated data for two groups under no difference
    is_zero1 <- rbinom(n_per_group, 1, pi0)
    is_zero2 <- rbinom(n_per_group, 1, pi0)
    group1 <- ifelse(is_zero1 == 1, 0, rbeta(n_per_group, a, b))
    group2 <- ifelse(is_zero2 == 1, 0, rbeta(n_per_group, a, b))
    
    dat <- data.frame(
      outcome = c(group1, group2),
      studyarm = factor(rep(c("trt", "pbo"), each = n_per_group))
    )
    
    # Run analyzeROOT on simulated data under null hypothesis
    res <- tryCatch(
      analyzeROOT(dat, outcome = "outcome", treatment_var = "studyarm", baseline_vars = NULL),
      error = function(e) NULL
    )
    
    if (!is.null(res)) {
      results$converged[i] <- TRUE
      results$p_value[i] <- res$pairwise_contrasts$P[1]
    }
  }
  
  # Summarize Type I error and convergence rate
  summary_table <- do.call(rbind, lapply(pi0_vals, function(p) {
    rows <- results[results$pi0 == p, ]
    valid <- rows[rows$converged, ]
    data.frame(
      pi0 = p,
      Type1_Error = round(mean(valid$p_value < 0.05, na.rm = TRUE), 3),
      Convergence_Rate = round(mean(rows$converged), 3),
      n_valid = sum(rows$converged)
    )
  }))
  
  rownames(summary_table) <- NULL
  print(summary_table)
}

# N=50 per group 
sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)

# N=100 per group 
sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)

# N=200 per group 
sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)

### Results
# > # N=50 per group 
#   > sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.060                1    5000
# 2 0.2       0.049                1    5000
# 3 0.4       0.053                1    5000
# 4 0.6       0.056                1    5000
# 5 0.8       0.051                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.062                1    5000
# 2 0.2       0.058                1    5000
# 3 0.4       0.055                1    5000
# 4 0.6       0.053                1    5000
# 5 0.8       0.051                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.052                1    5000
# 2 0.2       0.058                1    5000
# 3 0.4       0.057                1    5000
# 4 0.6       0.051                1    5000
# 5 0.8       0.051                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.050                1    5000
# 2 0.2       0.061                1    5000
# 3 0.4       0.054                1    5000
# 4 0.6       0.053                1    5000
# 5 0.8       0.055                1    5000
# > 
#   > # N=100 per group 
#   > sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.051                1    5000
# 2 0.2       0.055                1    5000
# 3 0.4       0.047                1    5000
# 4 0.6       0.052                1    5000
# 5 0.8       0.049                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.048                1    5000
# 2 0.2       0.048                1    5000
# 3 0.4       0.050                1    5000
# 4 0.6       0.050                1    5000
# 5 0.8       0.045                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.051                1    5000
# 2 0.2       0.050                1    5000
# 3 0.4       0.053                1    5000
# 4 0.6       0.053                1    5000
# 5 0.8       0.052                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.057                1    5000
# 2 0.2       0.057                1    5000
# 3 0.4       0.054                1    5000
# 4 0.6       0.056                1    5000
# 5 0.8       0.056                1    5000
# > 
#   > # N=200 per group 
#   > sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.052                1    5000
# 2 0.2       0.046                1    5000
# 3 0.4       0.053                1    5000
# 4 0.6       0.049                1    5000
# 5 0.8       0.052                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.052                1    5000
# 2 0.2       0.051                1    5000
# 3 0.4       0.047                1    5000
# 4 0.6       0.047                1    5000
# 5 0.8       0.051                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.052                1    5000
# 2 0.2       0.057                1    5000
# 3 0.4       0.052                1    5000
# 4 0.6       0.052                1    5000
# 5 0.8       0.051                1    5000
# > sim_type1_nocovar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.050                1    5000
# 2 0.2       0.050                1    5000
# 3 0.4       0.052                1    4999
# 4 0.6       0.056                1    5000
# 5 0.8       0.050                1    5000

#############################################################################################
# simulate type-I error under various structural zeros and shape parameters with covariate  #
# baseline pain is predictive of outcome                                                    #
#############################################################################################

set.seed(73865)
sim_type1_covar <- function(nsim, pi0_vals, n_per_group, a, b) {
  results <- expand.grid(sim = 1:nsim, pi0 = pi0_vals)
  results$p_value <- NA_real_
  results$converged <- FALSE
  
  for (i in seq_len(nrow(results))) {
    pi0 <- results$pi0[i]
    
    baselinepain <- round(runif(n_per_group * 2, 4, 10))
    studyarm <- factor(rep(c("trt", "pbo"), each = n_per_group))
    
    # Structural zero model parameters
    beta_pi <- 0.4
    
    if (pi0 == 0) {
      prob_zero <- rep(0, length(baselinepain))
    } else {
      alpha_pi <- uniroot(
        function(alpha) mean(plogis(alpha + beta_pi * baselinepain)) - pi0,
        interval = c(-10, 10)
      )$root
      prob_zero <- plogis(alpha_pi + beta_pi * baselinepain)
    }
    
    # Simulate structural zeros
    is_zero <- rbinom(n_per_group * 2, 1, prob_zero)
    
    # Mean of beta model (no change)
    logit_mu <- 1.5 - 0.3 * baselinepain
    mu <- plogis(logit_mu)
    
    # Simulate outcome
    outcome <- ifelse(
      is_zero == 1,
      0,
      pmin(rbeta(n_per_group * 2, shape1 = mu * a * b, shape2 = (1 - mu) * a * b), 0.995)
    )
    dat <- data.frame(
      outcome = outcome,
      studyarm = studyarm,
      baselinepain = baselinepain
    )
    
    # Run analyzeROOT on simulated data under null hypothesis
    res <- tryCatch(
      analyzeROOT(dat, outcome = "outcome", treatment_var = "studyarm",  baseline_vars = "baselinepain"),
      error = function(e) NULL
    )
    
    if (!is.null(res)) {
      results$converged[i] <- TRUE
      results$p_value[i] <- res$pairwise_contrasts$P[1]
    }
  }
  
  # Summarize Type I error and convergence rate
  summary_table <- do.call(rbind, lapply(pi0_vals, function(p) {
    rows <- results[results$pi0 == p, ]
    valid <- rows[rows$converged, ]
    data.frame(
      pi0 = p,
      Type1_Error = round(mean(valid$p_value < 0.05, na.rm = TRUE), 3),
      Convergence_Rate = round(mean(rows$converged), 3),
      n_valid = sum(rows$converged)
    )
  }))
  
  rownames(summary_table) <- NULL
  print(summary_table)
}

# N=50 per group 
sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)

# N=100 per group 
sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)

# N=200 per group 
sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)

# > # N=50 per group 
#   > sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.058                1    5000
# 2 0.2       0.058                1    5000
# 3 0.4       0.060                1    5000
# 4 0.6       0.057                1    5000
# 5 0.8       0.047                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.056                1    5000
# 2 0.2       0.054                1    5000
# 3 0.4       0.058                1    5000
# 4 0.6       0.055                1    5000
# 5 0.8       0.044                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.057                1    5000
# 2 0.2       0.053                1    5000
# 3 0.4       0.053                1    5000
# 4 0.6       0.056                1    5000
# 5 0.8       0.046                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=50, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.064            0.895    4476
# 2 0.2       0.051            1.000    5000
# 3 0.4       0.053            1.000    5000
# 4 0.6       0.052            1.000    5000
# 5 0.8       0.033            1.000    5000
# > 
#   > # N=100 per group 
#   > sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.054                1    5000
# 2 0.2       0.049                1    5000
# 3 0.4       0.056                1    5000
# 4 0.6       0.053                1    5000
# 5 0.8       0.047                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.057                1    5000
# 2 0.2       0.052                1    5000
# 3 0.4       0.052                1    5000
# 4 0.6       0.050                1    5000
# 5 0.8       0.050                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.053                1    5000
# 2 0.2       0.056                1    5000
# 3 0.4       0.053                1    5000
# 4 0.6       0.050                1    5000
# 5 0.8       0.040                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=100, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.061            0.741    3707
# 2 0.2       0.061            1.000    5000
# 3 0.4       0.054            1.000    5000
# 4 0.6       0.048            1.000    5000
# 5 0.8       0.046            1.000    5000
# > 
#   > # N=200 per group 
#   > sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.051                1    5000
# 2 0.2       0.050                1    5000
# 3 0.4       0.049                1    5000
# 4 0.6       0.055                1    5000
# 5 0.8       0.050                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 1, b = 3)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.055                1    5000
# 2 0.2       0.051                1    5000
# 3 0.4       0.056                1    5000
# 4 0.6       0.054                1    5000
# 5 0.8       0.045                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 2, b = 2)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.055                1    5000
# 2 0.2       0.054                1    5000
# 3 0.4       0.053                1    5000
# 4 0.6       0.054                1    5000
# 5 0.8       0.051                1    5000
# > sim_type1_covar(nsim = 5000, n_per_group=200, pi0_vals = seq(0, 0.8, 0.2), a = 0.5, b = 0.5)
# pi0 Type1_Error Convergence_Rate n_valid
# 1 0.0       0.056            0.672    3362
# 2 0.2       0.060            1.000    5000
# 3 0.4       0.056            1.000    5000
# 4 0.6       0.050            1.000    5000
# 5 0.8       0.047            1.000    5000

