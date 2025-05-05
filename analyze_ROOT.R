library(gamlss)

analyzeROOT <- function(data, outcome, treatment_var, baseline_vars) {
  
  # Fit zero-inflated beta regression model
  vars <- c(treatment_var, baseline_vars)
  formula_mu <- reformulate(vars, response = outcome)
  formula_pi <- reformulate(vars)
  
  fit <- gamlss(
    formula = formula_mu,
    sigma.formula = ~ 1,   
    nu.formula = formula_pi,
    family = BEZI,
    data = data,
    control = gamlss.control(trace = FALSE)
  )
  
  # Set covariates for estimation of adjusted marginal means
  trt_levels <- levels(data[[treatment_var]])
  
  if (is.null(baseline_vars)) {
    newdata <- data.frame(treatment = factor(trt_levels, levels = trt_levels))
    colnames(newdata) <- treatment_var
  } else {
    baseline_means <- sapply(data[baseline_vars], mean, na.rm = TRUE)
    newdata <- as.data.frame(matrix(rep(baseline_means, times = length(trt_levels)),
                                    nrow = length(trt_levels), byrow = TRUE))
    colnames(newdata) <- baseline_vars
    newdata[[treatment_var]] <- factor(trt_levels, levels = trt_levels)
    newdata <- newdata[, c(treatment_var, baseline_vars)]
  }
  
  
  # Design matrices
  X_mu <- model.matrix(reformulate(vars), data = newdata)
  X_pi <- model.matrix(reformulate(vars), data = newdata)
  
  # Extract and relabel coefficients 
  # (gamlss refers to logistic submodel coefficients as 'nu', rather than 'pi')
  beta_mu <- coef(fit, what = "mu")
  beta_pi <- coef(fit, what = "nu")
  names(beta_mu) <- paste0("mu_", names(beta_mu))
  names(beta_pi) <- paste0("pi_", names(beta_pi))
  beta_all <- c(beta_mu, beta_pi)
  
  # Full variance-covariance matrix (drop precision parameter)
  V_all <- vcov(fit)
  mu_names <- paste0("mu_", names(coef(fit, "mu")))
  sigma_names <- paste0("sigma_", names(coef(fit, "sigma")))
  pi_names <- paste0("pi_", names(coef(fit, "nu")))
  rownames(V_all) <- colnames(V_all) <- c(mu_names, sigma_names, pi_names)
  keep_names <- c(mu_names, pi_names)
  V <- V_all[keep_names, keep_names]
  
  # Linear predictors
  mu_hat <- plogis(as.numeric(X_mu %*% beta_mu))
  pi_hat <- plogis(as.numeric(X_pi %*% beta_pi))
  
  # Adjusted marginal means
  f <- mu_hat * (1 - pi_hat)
  
  # Gradients
  grad_mu <- (1 - pi_hat) * mu_hat * (1 - mu_hat) * X_mu
  grad_pi <- -mu_hat * pi_hat * (1 - pi_hat) * X_pi
  grad_list <- lapply(seq_along(mu_hat), function(i) c(grad_mu[i, ], grad_pi[i, ]))
  
  # SEs by delta method and CIs
  f_se <- sapply(grad_list, function(g) sqrt(t(g) %*% V %*% g))
  f_lower <- f - qnorm(0.975) * f_se
  f_upper <- f + qnorm(0.975) * f_se
  
  # Summary of adjusted marginal mean statistics
  marginal_df <- data.frame(
    TrtGrp = trt_levels,
    Estimate = round(f, 3),
    SE = round(f_se, 3),
    Lower95CI = round(f_lower, 3),
    Upper95CI = round(f_upper, 3)
  )
  
  # Pairwise contrasts
  contrast_df <- do.call(rbind, combn(seq_along(f), 2, function(idx) {
    i <- idx[1]; j <- idx[2]
    delta <- f[i] - f[j]
    var_i <- t(grad_list[[i]]) %*% V %*% grad_list[[i]]
    var_j <- t(grad_list[[j]]) %*% V %*% grad_list[[j]]
    cov_ij <- t(grad_list[[i]]) %*% V %*% grad_list[[j]]
    delta_var <- var_i + var_j - 2 * cov_ij
    delta_se <- sqrt(delta_var)
    z <- delta / delta_se
    p <- 2 * (1 - pnorm(abs(z)))
    data.frame(
      Contrast = paste(trt_levels[i], "-", trt_levels[j]),
      Estimate = round(delta, 3),
      SE = round(delta_se, 3),
      Lower95CI = round(delta - qnorm(0.975) * delta_se, 3),
      Upper95CI = round(delta + qnorm(0.975) * delta_se, 3),
      P = signif(p, digits=4)
    )
  }, simplify = FALSE))
  
  return(list(
    model = fit,
    marginal_means = marginal_df,
    pairwise_contrasts = contrast_df)
  )
}