library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(lubridate)
library(WeightIt)
library(balanceHD)
library(ebal)
library(cobalt)
library(R.utils)
source("utilities.R")

argv <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

EST_OUTCOME <- as.logical(get_attr_default(argv, "est_outcome", TRUE))
OUTCOME_CV <- as.logical(get_attr_default(argv, "outcome_cv", TRUE))
Y_LAMBDA <- as.numeric(get_attr_default(argv, "y_lambda", 1))
Y_ALPHA <- as.numeric(get_attr_default(argv, "y_alpha", 1))

EST_PROPENSITY <- as.logical(get_attr_default(argv, "est_propensity", TRUE))
PROP_CV <- as.logical(get_attr_default(argv, "prop_cv", TRUE))
T_LAMBDA <- as.numeric(get_attr_default(argv, "t_lambda", 1))
T_ALPHA <- as.numeric(get_attr_default(argv, "t_alpha", 1))

## Dot product between outcome and propensity coefficients
AB_DP  <- as.numeric(get_attr_default(argv, "ab_dp", 0.75))

estimand <- as.character(get_attr_default(argv, "estimand", "ATT"))

tau <- as.numeric(get_attr_default(argv, "tau", 0))
coef_setting <- as.numeric(get_attr_default(argv, "coef", 1))
mscale <- as.numeric(get_attr_default(argv, "mscale", 2))
escale <- as.numeric(get_attr_default(argv, "escale", 4))

eta_clip <- as.numeric(get_attr_default(argv, "eta", 0.1))

times <- as.numeric(get_attr_default(argv, "times", 50))
iters <- as.numeric(get_attr_default(argv, "iters", 50))


sigma2_y <- as.numeric(get_attr_default(argv, "sigma2_y", 1))

n <- as.numeric(get_attr_default(argv, "n", 500))
p <- as.numeric(get_attr_default(argv, "p", 1000))

use_vectorized <- as.logical(get_attr_default(argv, "vec", TRUE))


print(sprintf("Using mscale: %s", mscale))
print(sprintf("Using escale: %s", escale))
print(sprintf("Sample size (n): %s", n))
print(sprintf("Dimension (p): %s", p))
print(sprintf("Running for %s iters", iters))
print(sprintf("Drawing %s null space vectors", times))
print(sprintf("Estimating propensity? %s", EST_PROPENSITY))

eigen_debug <- FALSE
bias_debug <- FALSE
bias_times <- if(bias_debug) 1 else times

w2_scale_vec <- c(seq(1, -1, by=-.2))

results_array <- array(dim=c(iters, 11, length(w2_scale_vec)))
w2lim_true_vec <- numeric(iters)

eta_matrix <- matrix(NA, nrow=iters, ncol=length(w2_scale_vec))

dimnames(results_array) <- list(1:iters,
                                c("Regression", "IPW", "IPW_d", "IPW_d_oracle",
                                  "IPW_clip", "AIPW", "AIPW_d", "AIPW_clip",
                                  "Naive", "BalanceHD", "BalanceHD_weights_only"),
                                w2_scale_vec)

## Generate iters datasets
for(iter  in 1:iters) {

  ## ################
  ## Generate dataset
  ## #################

  ## coef_setting controls the number of non-zeros in the outcome model
  ## alpha are the coefficients of the outcome model
  coef_settings <- c(20, 50, 80)
  qa <- coef_settings[coef_setting]
  if(p < qa)
    qa <- p

  alpha <- c(rnorm(qa) / qa, rep(0, p - qa))
  alpha <- alpha / sqrt(sum(alpha^2))


  ## beta are the coefficients of the propensity model
  ## a random (dense) vector with inner product of AB_DP with alpha
  NC  <- rstiefel::NullC(alpha)
  beta  <- AB_DP*alpha + sqrt(1-AB_DP^2)*(NC  %*% rustiefel(p-1, 1))
  

  true_ate <- tau

######################
  ## Linear outcome model and logist treatment model
  ## Function from utilities.R
######################
  simdat <- gen_linearY_logisticT(n, p, tau, alpha, beta, mscale, escale, sigma2_y)

  X <- simdat$X
  T <- simdat$T
  Y <- simdat$Y
  m <- simdat$m
  e <- simdat$e

  ## #############################################
  ## Fit the outcome model using glmnet (see utilites.R)
  ## #############################################
  
  out_ests <- if(EST_OUTCOME) estimate_outcome(X, T, Y, estimand, alpha=Y_ALPHA, include_intercept=TRUE,
                                               coef_policy="lambda.min", pred_policy="lambda.1se")
              else list()

  alpha_hat <- get_attr_default(out_ests, "alpha_hat", alpha)
  outcome_intercept <- get_attr_default(out_ests, "intercept", 0)
  alpha_hat_normalized <- get_attr_default(out_ests, "alpha_hat_normalized",
                                           alpha / sqrt(sum(alpha^2)))
  mhat0 <- get_attr_default(out_ests, "mhat0", mean(X %*% alpha))
  mhat1 <- get_attr_default(out_ests, "mhat1", mean(Y[T==1]))
  tau_hat <- get_attr_default(out_ests, "tau_hat", tau)
  outcome_preds <- get_attr_default(out_ests, "preds", X %*% alpha)

  if(EST_PROPENSITY) {
    prop_ests <- estimate_propensity(X, T, cv=PROP_CV,
                                     T_lambda_min=T_LAMBDA,
                                     eta=eta_clip, alpha=T_ALPHA)
  } else {
    prop_ests <- list()
  }

  beta_hat <- get_attr_default(prop_ests, "beta_hat", beta * escale)
  beta_hat_normalized <- get_attr_default(prop_ests, "beta_hat_normalized", beta)
  escale_hat <- get_attr_default(prop_ests, "escale_hat", escale)
  ehat <- get_attr_default(prop_ests, "ehat", e)
  xb <- X %*% beta_hat_normalized


  ## etas for all glm fits
  ehat_all <- get_attr_default(prop_ests, "ehat_all", matrix(e, ncol=1))
  eta_hat <- pmin(1 - apply(ehat_all, 2, max), apply(ehat_all, 2, min))





  alpha_normalized <- alpha/sqrt(sum(alpha^2))
  beta_normalized <- beta/sqrt(sum(beta^2))

  ab_dot_prod_true <- as.numeric(t(alpha_normalized) %*% beta_normalized)
  mvecs_true <- cbind((alpha_normalized + beta_normalized) / sqrt(2 + 2 * ab_dot_prod_true),
  (alpha_normalized - beta_normalized) / sqrt(2 - 2 * ab_dot_prod_true))
  mvals_true <- c((ab_dot_prod_true + 1)/2, (ab_dot_prod_true - 1)/2)
  if(ab_dot_prod_true < 0){
    mvals_true <- mvals_true[2:1]
    mvecs_true <- mvecs_true[, 2:1]
  }


  ## ###################################
  ## Non-reduced or regularized estimands include regression estimator
  ## and simple difference in means
  ## ###############################
  results_array[iter, "Regression", ] <- tau_hat

  ## Naive difference in means
  naive <- mean(Y[T == 1]) - mean(Y[T == 0])
  results_array[iter, "Naive", ] <- naive

  ## If regularization completely collapses ehat,
  ## no further reduction is possible
  if(var(ehat) == 0) {

    results_array[iter, "IPW_d", ] <- ipw
    results_array[iter, "AIPW_d", ] <- aipw

    next
  }

  ## ####################################
  ## balanceHD (Wager and Athey)
  ## method.fit = "none" does weights only (no outcome model)
  ## target.pop=1 means ATT
  ## When not weights only, the inferred outcome model
  ## is identical to the one inferred by estimate_outcome
  ## Thus, the only difference between balanceHD and our reducer method
  ## Is the function that we balance on
  ## ####################################
  residual_balance <- residualBalance.ate(X, Y, T, target.pop = 1,
                                          alpha=Y_ALPHA)

  residual_balance_weights_only <- residualBalance.ate(X, Y, T, target.pop = 1,
                                                       alpha=Y_ALPHA, fit.method="none")

  ## ####################################
  ## Compute hyperbola for reduction
  ## ###################################

  ab_dot_prod <- as.numeric(t(alpha_hat_normalized) %*% beta_hat_normalized)
  mvecs <- cbind((alpha_hat_normalized + beta_hat_normalized)/sqrt(2 + 2 * ab_dot_prod),
  (alpha_hat_normalized - beta_hat_normalized)/sqrt(2 - 2 * ab_dot_prod))
  mvals <- c((ab_dot_prod + 1)/2, (ab_dot_prod - 1)/2)

  if(eigen_debug){
    ab_outer <- alpha_hat_normalized %*% t(beta_hat_normalized)
    spec <- eigen(0.5 * (ab_outer + t(ab_outer)))
    debug_mvals <- spec$values[c(1,p)]
    debug_mvecs <- spec$vectors[,c(1,p)]
  }

  if(is.na(ab_dot_prod))
    browser()

  ## Depending on correlation between alpha and beta, switch eigenspace labels to
  ## traverse the "closed" side of the hyperbola as we vary the w2 parameter
  if(ab_dot_prod < 0){
    mvals <- mvals[2:1]
    mvecs <- mvecs[,2:1]
  }

  ## ######################################################
  ## IPW_d_oracle
  ##
  ## Reduced-IPW estimator based on the true propensity score and true outcome model
  ## Reductions based on the true propensity score and outcome model
  ## ######################################################
  w2_lim_true <- -sqrt(abs(mvals_true[2]))
  w2lim_true_vec[iters] <- w2_lim_true
  for(j in 1:length(w2_scale_vec)) {

    w2scale <- w2_scale_vec[j]
    bias <- get_bias(T=T, Y=Y, X=X, xb=X %*% beta_normalized,
                     estimand=estimand,
                     mvecs=mvecs_true, mvals=mvals_true,
                     ab_dot_prod=ab_dot_prod_true,
                     escale=escale,
                     w2=w2scale*w2_lim_true,
                     w2lim=w2_lim_true,
                     times=bias_times,
                     DEBUG=bias_debug,
                     alpha_hat_normalized=alpha_normalized,
                     beta_hat_normalized=beta_normalized)

    ipw_d_oracle <- bias$bias1 - bias$bias0
    results_array[iter, "IPW_d_oracle", j] <- ipw_d_oracle
  }


  ## ######################################################
  ## Compute reduced estimands for grid of w2 value in [-1, 1]
  ## when j=1, gamma = alpha_hat
  ## when j=length(w2_scale_vec), gamma = beta_hat
  ## ######################################################
  for(j in 1:length(w2_scale_vec)) {
    w2scale <- w2_scale_vec[j]

    print(paste(w2scale, iter, sep=", "))

    ## w2 limits are determined by the eigenvalue corresponding to the "closed"
    ## side of the hyperbola.
    ## at w2 = w2_lim, d(X) = e(X); at w2 = -w2_lim, d(X) = mhat(X)
    w2_lim <- -sqrt(abs(mvals[2]))
    w2 <- w2scale * w2_lim

    if(estimand == "ATT") {
      residual <- 0 + (1-T) * (Y - outcome_preds)
    }
    else {
      warning("ATE not currently supported")
      residual <- Y - X %*% alpha_hat - T * tau_hat - outcome_intercept
    }

    ## ############# IPW_d ###################
    ## Compute IPW_d, by reweighting using reduction d(X)
    ## 'w2' is the tuning parameter for how similar the deconfounding score
    ## is to the propensity and prognostic score (w2=0 is equidistant)
    ## 'times' is the number of Monte Carlo samples to use to compute weighted estimates
    ## Larger values reduces Monte Carlo variance but increases computation time
    ## ######################################
    bias <- get_bias(T=T, Y=Y, X=X, xb=xb, estimand=estimand,
                     mvecs=mvecs, mvals=mvals,
                     ab_dot_prod=ab_dot_prod, escale=escale_hat,
                     w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug,
                     alpha_hat_normalized=alpha_hat_normalized, beta_hat_normalized=beta_hat_normalized)

    ipw_d <- bias$bias1 - bias$bias0

    ## the most extreme reduced propensity score value
    ## used to compute the comparable clipping estimator
    eta <- bias$eta

    ## ####################################
    ## Standard IPW and AIPW estimators
    ## As a fair comparison with our reducer estimator,
    ## use the glmfit based on the regularization setting
    ## for which the most extreme propensity score matches
    ## the most extreme reduced score e_d (See Figure 3)
    ## ##################################
    eta_matrix[iter, j] <- eta
    ehat_match <- ehat_all[, which.min(abs(eta_hat - eta))]

    ## Standard IPW
    ipw <- ipw_est(ehat_match, T, Y, estimand, hajek=TRUE)

    ## Standard AIPW
    aipw <- tau_hat +
      ipw_est(ehat_match, T, residual, estimand, hajek=TRUE)
    if(is.na(aipw))
      browser()

    ## ###############################
    ## Compute IPW and AIPW based on clipped propensity scores
    ## ###############################

    ehat_clip <- ehat
    ehat_clip[ehat < eta] <- eta
    ehat_clip[ehat > 1 - eta] <- 1 - eta

    ipw_clip <- ipw_est(ehat_clip, T, Y, estimand, hajek=TRUE)
    aipw_clip <- tau_hat +
      ipw_est(ehat_clip, T, residual, estimand, hajek=TRUE)



    ## ############# AIPW_d ######################################
    ## Compute negative regression bias by balancing on residuals
    ## 'w2' is the tuning parameter for how similar the deconfounding score
    ## is to the propensity and prognostic score (w2=0 is equidistant)
    ## 'times' is the number of Monte Carlo samples to use to compute weighted estimates
    ## Larger values reduces Monte Carlo variance but increases computation time
    ## ################################################
    bias <- get_bias(T=T, Y=residual, X=X, xb=xb, estimand=estimand,
                     mvecs=mvecs, mvals=mvals,
                     ab_dot_prod=ab_dot_prod, escale=escale_hat,
                     w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug,
                     alpha_hat_normalized=alpha_hat_normalized, beta_hat_normalized=beta_hat_normalized)

    ## Correct bias and compute AIPW_d
    aipw_d <- tau_hat + bias$bias1 - bias$bias0

    ## ###############################################
    
    results_array[iter, "IPW_d", j] <- ipw_d
    results_array[iter, "AIPW_d", j] <- aipw_d
    results_array[iter, "IPW", j] <- ipw
    results_array[iter, "IPW_clip", j] <- ipw_clip
    results_array[iter, "AIPW", j] <- aipw
    results_array[iter, "AIPW_clip", j] <- aipw_clip
    results_array[iter, "BalanceHD", j] <- residual_balance
    results_array[iter, "BalanceHD_weights_only", j] <- residual_balance_weights_only

    print(sprintf("Naive: %.3f, Reg: %.3f, IPW: %.3f, IPW-clip: %.3f, BalanceHD: %.3f, BalanceHD_weights_only: %.3f, IPW-d: %.3f, IPW-d-oracle: %.3f, AIPW: %.3f, AIPW-clip: %.3f, AIPW-d: %.3f",
                  naive, tau_hat, ipw, ipw_clip, residual_balance, residual_balance_weights_only,
                  ipw_d, results_array[iter, "IPW_d_oracle", j],
                  aipw, aipw_clip, aipw_d))
  }

}

## RMSE
sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

## MAD
apply(abs(results_array - true_ate), c(2, 3), function(x) median(x, na.rm=TRUE))

save(results_array, true_ate, w2lim_true_vec, eta_matrix, 
     file=sprintf("results/results_n%i_p%i_coef%.2f_escale%.1f_mscale%.1f_yalpha%i_talpha%i_estpropensity%s_%s_%s.RData",
                  n, p, coef_setting, escale, mscale, Y_ALPHA, T_ALPHA, EST_PROPENSITY, estimand, ab_dot_prod_true,
                  gsub(" ", "", now(), fixed=TRUE)))

## RMSE, Bias and Variance as in Figure 2 of paper
make_bias_var_plot(results_array, true_ate)
