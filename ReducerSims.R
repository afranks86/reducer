library(mvtnorm)
library(rstiefel)
library(glmnet)
library(lubridate)
library(R.utils)
source("utilities.R")

argv <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

EST_OUTCOME <- as.logical(get_attr_default(argv, "est_outcome", TRUE))
OUTCOME_CV <- as.logical(get_attr_default(argv, "outcome_cv", TRUE))
Y_LAMBDA <- as.numeric(get_attr_default(argv, "y_lambda", 0))
Y_ALPHA <- as.numeric(get_attr_default(argv, "y_alpha", 0))

EST_PROPENSITY <- as.logical(get_attr_default(argv, "est_propensity", TRUE))
PROP_CV <- as.logical(get_attr_default(argv, "prop_cv", TRUE))
T_LAMBDA <- as.numeric(get_attr_default(argv, "t_lambda", 1))
T_ALPHA <- as.numeric(get_attr_default(argv, "t_alpha", 1))

estimand <- as.character(get_attr_default(argv, "estimand", "ATT"))

tau <- as.numeric(get_attr_default(argv, "tau", 1))
coef_setting <- as.numeric(get_attr_default(argv, "coef", 0))
mscale <- as.numeric(get_attr_default(argv, "mscale", 5))
escale <- as.numeric(get_attr_default(argv, "escale", 4))

eta_clip <- as.numeric(get_attr_default(argv, "eta", 0.1))

times <- as.numeric(get_attr_default(argv, "times", 50))
iters <- as.numeric(get_attr_default(argv, "iters", 50))

sigma2_y <- as.numeric(get_attr_default(argv, "sigma2_y", 4))

n <- as.numeric(get_attr_default(argv, "n", 100))
p <- as.numeric(get_attr_default(argv, "p", 100))

use_vectorized <- as.logical(get_attr_default(argv, "vec", TRUE))
get_bias <- if(use_vectorized) get_bias_vec else get_bias_old


print(sprintf("Using mscale: %s", mscale))
print(sprintf("Using escale: %s", escale))
print(sprintf("Sample size (n): %s", n))
print(sprintf("Dimension (p): %s", p))
print(sprintf("Running for %s iters", iters))
print(sprintf("Drawing %s null space vectors", times))
print(sprintf("Estimating propensity? %s", EST_PROPENSITY))

## w2_scale_vec <- c(seq(0.95, 0.5, by=-0.05))
## w2_scale_vec <- c(w2_scale_vec, 0, -rev(w2_scale_vec))
eigen_debug <- FALSE
bias_debug <- FALSE
bias_times <- if(bias_debug) 1 else times

w2_scale_vec <- c(seq(1, -1, by=-.2))

results_array <- array(dim=c(iters, 9, length(w2_scale_vec)))
w2lim_true_vec <- numeric(iters)

eta_matrix <- matrix(NA, nrow=iters, ncol=length(w2_scale_vec))

dimnames(results_array) <- list(1:iters,
                                c("Regression", "IPW", "IPW_d", "IPW_d_oracle",
                                  "IPW_clip", "AIPW", "AIPW_d", "AIPW_clip",
                                  "Naive"),
                                w2_scale_vec)

for(iter  in 1:iters) {

    ## ################
    ## Generate dataset
    ## #################
    
    if(coef_setting == 1){
        #alpha <-  1 / (1 + (23 * (0:(length(beta.main) - 1))) %% length(beta.main))
        alpha <- 1 / (1:p)
        alpha <- alpha / sqrt(sum(alpha^2))
        #alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        #beta <- c(1, 0, 1, rep(0, p-3))/sqrt(2)
        beta <- c(rep(1/40, 50), rep(0, p-50))
    } else {
        alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        beta <- c(-1, 0.75, 1, rep(0, p-3))/sqrt(3)
    }
    
    true_ate <- tau
    
    simdat <- gen_linearY_logisticT(n, p, tau, alpha, beta, mscale, escale, sigma2_y)
    
    X <- simdat$X
    T <- simdat$T
    Y <- simdat$Y
    m <- simdat$m
    e <- simdat$e
    
    ## ################
    ## Set nuisance parameters
    ## #################
    
    out_ests <- if(EST_OUTCOME) estimate_outcome(X, T, Y, estimand, cv=OUTCOME_CV, Y_lambda_min=Y_LAMBDA, alpha=Y_ALPHA)
                else list()
   
    alpha_hat <- get_attr_default(out_ests, "alpha_hat", alpha)
    alpha_hat_normalized <- get_attr_default(out_ests, "alpha_hat_normalized",
                                             alpha / sqrt(sum(alpha^2)))
    mhat0 <- get_attr_default(out_ests, "mhat0", X %*% alpha)
    mhat1 <- get_attr_default(out_ests, "mhat1", X %*% alpha + tau)
    tau_hat <- get_attr_default(out_ests, "tau_hat", tau)
    
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
    

    ## ##################
    ## Non-reduced estimands
    ## ##################
    
    ## IPW d using true prognostic score and true propensity score
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
    ## Compute hyperbola for reduction
    ## ##################
    
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
   
    ## Depending on correlation between alpha and beta, switch eigenspace labels to
    ## traverse the "closed" side of the hyperbola as we vary the w2 parameter
    if(ab_dot_prod < 0){
        mvals <- mvals[2:1]
        mvecs <- mvecs[,2:1]
    }

    for(j in 1:length(w2_scale_vec)) {
        w2scale <- w2_scale_vec[j]
        
        print(paste(w2scale, iter, sep=", "))

        # w2 limits are determined by the eigenvalue corresponding to the "closed"
        # side of the hyperbola.
        # at w2 = w2_lim, d(X) = e(X); at w2 = -w2_lim, d(X) = mhat(X)
        w2_lim <- -sqrt(abs(mvals[2]))
        w2 <- w2scale * w2_lim
        
        if(estimand == "ATT") {
            residual <- 0 + (1-T) * (Y - X %*% alpha_hat)
        }
        else {
            residual <- Y - X %*% alpha_hat - T * tau_hat
        }
        
        ## Compute IPW_d, using reduced d 
        ## w2 is the tuning parameter for how similar d is to to ehat vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger should reduce variance)
        bias <- get_bias(T=T, Y=Y, X=X, xb=xb, estimand=estimand,
                         mvecs=mvecs, mvals=mvals,
                         ab_dot_prod=ab_dot_prod, escale=escale_hat,
                         w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug,
                         alpha_hat_normalized=alpha_hat_normalized, beta_hat_normalized=beta_hat_normalized)

        ipw_d <- bias$bias1 - bias$bias0
        eta <- bias$eta
        cache_edd <- bias$e_dd
        
        ## Find glmnet fit that matches most extreme pscore
        eta_matrix[iter, j] <- eta
        ehat_match <- ehat_all[, which.min(abs(eta_hat - eta))]

        ## Standard IPW
        ipw <- ipw_est(ehat_match, T, Y, estimand, hajek=TRUE)
        
        ## Standard AIPW
        aipw <- tau_hat +
            ipw_est(ehat_match, T, residual, estimand, hajek=TRUE)

        ## Compute clipped estimators
        ehat_clip <- ehat
        ehat_clip[ehat < eta] <- eta
        ehat_clip[ehat > 1 - eta] <- 1 - eta

        ipw_clip <- ipw_est(ehat_clip, T, Y, estimand, hajek=TRUE)
        aipw_clip <- tau_hat +
            ipw_est(ehat_clip, T, residual, estimand, hajek=TRUE)

        ## Compute negative regression bias by balancing on residuals
        ## w2 is the tuning parameter for how similar to e vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger shoudl reduce variance)
        #bias <- get_bias(T=T, Y=residual, X=X, xb=xb, estimand=estimand,
        #                 mvecs=mvecs, mvals=mvals,
        #                 ab_dot_prod=ab_dot_prod, escale=escale_hat,
        #                 w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug,
        #                 alpha_hat_normalized=alpha_hat_normalized, beta_hat_normalized=beta_hat_normalized)
        bias <- compute_bias(T=T, Y=residual, e_dd=cache_edd, estimand=estimand)

        ## Correct bias and compute AIPW_d
        aipw_d <- tau_hat + bias$bias1 - bias$bias0

        results_array[iter, "IPW_d", j] <- ipw_d
        results_array[iter, "AIPW_d", j] <- aipw_d
        results_array[iter, "IPW", j] <- ipw
        results_array[iter, "IPW_clip", j] <- ipw_clip
        results_array[iter, "AIPW", j] <- aipw
        results_array[iter, "AIPW_clip", ] <- aipw_clip


        print(sprintf("Naive: %.3f, Reg: %.3f, IPW: %.3f, IPW-clip: %.3f, IPW-d: %.3f, IPW-d-oracle: %.3f, AIPW: %.3f, AIPW-clip: %.3f, AIPW-d: %.3f",
                      naive, tau_hat, ipw, ipw_clip,
                      ipw_d, results_array[iter, "IPW_d_oracle", j],
                      aipw, aipw_clip, aipw_d))
    }

}

sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

apply(abs(results_array - true_ate), c(2, 3), function(x) median(x, na.rm=TRUE))

save(results_array, true_ate, w2lim_true_vec, eta_matrix, 
     file=sprintf("results/results_n%i_p%i_coef%i_escale%.1f_mscale%.1f_yalpha%i_estpropensity%s_%s_%s.RData",
                  n, p, coef_setting, escale, mscale, Y_ALPHA, EST_PROPENSITY, estimand,
                  gsub(" ", "", now(), fixed=TRUE)))
