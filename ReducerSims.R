library(mvtnorm)
library(rstiefel)
library(glmnet)
library(lubridate)
library(R.utils)
source("utilities.R")

argv <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

EST_OUTCOME <- as.logical(get_attr_default(argv, "est_outcome", TRUE))
OUTCOME_CV <- as.logical(get_attr_default(argv, "outcome_cv", TRUE))
Y_LAMBDA <- as.numeric(get_attr_default(argv, "y_lambda", 115))
Y_ALPHA <- as.numeric(get_attr_default(argv, "y_alpha", 1))

EST_PROPENSITY <- as.logical(get_attr_default(argv, "est_propensity", FALSE))
PROP_CV <- as.logical(get_attr_default(argv, "prop_cv", TRUE))
T_LAMBDA <- as.numeric(get_attr_default(argv, "t_lambda", 115))
T_ALPHA <- as.numeric(get_attr_default(argv, "t_alpha", 1))

coef_setting <- as.numeric(get_attr_default(argv, "coef", 1))
mscale <- as.numeric(get_attr_default(argv, "mscale", 15))
escale <- as.numeric(get_attr_default(argv, "escale", 3))

times <- as.numeric(get_attr_default(argv, "times", 50))
iters <- as.numeric(get_attr_default(argv, "iters", 50))

n <- as.numeric(get_attr_default(argv, "n", 1000))
p <- as.numeric(get_attr_default(argv, "p", 1000))

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

w2_scale_vec <- c(seq(1, -1, by=-.1))

true_ate_vec <- tau_hat_vec <- ipw_vec <- ipw_d_vec <- aipw_vec <- aipw_d_vec <- naive_vec <- numeric(iters)
results_array <- array(dim=c(iters, 6, length(w2_scale_vec)))


for(iter  in 1:iters) {

    ## ################
    ## Generate dataset
    ## #################
    
    if(coef_setting == 1){
        alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        beta <- c(1, 0, 1, rep(0, p-3))/sqrt(2)
    } else {
        alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        beta <- c(-1, 0.75, 1, rep(0, p-3))/sqrt(3)
    }
    cor(alpha, beta)
    
    true_ate <- tau <- 5
    
    simdat <- gen_linearY_logisticT(n, p, tau, alpha, beta, mscale, escale)
    
    X <- simdat$X
    T <- simdat$T
    Y <- simdat$Y
    m <- simdat$m
    e <- simdat$e
    
    ## ################
    ## Set nuisance parameters
    ## #################
    
    out_ests <- if(EST_OUTCOME) estimate_outcome(X, T, Y, cv=OUTCOME_CV, Y_lambda_min=Y_LAMBDA, alpha=Y_ALPHA)
                else list()
   
    alpha_hat <- get_attr_default(out_ests, "alpha_hat", alpha)
    alpha_hat_normalized <- get_attr_default(out_ests, "alpha_hat_normalized",
                                             alpha / sqrt(sum(alpha^2)))
    mhat0 <- get_attr_default(out_ests, "mhat0", X %*% alpha)
    mhat1 <- get_attr_default(out_ests, "mhat1", X %*% alpha + tau)
    tau_hat <- get_attr_default(out_ests, "tau_hat", tau)
    
    prop_ests <- if(EST_PROPENSITY) estimate_propensity(X, T, cv=PROP_CV, T_lambda_min=T_LAMBDA, alpha=T_ALPHA)
                 else list()
    
    beta_hat <- get_attr_default(prop_ests, "beta_hat", beta * escale)
    beta_hat_normalized <- get_attr_default(prop_ests, "beta_hat_normalized", beta)
    escale_hat <- get_attr_default(prop_ests, "escale_hat", escale)
    ehat <- get_attr_default(prop_ests, "ehat", e)
    xb <- X %*% beta_hat_normalized
   
    ## ################
    ## Compute hyperbola
    ## #################
    
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
   
    # Depending on correlation between alpha and beta, switch eigenspace labels to
    # traverse the "closed" side of the hyperbola as we vary the w2 parameter
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

        residual <- Y - X %*% alpha_hat - T * tau_hat
        
        ## Naive difference in means
        ate_naive <- mean(Y[T == 1]) - mean(Y[T == 0])
        
        ## Standard IPW
        ate_ipw <- ipw_est(ehat, T, Y, hajek=TRUE)
        
        ## Compute IPW_d, using reduced d 
        ## w2 is the tuning parameter for how similar d is to to ehat vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger should reduce variance)
        bias <- get_bias(T=T, Y=Y, X=X, xb=xb, mvecs=mvecs, mvals=mvals,
                         ab_dot_prod=ab_dot_prod, escale=escale_hat,
                         w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug,
                         alpha_hat_normalized=alpha_hat_normalized, beta_hat_normalized=beta_hat_normalized)

        ate_ipw_d <- bias$bias1 - bias$bias0

        ## Compute standard AIPW using known true e  
        #mu_treat <- mean(X %*% alpha_hat + tau_hat) + sum(1 / e[T==1] * residual[T==1]) / sum(1 / e[T==1])
        #mu_ctrl <- mean(X %*% alpha_hat) + sum(1/ (1-e[T==0]) * residual[T==0]) / sum(1 / (1-e[T==0]))
        #ate_aipw <- mu_treat - mu_ctrl
        ate_aipw <- mean(X %*% alpha_hat + tau_hat) - mean(X %*% alpha_hat) + ipw_est(ehat, T, residual, hajek=TRUE)
        
        ## Compute negative regression bias by balancing on residuals
        ## w2 is the tuning parameter for how similar to e vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger shoudl reduce variance)
        bias <- get_bias(T=T, Y=residual, X=X, xb=xb, mvecs=mvecs, mvals=mvals,
                         ab_dot_prod=ab_dot_prod, escale=escale_hat,
                         w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug,
                         alpha_hat_normalized=alpha_hat_normalized, beta_hat_normalized=beta_hat_normalized)
        ## Correct bias and compute AIPW_d
        mu_treat <- mean(X %*% alpha_hat + tau_hat) + bias$bias1
        mu_ctrl <- mean(X %*% alpha_hat) + bias$bias0
        ate_aipw_d <- mu_treat - mu_ctrl

        print(sprintf("Naive: %.3f, Reg: %.3f, IPW: %.3f, IPW-d: %.3f, AIPW: %.3f, AIPW-d: %.3f",
                      ate_naive, tau_hat, ate_ipw, ate_ipw_d, ate_aipw, ate_aipw_d))
        
        true_ate_vec[iter] <- true_ate
        tau_hat_vec[iter] <- tau_hat
        ipw_vec[iter] <- ate_ipw
        ipw_d_vec[iter] <- ate_ipw_d
        aipw_vec[iter] <- ate_aipw
        aipw_d_vec[iter] <- ate_aipw_d
        naive_vec[iter] <- ate_naive

        results_array[iter, 1, j] <- tau_hat
        results_array[iter, 2, j] <- ate_ipw
        results_array[iter, 3, j] <- ate_ipw_d
        results_array[iter, 4, j] <- ate_aipw
        results_array[iter, 5, j] <- ate_aipw_d
        results_array[iter, 6, j] <- ate_naive
        
    }

                                        # rmse
    sqrt(apply((rbind(tau_hat_vec, ipw_vec, ipw_d_vec, aipw_vec, aipw_d_vec) - true_ate_vec)^2, 1, 
               function(x) mean(x, na.rm=TRUE)))

                                        # mean absolute error
    apply(abs(rbind(tau_hat_vec, ipw_vec, ipw_d_vec, aipw_vec, aipw_d_vec) - true_ate_vec), 1, 
          function(x) median(x, na.rm=TRUE))

                                        # bias
    apply(rbind(tau_hat_vec, ipw_vec, ipw_d_vec, aipw_vec, aipw_d_vec) - true_ate_vec, 1, 
          function(x) mean(x, na.rm=TRUE))

                                        # neg/pos bias fraction
    apply((sign(rbind(tau_hat_vec, ipw_vec, ipw_d_vec, aipw_vec, aipw_d_vec) - true_ate_vec)+1)/2, 1, 
          function(x) mean(x, na.rm=TRUE))

    #save(results_array,
    #     file=sprintf("results_negn%i_p%i_escale%.1f_wmax%.2f_wmin%.2f_%s.RData",
    #                  n, p ,escale,
    #                  max(w2_scale_vec), min(w2_scale_vec),
    #                  today()))
}

dimnames(results_array) <- list(1:iters, c("Regression", "IPW", "IPW_d", "AIPW", "AIPW_d", "Naive"), w2_scale_vec)
sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

apply(abs(results_array - true_ate), c(2, 3), function(x) median(x, na.rm=TRUE))

save(results_array, true_ate,
     file=sprintf("results/results_n%i_p%i_coef%i_escale%.1f_mscale%.1f_yalpha%i_estpropensity%s_%s.RData",
                  n, p, coef_setting, escale, mscale, Y_ALPHA, EST_PROPENSITY,
                  gsub(" ", "", now(), fixed=TRUE)))