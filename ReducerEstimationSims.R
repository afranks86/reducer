library(balanceHD)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(lubridate)
source("utilities.R")

## w2_scale_vec <- c(seq(0.95, 0.5, by=-0.05))
## w2_scale_vec <- c(w2_scale_vec, 0, -rev(w2_scale_vec))
eigen_debug <- FALSE
bias_debug <- FALSE
bias_times <- if(bias_debug) 1 else 50

w2_scale_vec <- c(seq(1, -1, by=-.1))
iters <- 20

true_ate_vec <- tau_hat_vec <- ipw_vec <- ipw_d_vec <- aipw_vec <- aipw_d_vec <- numeric(iters)
results_array <- array(dim=c(iters, 5, length(w2_scale_vec)))

for(iter  in 1:iters) {

    ## ################
    ## Generate dataset
    ## #################
    
    n = 1000
    p = 1000
    
    alpha <- rustiefel(p, 1) 
    beta <- rustiefel(p, 1) 
    
    alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
    beta <- c(1, 0, 1, rep(0, p-3))/sqrt(2)
    cor(alpha, beta)
    
    ## X <- rmvnorm(n, rep(0, p), diag(1, p))
    X <- matrix(rnorm(n * p), nrow=n, ncol=p)
    mscale <- -15
    
    m <- mscale * X %*% alpha
    xb <- X %*% beta 
    escale <- 3
    
    e <- exp(escale*xb)/(1 + exp(escale*xb))

    tau = 5
    T <- rbinom(n, 1, e)
    Y <- m +  tau * T + rnorm(n, 0, 1)
    
    true_ate <- tau


    cvglm <- cv.glmnet(cbind(T, X), Y, family="gaussian", 
                       alpha=0, penalty.factor = c(0, 0, rep(1, p)), intercept=FALSE)
    Y_lambda_min <- cvglm$lambda.min

    
    glm_coefs <- coef(glmnet(cbind(T, X), Y, family="gaussian", 
                             alpha=0, penalty.factor = c(0, rep(1, p)),intercept=FALSE, 
                             lambda=Y_lambda_min))

    cvglm <- cv.glmnet(cbind(X), T, family="binomial", 
                    alpha=0, penalty.factor = rep(1, p),intercept=FALSE)
    T_lambda_min <- cvglm$lambda.min

    T_glm_coefs <- coef(glmnet(cbind(X), T, family="binomial", 
                               alpha=0, penalty.factor = rep(1, p),intercept=FALSE, 
                               lambda=T_lambda_min))[-1]

    ehat = predict(glmnet(cbind(X), T, family="binomial", 
                          alpha=0, penalty.factor = rep(1, p),intercept=FALSE, 
                          lambda=T_lambda_min), type="response", newx=X)

    alpha_hat <- glm_coefs[-c(1:2)]
    beta_hat <- T_glm_coefs
    tau_hat <- glm_coefs[2]    
    
    alpha_hat_norm <- sqrt(sum(alpha_hat^2))
    beta_hat_norm <- sqrt(sum(beta_hat^2))
    
    alpha_hat_normalized <- alpha_hat / alpha_hat_norm
    beta_hat_normalized <- beta_hat / beta_hat_norm

    ab_dot_prod <- as.numeric(t(alpha_hat_normalized) %*% beta_hat_normalized) 
    if(eigen_debug){
      ab_outer <- alpha_hat_normalized %*% t(beta_hat_normalized)
      spec <- eigen(0.5 * (ab_outer + t(ab_outer)))
      debug_mvals <- spec$values[c(1,p)]
      debug_mvecs <- spec$vectors[,c(1,p)]
    }

    mvecs <- cbind((alpha_hat_normalized + beta_hat_normalized)/sqrt(2 + 2 * ab_dot_prod),
    (alpha_hat_normalized - beta_hat_normalized)/sqrt(2 - 2 * ab_dot_prod))
    mvals <- c((ab_dot_prod + 1)/2, (ab_dot_prod - 1)/2)
   
    # Depending on correlation between alpha and beta, switch eigenspace labels to
    # traverse the "closed" side of the hyperbola as we vary the w2 parameter
    if(ab_dot_prod < 0){
        mvals <- mvals[2:1]
        mvecs <- mvecs[, 2:1]
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
        
        ## Standard IPW with known true e
        ate_ipw <- ipw_est(ehat, T, Y, hajek=TRUE)
        
        ## Compute IPW_d, using reduced d 
        ## w2 is the tuning parameter for how similar to e vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger should reduce variance)
        bias <- get_bias(T=T, Y=Y, X=X, xb=X %*% beta_hat, mvecs=mvecs, mvals=mvals,
                         ab_dot_prod=ab_dot_prod, escale=beta_hat_norm,
                         w2=w2, w2lim=w2_lim, times=bias_times, DEBUG=bias_debug)

        ate_ipw_d <- bias$bias1 - bias$bias0

        ## Compute standard AIPW using known true e  
        mu_treat <- mean(X %*% alpha_hat + tau_hat) + sum(1 / ehat[T==1] * residual[T==1]) / sum(1 / ehat[T==1])
        mu_ctrl <- mean(X %*% alpha_hat) + sum(1/ (1-ehat[T==0]) * residual[T==0]) / sum(1 / (1-ehat[T==0]))
        ate_aipw <- mu_treat - mu_ctrl
        
        ## Compute negative regression bias by balancing on residuals
        ## w2 is the tuning parameter for how similar to e vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger shoudl reduce variance)
        bias <- get_bias(T=T, Y=residual, X=X, xb=X %*% beta_hat,
                         mvecs=mvecs, mvals=mvals,
                         ab_dot_prod=ab_dot_prod, w2=w2, escale=beta_hat_norm,
                         w2lim=w2_lim, times=bias_times, DEBUG=bias_debug)

        ## Correct bias and compute AIPW_d
        mu_treat <- mean(X %*% alpha_hat + tau_hat) + bias$bias1
        mu_ctrl <- mean(X %*% alpha_hat) + bias$bias0
        ate_aipw_d <- mu_treat - mu_ctrl

        print(sprintf("Reg: %.3f, IPW: %.3f, IPW-d: %.3f, AIPW: %.3f, AIPW-d: %.3f",
                      tau_hat, ate_ipw, ate_ipw_d, ate_aipw, ate_aipw_d))
        
        true_ate_vec[iter] <- true_ate
        tau_hat_vec[iter] <- tau_hat
        ipw_vec[iter] <- ate_ipw
        ipw_d_vec[iter] <- ate_ipw_d
        aipw_vec[iter] <- ate_aipw
        aipw_d_vec[iter] <- ate_aipw_d

        results_array[iter, 1, j] <- tau_hat
        results_array[iter, 2, j] <- ate_ipw
        results_array[iter, 3, j] <- ate_ipw_d
        results_array[iter, 4, j] <- ate_aipw
        results_array[iter, 5, j] <- ate_aipw_d
        
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

    save(results_array,
         file=sprintf("beta_est_n%i_p%i_escale%.1f_wmax%.2f_wmin%.2f_%s.RData",
                      n, p ,escale,
                      max(w2_scale_vec), min(w2_scale_vec),
                      today()))
}

dimnames(results_array) <- list(1:iters, c("Regression", "IPW", "IPW_d", "AIPW", "AIPW_d"), w2_scale_vec)
sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

apply(abs(results_array - true_ate), c(2, 3), function(x) median(x, na.rm=TRUE))

save(results_array,
     file=sprintf("beta_est_n_n%i_p%i_escale%.1f_wmax%.2f_wmin%.2f_%s.RData",
                  n, p , escale,
                  max(w2_scale_vec), min(w2_scale_vec),
                  today()))

library(tidyverse)
library(ggridges)
library(superheat)
library(patchwork)
as_tibble(results_array[, , "0.6"]) %>%
    gather(key=Type, value=Estimate) %>%
    ggplot() +
    geom_density_ridges(aes(x=Estimate, y=Type, fill=Type), stat="binline") +
    geom_vline(xintercept=5) + theme_bw()


rmse_mat <- sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))
bias_mat <- apply(results_array - true_ate, c(2, 3), function(x) mean(x, na.rm=TRUE))
var_mat <- rmse_mat^2 - bias_mat^2


tib <- as_tibble(t(rmse_mat))
tib$value <- as.numeric(colnames(rmse_mat))
rmse_plot <- tib %>% gather(key=Type, value=RMSE, -value) %>%
    ggplot() + geom_line(aes(x=value, y=RMSE, col=Type)) +
    theme_bw() + theme(legend.position="none")

tib <- as_tibble(t(bias_mat))
tib$value <- as.numeric(colnames(bias_mat))
bias_plot <- tib %>% gather(key=Type, value=Bias, -value) %>%
    ggplot() + geom_line(aes(x=value, y=Bias, col=Type)) +
    theme_bw() + theme(legend.position="none")

tib <- as_tibble(t(var_mat))
tib$value <- as.numeric(colnames(bias_mat))
var_plot <- tib %>% gather(key=Type, value=Var, -value) %>%
    ggplot() + geom_line(aes(x=value, y=Var, col=Type)) + theme_bw() 


rmse_plot + bias_plot + var_plot

cor(results_array[, "IPW_d",  ])
superheat::superheat(cor(results_array[, "AIPW_d",  ]))

cor(results_array[, ,  "0.9"])
cor(results_array[, ,  "0.7"])



