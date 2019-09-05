library(balanceHD)
library(mvtnorm)
library(rstiefel)
library(glmnet)
source("utilities.R")

w2_scale_vec <- c(1, seq(0.95, 0.5, by=-0.05))
w2_scale_vec <- c(w2_scale_vec, 0, -rev(w2_scale_vec))
iters <- 100

true_ate_vec <- tau_hat_vec <- ipw_vec <- ipw_d_vec <- aipw_vec <- aipw_d_vec <- numeric(iters)
results_array <- array(dim=c(iters, 5, length(w2_scale_vec)))

for(j in 1:length(w2_scale_vec)) {
    w2scale <- w2_scale_vec[j]

    for(iter  in 1:iters) {
        print(paste(w2scale, iter, sep=", "))
        n = 1000
        p = 1000
        
        alpha <- rustiefel(p, 1) 
        beta <- rustiefel(p, 1) 
        
        alpha <- c(1, -1, 0, rep(0, p-3))/sqrt(2)
        beta <- c(1, 0, 1, rep(0, p-3))/sqrt(2)
        cor(alpha, beta)
        
        X <- rmvnorm(n, rep(0, p), diag(1, p))
        mscale <- 15
        
        m <- mscale * X %*% alpha
        xb <- X %*% beta 
        escale <- 1/2
        
        e <- exp(escale*xb)/(1 + exp(escale*xb))

        tau = 5
        T <- rbinom(n, 1, e)
        Y <- m +  tau * T + rnorm(n, 0, 1)
        
        true_ate <- tau
        
        ## Y0 <- Y[T==0]
        ## Y1 <- Y[T==1]
        
        ## mean(Y0)
        ## mean(Y0 / (1 - e[T==0]))
        
        ## e_m <- predict(glm(T ~ m))
        ## mean(T * Y / e_m) - mean((1-T) * Y / (1-e_m))
        ## mean(T * Y / e) - mean((1-T) * Y / (1-e))
        
        
        Y_lambda_min <- cv.glmnet(cbind(T, X), Y, family="gaussian", 
                                  alpha=0, penalty.factor = c(0, rep(1, p)), intercept=FALSE)$lambda.min
        glm_coefs <- coef(glmnet(cbind(T, X), Y, family="gaussian", 
                                 alpha=0, penalty.factor = c(0, rep(1, p)),intercept=FALSE, 
                                 lambda=Y_lambda_min))
        
        ehat <- predict(glmnet(X, T, family="binomial", alpha=0), newx=X, type="response")
        
        alpha_hat <- glm_coefs[-c(1:2)]
        tau_hat <- glm_coefs[2]    
        
        alpha_hat_norm <- sqrt(sum(alpha_hat^2))
        alpha_hat_normalized <- alpha_hat / alpha_hat_norm
        
        M <- (alpha_hat_normalized %*% t(beta) + beta %*% t(alpha_hat_normalized))/2
        eig <- eigen(M)
        mvecs <- eig$vectors[, c(1,p)]
        mvals <- eig$values[c(1,p)]
        if(abs(mvals[2]) > mvals[1]) {
            mvals <- -mvals[2:1]
            mvecs <- mvecs[, 2:1]
        }
        
        ab_dot_prod <- abs(t(alpha_hat_normalized) %*% beta)
        
        w2_beta_lim <- (t(beta) %*% mvecs)[2]
        w2_alpha_lim <- (t(alpha_hat_normalized) %*% mvecs)[2]
        
        w2 <- w2scale*w2_lim

        ## udv <- svd(X %*% N %*% t(N))
        ## v <- udv$v[, 1]
        ## var(X %*% N %*% rustiefel(p-2, 1))
        ## g <- w1*mvecs[, 1] + lw2 * mvecs[, 2] + sqrt(1-w1^2 - w2^2) * v
        ## g <- w1*mvecs[, 1] + w2 * mvecs[, 2] + sqrt(1-w1^2 - w2^2) * N %*% rustiefel(p-2, 1)
        
        ## d <- Re(X %*% g)
        ## res <- glm(T ~ d, family="binomial")
        ## e_d <- predict(res, type="response")
        ## sum(1 / e_d[T==1] * Y[T==1]) / sum(1 / e_d[T==1]) - sum(1 / (1-e_d[T==0]) * Y[T==0]) / sum(1 / (1-e_d[T==0])) 

        N <- NullC(mvecs)
        if(abs(w2) == abs(w2_lim))
            g <- w1*mvecs[, 1] + w2 * mvecs[, 2]
        else 
            g <- w1*mvecs[, 1] + w2 * mvecs[, 2] + sqrt(1-w1^2 - w2^2) * N %*% rustiefel(p-2, 1)
        
        if(w2==0 & (abs(t(g) %*% beta) - abs(t(g) %*% alpha_hat_normalized) > 1e-10))
            browser()
        
        d <- Re(X %*% g)
        d_a <- X %*% alpha_hat_normalized
        
        res <- glm(T ~ d, family="binomial")
        e_d <- predict(res, type="response")
        
        residual <- Y - X %*% alpha_hat - T * tau_hat
        
        ## Standard IPW with known true e
        ate_ipw <- sum(1 / e[T==1] * Y[T==1]) / sum(1 / e[T==1]) -  sum(1 / (1-e[T==0]) * Y[T==0]) / sum(1 / (1-e[T==0]))
        
        ## Compute IPW_d, using reduced d 
        ## w2 is the tuning parameter for how similar to e vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger should reduce variance)
        bias <- get_bias(T=T, residual=Y, alpha_hat_normalized=alpha_hat_normalized, 
                         beta=beta, tau_hat=tau_hat, w2=w2, times=1000)
        ate_ipw_d <- bias$bias1 - bias$bias0

        ## Compute standard AIPW using known true e  
        mu_treat <- mean(X %*% alpha_hat + tau_hat) + sum(1 / e[T==1] * residual[T==1]) / sum(1 / e[T==1])
        mu_ctrl <- mean(X %*% alpha_hat) + sum(1/ (1-e[T==0]) * residual[T==0]) / sum(1 / (1-e[T==0]))
        ate_aipw <- mu_treat - mu_ctrl
        
        ## Compute negative regression bias by balancing on residuals
        ## w2 is the tuning parameter for how similar to e vs mhat,
        ## times is the number of reductions to use to compute weighted estimates (larger shoudl reduce variance)
        bias <- get_bias(T=T, residual=residual, alpha_hat_normalized=alpha_hat_normalized, 
                         beta=beta, tau_hat=tau_hat, w2=w2, times=1000)
        
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

}

dimnames(results_array) <- list(1:iters, c("Regression", "IPW", "IPW_d", "AIPW", "AIPW_d"), w2_scale_vec)
sqrt(apply((results_array - true_ate)^2, c(2, 3), function(x) mean(x, na.rm=TRUE)))

apply(abs(results_array - true_ate), c(2, 3), function(x) median(x, na.rm=TRUE))

save(results_array, file=sprintf("results_n%i_p%i_escale%.1f.RData", n, p ,escale))

library(tidyverse)
library(ggridges)
as_tibble(results_array[, , "0.5"]) %>%
    gather(key=Type, value=Estimate) %>%
    ggplot() +
    geom_density_ridges(aes(x=Estimate, y=Type, fill=Type), stat="binline") +
    geom_vline(xintercept=5) + theme_bw()

rmse_mat <- apply(abs(results_array - true_ate), c(2, 3), function(x) median(x, na.rm=TRUE))


tibble(rmse_mat) 





