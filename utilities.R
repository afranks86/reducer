shrinkage_plot <- function(a, b) {
    plot(x=a, y=rep(1, length(a)), ylim=c(-0.1, 1.1))
    points(x=b, y=rep(0, length(b)))
    segments(x0=a, x1=b, y0=rep(1, length(a)), y1=rep(0, length(b)))
    abline(v=1/2, lty=2)
    abline(h=1)
    abline(h=0)
}

gen_linearY_logisticT <- function(n, p, tau, alpha, beta, mscale, escale){
    X <- matrix(rnorm(n * p), nrow=n, ncol=p)
    
    m <- mscale * X %*% alpha
    xb <- X %*% beta 
    
    e <- logistic(escale * xb)

    T <- rbinom(n, 1, e)
    Y <- m +  tau * T + rnorm(n, 0, 1)
    
    list(X=X, T=T, Y=Y, m=m, xb=xb, e=e)
}

estimate_outcome <- function(X, T, Y, cv=TRUE, Y_lambda_min=115, alpha=0){
    if(cv){
        cvglm <- cv.glmnet(cbind(T, X), Y, family="gaussian", 
                           alpha=alpha, penalty.factor = c(0, 0, rep(1, p)), intercept=FALSE)
        Y_lambda_min <- cvglm$lambda.min
    }
    
    outcome_fit <- glmnet(cbind(T, X), Y, family="gaussian", 
                       alpha=alpha, penalty.factor = c(0, rep(1, p)), intercept=FALSE, 
                       lambda=Y_lambda_min)
    alpha_hat <- coef(outcome_fit)[-c(1,2)]
    alpha_hat_normalized <- alpha_hat / sqrt(sum(alpha_hat^2))
    tau_hat <- coef(outcome_fit)[2]
    mhat0 <- predict(outcome_fit, cbind(0, X))
    mhat1 <- predict(outcome_fit, cbind(1, X))
    
    list(alpha_hat=alpha_hat,
         alpha_hat_normalized=alpha_hat_normalized,
         tau_hat=tau_hat,
         mhat0=mhat0,
         mhat1=mhat1)
}

estimate_propensity <- function(X, T, cv=TRUE, T_lambda_min=115, alpha=0){ 
    if(cv){
        cvglm <- cv.glmnet(cbind(X), T, family="binomial", 
                           alpha=alpha, penalty.factor = rep(1, p),intercept=FALSE)
        T_lambda_min <- cvglm$lambda.min
    }
    
    propensity_fit <- glmnet(cbind(X), T, family="binomial", 
                          alpha=alpha, penalty.factor = rep(1, p),intercept=FALSE, 
                          lambda=T_lambda_min)
    beta_hat <- coef(propensity_fit)[-1]
    escale_hat <- sqrt(sum(beta_hat^2))
    beta_hat_normalized <- beta_hat / escale_hat
    ehat <- predict(propensity_fit, type="response", newx=X)
    
    list(beta_hat=beta_hat,
         beta_hat_normalized=beta_hat_normalized,
         escale_hat=escale_hat,
         ehat=ehat)
}

# Grabs named attribute from list if not na; else returns default
get_attr_default <- function(thelist, attrname, default){
    if(!is.null(thelist[[attrname]])) thelist[[attrname]] else default
}

compute_gammas <- function(ab_dot_prod, mvals, mvecs, w2, times){
    p <- dim(mvecs)[1]
    # Sample random basis vectors from the null space
    null_vecs <- NullC(mvecs)[, 1:times]#sample(p, size=times, replace=FALSE)]
    # Compute coefficients for span(alpha, beta) eigenvectors
    w1 <- sqrt(as.numeric((ab_dot_prod - mvals[2] * w2^2)/mvals[1]))
    
    # Stop if magnitudes of w1 and w2 are too large, and likely not due to numerical error.
    # Inserted so that max statement in next line does not silence bugs.
    stopifnot(1 - w1^2 - w2^2 > -1e-6)
    gammas <- w1 * mvecs[, 1] + w2 * mvecs[, 2] + sqrt(max(1 - w1^2 - w2^2, 0)) * null_vecs
    gammas
}

get_bias_vec <- function(T, Y, X, xb, mvecs, mvals, ab_dot_prod, escale, w2=0, w2lim, times=10,
                      DEBUG=FALSE, alpha_hat_normalized=NA, beta_hat_normalized=NA) {
    # Compute coefficient vectors gamma corresponding to reductions
    gammas <- compute_gammas(ab_dot_prod, mvals, mvecs, w2, times)
    # Compute corresponding reductions d
    dd <- X %*% gammas
    dd_norms <- sqrt(colSums(dd^2))
    dd_normalized <- t(t(dd) / dd_norms)
    # Magnitudes of projected normalized linearized propensity onto dd
    # There is a c1 and c2 for each d
    c1s <- t(xb / sqrt(sum(xb^2))) %*% dd_normalized
    # Stop if magnitude of c1 is too large, and likely not due to numerical error.
    # Inserted so that max statement in next line does not silence bugs.
    stopifnot(all(1 - c1s^2 > -1e-6))
    c2s <- sqrt(pmax(1-c1s^2, 0))
    
    # For each d, we will be averaging over a bunch of normal samples
    normal_samples <- rnorm(500, sd=1)
    # Matrix of projected propensity scores
    e_dd <- matrix(NA, nr=dim(dd)[1], nc=dim(dd)[2])
    for(i in 1:times){
        # For each component of dd, add each value of normal_sample and evaluate logistic
        # Implement this as an outer sum between a column of dd and normal_samples
        integrand <- invlogit(escale * (outer(c1s[i] * dd[,i], c2s[i] * normal_samples, FUN="+")), a=0, b=1)
        e_dd[,i] <- rowMeans(integrand)
    }
   
    # Compute inverse weighted group means -- uses Hajek estimator
    trt_wt <- 1 / e_dd
    ctrl_wt <- 1 / (1-e_dd)
    trt_wt_means <- colSums((trt_wt * as.vector(Y))[T == 1,,drop=FALSE]) /
        colSums(trt_wt[T == 1,,drop=FALSE])
    ctrl_wt_means <- colSums((ctrl_wt * as.vector(Y))[T == 0,,drop=FALSE]) /
        colSums(ctrl_wt[T == 0,,drop=FALSE])
    list(bias0=mean(ctrl_wt_means), bias1=mean(trt_wt_means))
    
    # Debugging leftovers
    #list(gammas=gammas, e_dd=e_dd, ests=ests, integrand=save_integrand,
    #     c1s=c1s, dd=dd, dd_normalized=dd_normalized)
}

get_bias_old <- function(T, Y, X, xb, mvecs, mvals, ab_dot_prod, escale, w2=0, w2lim, times=10,
                     DEBUG=FALSE, alpha_hat_normalized=NA, beta_hat_normalized=NA) {
    
    N <- NullC(mvecs)
    
    # Debugging leftovers 
    #reductions <- get_bias_vec(T, Y, X, xb, mvecs, mvals, ab_dot_prod, escale, w2, w2lim, times)
    #gg <- reductions$gammas
    #e_dd <- reductions$e_dd
    #ests <- reductions$ests
    #debug_integrand <- reductions$integrand

    w1 <- sqrt(as.numeric((ab_dot_prod - mvals[2] * w2^2)/mvals[1]))
    bias0 <- bias1 <- 0
    normal_samples <- rnorm(500, sd=1)

    for(i in 1:times) {
        # Stop if magnitudes of w1 and w2 are too large, and likely not due to numerical error.
        # Inserted so that max statement in next line does not silence bugs.
        stopifnot(1 - w1^2 - w2^2 > -1e-6)
        g <- w1 * mvecs[, 1] + w2 * mvecs[, 2] + sqrt(max(1 - w1^2 - w2^2, 0)) * N[, i]# %*% rustiefel(p-2, 1)
        d <- Re(X %*% g)
        
        if(DEBUG){
            if(is.na(alpha_hat_normalized) || is.na(beta_hat_normalized))
                stop("Debug requires alpha_hat_normalized and beta_hat_normalized vectors.")
            mhatX <- X %*% alpha_hat_normalized
            eX <- X %*% beta_hat_normalized
            
            print(sprintf("Partial correlation: %s",
                          round(cor(lm(mhatX ~ d)$resid, lm(eX ~ d)$resid), 4)))
            
            print(sprintf("Hyperbola condition: %s",
                          round((alpha_hat_normalized %*% as.vector(g)) *
                                  (as.vector(g) %*% beta_hat_normalized) -
                                    (alpha_hat_normalized %*% beta_hat_normalized), 4)))
           
            if(w2 == -sqrt(abs(mvals[2])))
              print(sprintf("Checking that w2 = min corresponds to e(X): %s",
                            round(cor(g, beta_hat_normalized), 4)))
            if(w2 == sqrt(abs(mvals[2])))
              print(sprintf("Checking that w2 = max corresponds to mhat(X): %s",
                            round(cor(g, alpha_hat_normalized), 4)))
            if(w2 == 0)
              print(sprintf("Checking that w2 = 0 correspond to bisector: %s",
                            round(sum(alpha_hat_normalized %*% as.vector(g) -
                                          beta_hat_normalized %*% as.vector(g)), 4)))
        }

        c1 <- as.numeric(t(xb/sqrt(sum(xb^2))) %*% d/sqrt(sum(d^2)))
        # Stop if magnitude of c1 is too large, and likely not due to numerical error.
        # Inserted so that max statement in next line does not silence bugs.
        stopifnot(1 - c1^2 > -1e-6)
        c2 <- sqrt(max(1-c1^2, 0))

        e_d <- sapply(d, function(di) {
            mean(invlogit(escale*(c1*di + c2*normal_samples), a=0, b=1))
        })
        
        ## res <- glm(T ~ d - 1, family="binomial")
        ## e_d <- predict(res, type="response")
        
        bias0 <- bias0 + sum(1 / (1-e_d[T==0]) * Y[T==0]) / sum(1 / (1-e_d[T==0]))
        bias1 <- bias1 + sum(1 / e_d[T==1] * Y[T==1]) / sum(1 / e_d[T==1])

        if(is.na(bias0) | is.na(bias1))
            browser()
        #biases[i] <- mybias1 - mybias0
        
    }

    #print(summary(biases))
    list(bias0=bias0/times, bias1=bias1/times)
}

invlogit <- function(x, a, b) {
    exp(a+x*b) / (1+exp(a+x*b))
}

logistic <- function(x){
    exp(x) / (1 + exp(x))
    #1 / (1 + exp(-x))
}

ipw_est <- function(e, T, Y, hajek=FALSE){
  wts <- (T / e) - ((1 - T) / (1 - e))
  if(hajek) {
    # Note that negative weights for control assignments are flipped to positive by normalization
    sum(wts[T == 1] * Y[T == 1]) / sum(wts[T == 1]) - sum(wts[T == 0] * Y[T == 0]) / sum(wts[T == 0])
  } else {
    mean(wts * Y)
  }
}