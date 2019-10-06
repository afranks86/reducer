shrinkage_plot <- function(a, b) {
    plot(x=a, y=rep(1, length(a)), ylim=c(-0.1, 1.1))
    points(x=b, y=rep(0, length(b)))
    segments(x0=a, x1=b, y0=rep(1, length(a)), y1=rep(0, length(b)))
    abline(v=1/2, lty=2)
    abline(h=1)
    abline(h=0)
}

gen_linearY_logisticT <- function(n, p, tau, alpha, beta, mscale, escale, sigma2_y){
    X <- matrix(rnorm(n * p), nrow=n, ncol=p)
    
    m <- mscale * X %*% alpha
    xb <- X %*% beta 
    
    e <- logistic(escale * xb)

    T <- rbinom(n, 1, e)
    Y <- m +  tau * T + rnorm(n, 0, sigma2_y)
    
    list(X=X, T=T, Y=Y, m=m, xb=xb, e=e)
}

estimate_outcome <- function(X, T, Y, estimand, cv=TRUE, Y_lambda_min=115, alpha=0){
 
    if(estimand == "ATC") {
        Xfit <- X[T==0, ]
        Yfit <- Y[T==0]
        Xpred <- X[T==1, ]
    } else if (estimand == "ATT"){
        Xfit <- X[T==1, ]
        Yfit <- Y[T==1]
        Xpred <- X[T==0, ]
    } else {
        Xfit <- cbind(T, X)
        Yfit <- Y
    }
    
    if(cv){
        if(estimand == "ATC" | estimand == "ATT") {
            penalty.factor <- rep(1, p)
            cvglm <- cv.glmnet(Xfit, Yfit, family="gaussian", 
                               alpha=alpha,
                               penalty.factor=penalty.factor,
                               intercept=FALSE)
            Y_lambda_min <- cvglm$lambda.1se

        } else {

            ## For ATE regression doesn't penalize T coefficient
            penalty.factor <- c(0, rep(1, p))
            cvglm <- cv.glmnet(Xfit, Yfit, family="gaussian", 
                               alpha=alpha, penalty.factor = penalty.factor,
                               intercept=FALSE)
            Y_lambda_min <- cvglm$lambda.1se
        }


    }
    
    outcome_fit <- glmnet(Xfit, Yfit, family="gaussian", 
                          alpha=alpha, penalty.factor = penalty.factor,
                          intercept=FALSE, 
                          lambda=Y_lambda_min)

    if(estimand == "ATC" | estimand == "ATT") {
        alpha_hat <- coef(outcome_fit)[-1]
        alpha_hat_normalized <- alpha_hat / sqrt(sum(alpha_hat^2))

        mhat0 <- predict(outcome_fit, cbind(X[T==0, ]))
        mhat1 <- predict(outcome_fit, cbind(X[T==1, ]))
        
    }  else {
        alpha_hat <- coef(outcome_fit)[-c(1,2)]
        alpha_hat_normalized <- alpha_hat / sqrt(sum(alpha_hat^2))
        mhat0 <- predict(outcome_fit, cbind(0, X))
        mhat1 <- predict(outcome_fit, cbind(1, X))        
    }
    tau_hat <- mean(mhat1) - mean(mhat0)
    
    list(alpha_hat=alpha_hat,
         alpha_hat_normalized=alpha_hat_normalized,
         tau_hat=tau_hat,
         mhat0=mhat0,
         mhat1=mhat1,
         lambda=Y_lambda_min)
}

estimate_propensity <- function(X, T, cv=TRUE, T_lambda_min=115,
                                eta=0.1, alpha=0){ 
    if(cv){
        cvglm <- cv.glmnet(cbind(X), T, family="binomial", 
                           alpha=alpha, penalty.factor = rep(1, p), intercept=FALSE)
        T_lambda_min <- cvglm$lambda.1se
    }

    if(eta > 1 | eta < 0) {
        eta <- 0
        warning(sprintf("eta  = %f", eta))
    }

    propensity_fit <- glmnet(cbind(X), T, family="binomial", 
                          alpha=alpha, penalty.factor = rep(1, p),intercept=FALSE, 
                          lambda=T_lambda_min)
    beta_hat <- coef(propensity_fit)[-1]
    escale_hat <- sqrt(sum(beta_hat^2))

    beta_hat_normalized <- if(escale_hat > 0) {
                               beta_hat / escale_hat
                           } else {
                               beta_hat 
                           }
    ehat <- predict(propensity_fit, type="response", newx=X)

    ## Clipping
    ehat_clip <- ehat
    ehat_clip[ehat < eta] <- eta
    ehat_clip[ehat > 1 - eta] <- 1 - eta
    
    list(beta_hat=beta_hat,
         beta_hat_normalized=beta_hat_normalized,
         escale_hat=escale_hat,
         ehat=ehat,
         ehat_clip=ehat_clip,
         lambda=T_lambda_min)
}

## Grabs named attribute from list if not na; else returns default
get_attr_default <- function(thelist, attrname, default){
    if(!is.null(thelist[[attrname]])) thelist[[attrname]] else default
}

compute_gammas <- function(ab_dot_prod, mvals, mvecs, w2, times){
    p <- dim(mvecs)[1]
    ## Sample random basis vectors from the null space
    null_vecs <- NullC(mvecs)[, 1:times]#sample(p, size=times, replace=FALSE)]
    ## Compute coefficients for span(alpha, beta) eigenvectors
    w1 <- sqrt(as.numeric((ab_dot_prod - mvals[2] * w2^2)/mvals[1]))
    
    ## Stop if magnitudes of w1 and w2 are too large, and likely not due to numerical error.
    ## Inserted so that max statement in next line does not silence bugs.
    stopifnot(1 - w1^2 - w2^2 > -1e-6)
        
    gammas <- w1 * mvecs[, 1] + w2 * mvecs[, 2] + sqrt(max(1 - w1^2 - w2^2, 0)) * null_vecs
    gammas
}

get_bias_vec <- function(T, Y, X, xb, estimand,
                         mvecs, mvals, ab_dot_prod, escale, w2=0, w2lim, times=10,
                         DEBUG=FALSE, alpha_hat_normalized=NA, beta_hat_normalized=NA) {
    
    ### Compute coefficient vectors gamma corresponding to reductions

    gammas <- compute_gammas(ab_dot_prod, mvals, mvecs, w2, times)

    ## Compute corresponding reductions d
    dd <- X %*% gammas
    dd_norms <- sqrt(colSums(dd^2))
    dd_normalized <- t(t(dd) / dd_norms)

    ## Magnitudes of projected normalized linearized propensity onto dd
    ## There is a c1 and c2 for each d
    c1s <- t(xb / sqrt(sum(xb^2))) %*% dd_normalized

    ## Stop if magnitude of c1 is too large, and likely not due to numerical error.
    ## Inserted so that max statement in next line does not silence bugs.
    stopifnot(all(1 - c1s^2 > -1e-6))
    c2s <- sqrt(pmax(1-c1s^2, 0))
    
    ## For each d, we will be averaging over a bunch of normal samples
    normal_samples <- rnorm(500, sd=1)
    ## Matrix of projected propensity scores
    e_dd <- matrix(NA, nr=dim(dd)[1], nc=dim(dd)[2])
    for(i in 1:times){
        ## For each component of dd, add each value of normal_sample and evaluate logistic
        ## Implement this as an outer sum between a column of dd and normal_samples
        integrand <- invlogit(escale * (outer(c1s[i] * dd[,i], c2s[i] * normal_samples, FUN="+")), a=0, b=1)
        e_dd[,i] <- rowMeans(integrand)
    }
   
    ## Compute inverse weighted group means -- uses Hajek estimator

    if(estimand == "ATT") {
        trt_wt <- matrix(1, nrow=nrow(e_dd), ncol=ncol(e_dd))
        ctrl_wt <- e_dd / (1-e_dd)
    }  else if (estimand == "ATC") {
        trt_wt <- (1 - e_dd) / e_dd
        ctrl_wt <- matrix(1, nrow=nrow(e_dd), ncol=ncol(e_dd))
    } else {
        trt_wt <- 1 / e_dd
        ctrl_wt <- 1 / (1 - e_dd) 
    }
    
    trt_wt_means <- colSums((trt_wt * as.vector(Y))[T == 1,,drop=FALSE]) /
        colSums(trt_wt[T == 1,,drop=FALSE])
    ctrl_wt_means <- colSums((ctrl_wt * as.vector(Y))[T == 0,,drop=FALSE]) /
        colSums(ctrl_wt[T == 0,,drop=FALSE])
    list(bias0=mean(ctrl_wt_means), bias1=mean(trt_wt_means))
    
}

invlogit <- function(x, a, b) {
    exp(a+x*b) / (1+exp(a+x*b))
}

logistic <- function(x){
    exp(x) / (1 + exp(x))
    #1 / (1 + exp(-x))
}

ipw_est <- function(e, T, Y, estimand, hajek=FALSE){

    if(estimand == "ATT") {
        wts <- T - (1-T) * e / (1-e) 
    } else if (estimand == "ATC") {
        wts <- T * (1-e) / e - (1-T)
    } else {
        wts <- (T / e) - ((1 - T) / (1 - e))
    }

  if(hajek) {
    ## Note that negative weights for control assignments are flipped to positive by normalization
    sum(wts[T == 1] * Y[T == 1]) / sum(wts[T == 1]) - sum(wts[T == 0] * Y[T == 0]) / sum(wts[T == 0])
  } else {
    mean(wts * Y)
  }
}
