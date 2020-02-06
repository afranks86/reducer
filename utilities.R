shrinkage_path <- function(tib, nms=NULL) {

    df <- tib
    df$observation = 1:nrow(df)

    df %>% gather(w, `Propensity Score`, -observation) %>%
        mutate(w = factor(w, levels=colnames(df))) %>%
        arrange(desc(w)) %>% 
        ggplot(aes(y = `Propensity Score`, x = w)) +
        geom_vline(xintercept=ncol(df)/2, col="red", size=2) +        
        geom_point() +
        geom_path(aes(group=observation)) +
        theme_bw(base_size=16) +
        scale_x_discrete(labels=nms)
    
}

double_shrinkage_plot <- function(baseline_e, e_top, e_bottom, nms=NULL) {


    df <- tibble(baseline=baseline_e, top=e_top, bottom=e_bottom)

    if(!is.null(nms))
        colnames(df) <- nms
    else
        nms <- colnames(df)
    
    df$observation = 1:nrow(df)

    df %>% gather(Type, `Propensity Score`, -observation) %>%
        mutate(Type = factor(Type, levels=nms[c(2, 1, 3)])) %>%
        arrange(desc(Type)) %>% 
        ggplot(aes(y = `Propensity Score`, x = Type)) +
        geom_point() +
        geom_path(aes(group=observation)) +
        theme_bw()
 
}


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


estimate_outcome <- function(X, T, Y, estimand, alpha=0,
                             include_intercept=TRUE,
                             scale_X = FALSE,
                             pred_policy = 'lambda.1se',
                             coef_policy = 'lambda.min') {
  if (scale_X) {
    scl = apply(X, 2, sd, na.rm = TRUE)
    is.binary = apply(X, 2, function(xx) sum(xx == 0) + sum(xx == 1) == length(xx))
    scl[is.binary] = 1
    Xscl = scale(X, center = FALSE, scale = scl)
  } else {
    Xscl = X
  }

  ## estimand <- "ATE"
  if (estimand == "ATT"){
    Xfit <- X[T==0, ]
    Yfit <- Y[T==0]
    Xpred <- X[T==1, ]

  } else {
    Xfit <- cbind(T, X)
    Yfit <- Y
    stop("Not supported in this branch")
  }

  cvglm <- glmnet::cv.glmnet(Xfit, Yfit, alpha = alpha,
                             intercept=include_intercept)
  balance_target  <-  colMeans(X[T==1,])
  lambdas <- list('lambda.min' = cvglm$lambda.min,
                  'lambda.1se' = cvglm$lambda.1se,
                  'undersmooth' = cvglm$lambda[max(which(cvglm$cvlo < min(cvglm$cvm)))])
  pred_lam <- lambdas[[pred_policy]]
  coef_lam <- lambdas[[coef_policy]]
  mu_pred <- predict(cvglm,
                     newx = matrix(balance_target, 1,
                                   length(balance_target)),
                     s=pred_lam)
  fitted_values <- predict(cvglm, newx=X, s=pred_lam)

  tau_hat  <- mean(Y[T==1]) - mu_pred

  intercept   <- coef(cvglm, s=coef_lam)[1]
  alpha_hat <- coef(cvglm, s=coef_lam)[-1]
  alpha_hat_normalized <- alpha_hat / sqrt(sum(alpha_hat^2))


  mhat0  <- mu_pred
  mhat1  <- mean(Y[T==1])

  list(alpha_hat=alpha_hat,
       alpha_hat_normalized=alpha_hat_normalized,
       tau_hat=tau_hat,
       intercept = intercept,
       mhat0=mhat0,
       mhat1=mhat1,
       preds=fitted_values)
}

estimate_propensity <- function(X, T, cv=TRUE, T_lambda_min=115,
                                eta=0.1, alpha=0){ 
    if(cv){
        cvglm <- cv.glmnet(cbind(X), T, family="binomial", 
                           alpha=alpha, penalty.factor = rep(1, p), intercept=FALSE)
        T_lambda_min <- cvglm$lambda.min
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


    propensity_fit_all <- glmnet(cbind(X), T, family="binomial", 
                             alpha=alpha, penalty.factor = rep(1, p),intercept=FALSE, 
                             lambda=cvglm$lambda)
    ehat_all <- predict(propensity_fit_all, type="response", newx=X)
    
    list(beta_hat=beta_hat,
         beta_hat_normalized=beta_hat_normalized,
         escale_hat=escale_hat,
         ehat=ehat,
         ehat_all = ehat_all,
         lambda=T_lambda_min)
}

## Grabs named attribute from list if not na; else returns default
get_attr_default <- function(thelist, attrname, default){
    if(!is.null(thelist[[attrname]])) thelist[[attrname]] else default
}

compute_gammas <- function(ab_dot_prod, mvals, mvecs, w2, times){
    p <- dim(mvecs)[1]

    ## Sample random basis vectors from the null space
    null_vecs <- NullC(mvecs)[, sample(p-ncol(mvecs), size=min(p-ncol(mvecs), times), replace=FALSE)]
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
    
    ## Compute many coefficient vectors gamma corresponding to reductions
    gammas <- compute_gammas(ab_dot_prod, mvals, mvecs, w2, times)

    ## Compute many corresponding reductions d
    dd <- X %*% gammas
    dd_norms <- sqrt(colSums(dd^2))
    dd_normalized <- t(t(dd) / dd_norms)

    ## Magnitudes of projected normalized linearized propensity onto dd
    ## There is a c1 and c2 for each d
    c1s <- t(xb / sqrt(sum(xb^2))) %*% dd_normalized

    ## Stop if magnitude of c1 is too large, and likely not due to numerical error.
    ## Inserted so that max statement in next line does not silence bugs.
    
    ## stopifnot(all(1 - c1s^2 >= -1e-8)) ## changed this to 0, why can' c1s be 0?
    c2s <- sqrt(pmax(1-c1s^2, 0))
    
    ## For each d, we will be averaging over a bunch of normal samples
    normal_samples <- rnorm(500, sd=1)
    ## Matrix of projected propensity scores
    e_dd <- matrix(NA, nr=dim(dd)[1], nc=dim(dd)[2])
    for(i in 1:ncol(e_dd)){
        ## For each component of dd, add each value of normal_sample and evaluate logistic
        ## Implement this as an outer sum between a column of dd and normal_samples
        integrand <- invlogit(escale * (outer(c1s[i] * dd[,i], c2s[i] * normal_samples, FUN="+")), a=0, b=1)
        e_dd[, i] <- rowMeans(integrand)
    }
   
    ## Compute extreme reduced propensity score values 
    eta <- min(mean(apply(e_dd, 2, min)), 1 - mean(apply(e_dd, 2, max)))
    
    ## Compute inverse weighted group means -- uses Hajek estimator
    ## First, comput weights for treated and control groups. There will be a vector of weights for each
    ## distinct reduction d.
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
   
    ## Compute vectors of weighted means for treated and control groups. 
    trt_wt_means <- colSums((trt_wt * as.vector(Y))[T == 1,,drop=FALSE]) /
        colSums(trt_wt[T == 1,,drop=FALSE])
    ctrl_wt_means <- colSums((ctrl_wt * as.vector(Y))[T == 0,,drop=FALSE]) /
        colSums(ctrl_wt[T == 0,,drop=FALSE])
    
    ## Return the means of these differently-weighted means. Since each is unbiased, the mean is also
    ## unbiased.
    list(bias0=mean(ctrl_wt_means), bias1=mean(trt_wt_means), eta=eta, e_dd=e_dd)
    
}

## TODO: Refactor this and get_bias. Have get_bias only return weights, then call get_bias.
compute_bias <- function(T, Y, e_dd, estimand="ATT"){
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

make_bias_var_plot  <- function() {



}
