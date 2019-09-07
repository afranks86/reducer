shrinkage_plot <- function(a, b) {
    plot(x=a, y=rep(1, length(a)), ylim=c(-0.1, 1.1))
    points(x=b, y=rep(0, length(b)))
    segments(x0=a, x1=b, y0=rep(1, length(a)), y1=rep(0, length(b)))
    abline(v=1/2, lty=2)
    abline(h=1)
    abline(h=0)
}

get_bias <- function(T, Y, X, xb, mvecs, mvals, ab_dot_prod, w2=0, w2lim, times=10, DEBUG=FALSE) {
    
    N <- NullC(mvecs)

    w1 <- sqrt(as.numeric((ab_dot_prod - mvals[2] * w2^2)/mvals[1]))
    bias0 <- bias1 <- 0
    normal_samples <- rnorm(500, sd=1)
                                        #biases <- rep(NA, times)

    for(i in 1:times) {
        # Stop if magnitudes of w1 and w2 are too large, and likely not due to numerical error.
        # Inserted so that max statement in next line does not silence bugs.
        stopifnot(1 - w1^2 - w2^2 > -1e-6)
        g <- w1 * mvecs[, 1] + w2 * mvecs[, 2] + sqrt(max(1 - w1^2 - w2^2, 0)) * N %*% rustiefel(p-2, 1)
        d <- Re(X %*% g)
        
        if(DEBUG){
          mhatX <- X %*% alpha_hat_normalized
          eX <- X %*% beta
          
          print(sprintf("Partial correlation: %s", round(cor(lm(mhatX ~ d)$resid, lm(eX ~ d)$resid), 4)))
          
          print(sprintf("Hyperbola condition: %s",
                        round((alpha_hat_normalized %*% as.vector(g)) *
                                (as.vector(g) %*% beta) - (alpha_hat_normalized %*% beta), 4)))
          
          if(w2 == 0)
            print(sprintf("Checking that w2 = 0 correspond to bisector: %s",
                          round(sum(alpha_hat_normalized %*% as.vector(g) - beta %*% as.vector(g)), 4)))
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

ipw_est <- function(e, T, Y, hajek=FALSE){
  wts <- (T / e) - ((1 - T) / (1 - e))
  if(hajek) {
    # Note that negative weights for control assignments are flipped to positive by normalization
    sum(wts[T == 1] * Y[T == 1]) / sum(wts[T == 1]) - sum(wts[T == 0] * Y[T == 0]) / sum(wts[T == 0])
  } else {
    mean(wts * Y)
  }
}
