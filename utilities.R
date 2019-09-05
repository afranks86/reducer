shrinkage_plot <- function(a, b) {
    plot(x=a, y=rep(1, length(a)), ylim=c(-0.1, 1.1))
    points(x=b, y=rep(0, length(b)))
    segments(x0=a, x1=b, y0=rep(1, length(a)), y1=rep(0, length(b)))
    abline(v=1/2, lty=2)
    abline(h=1)
    abline(h=0)
}

get_bias <- function(T, residual, alpha_hat_normalized, beta, tau_hat, w2=0, times=10) {
    
    M <- (alpha_hat_normalized %*% t(beta) + beta %*% t(alpha_hat_normalized))/2
    eig <- eigen(M)
    mvecs <- eig$vectors[, c(1,p)]
    mvals <- eig$values[c(1,p)]
    if(abs(mvals[2]) > mvals[1]) {
        mvals <- -mvals[2:1]
        mvecs <- mvecs[, 2:1]
    }
    
    const <- abs(t(alpha_hat_normalized) %*% beta)
    # print(const)

    w2_beta_lim <- (t(beta) %*% mvecs)[2]
    
    N <- NullC(mvecs)
    if(abs(w2) == abs(w2_lim)) {
        g <- w1*mvecs[, 1] + w2 * mvecs[, 2]
        d <- Re(X %*% g)
        res <- glm(T ~ d, family="binomial")
        e_d <- predict(res, type="response")
                
        bias0 <- sum(1 / (1-e_d[T==0]) * residual[T==0]) / sum(1 / (1-e_d[T==0]))
        bias1 <- sum(1 / e_d[T==1] * residual[T==1]) / sum(1 / e_d[T==1])
        times <- 1
    }
    else  {
        g <- w1*mvecs[, 1] + w2 * mvecs[, 2] + sqrt(1-w1^2 - w2^2) * N %*% rustiefel(p-2, 1)
        
        
        w1 <- sqrt(as.numeric((const - mvals[2] * w2^2)/mvals[1]))
        NC <- NullC(mvecs)
        
        bias0 <- bias1 <- 0
        for(i in 1:times) {
            g <- w1*mvecs[, 1] + w2 * mvecs[, 2] + sqrt(1-w1^2 - w2^2) * NC %*% rustiefel(p-2, 1)
            
            if(w2==0 & (abs(t(g) %*% beta) - abs(t(g) %*% alpha_hat_normalized) > 1e-10))
                browser()
            
            d <- Re(X %*% g)
            res <- glm(T ~ d, family="binomial")
            e_d <- predict(res, type="response")
            
            bias0 <- bias0 + sum(1 / (1-e_d[T==0]) * residual[T==0]) / sum(1 / (1-e_d[T==0]))
            bias1 <- bias1 + sum(1 / e_d[T==1] * residual[T==1]) / sum(1 / e_d[T==1])

            if(is.na(bias0) | is.na(bias1))
                browser()
            
        }
    }
    list(bias0=bias0/times, bias1=bias1/times)
}
