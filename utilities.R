shrinkage_plot <- function(a, b) {
    plot(x=a, y=rep(1, length(a)), ylim=c(-0.1, 1.1))
    points(x=b, y=rep(0, length(b)))
    segments(x0=a, x1=b, y0=rep(1, length(a)), y1=rep(0, length(b)))
    abline(v=1/2, lty=2)
    abline(h=1)
    abline(h=0)
}

get_bias <- function(T, Y, X, xb, mvecs, mvals, ab_dot_prod, w2=0, w2lim, times=10) {
    
    N <- NullC(mvecs)
    if(abs(w2) == abs(w2lim)) {

        w1 <- sqrt(as.numeric((ab_dot_prod - mvals[2] * w2^2)/mvals[1]))
        g <- w1*mvecs[, 1] + w2 * mvecs[, 2]
        d <- Re(X %*% g)
        
        ## res <- glm(T ~ d, family="binomial")
        ## e_old <- predict(res, type="response")
        
        c1 <- as.numeric(t(xb/sqrt(sum(xb^2))) %*% d/sqrt(sum(d^2)))
        c2 <- sqrt(1-c1^2)
        normal_samples <- rnorm(10000, sd=1)
        e_d <- sapply(d, function(di) {
            mean(invlogit(escale*(c1*di + c2*normal_samples), a=0, b=1))
        })
                
        bias0 <- sum(1 / (1-e_d[T==0]) * Y[T==0]) / sum(1 / (1-e_d[T==0]))
        bias1 <- sum(1 / e_d[T==1] * Y[T==1]) / sum(1 / e_d[T==1])
        times <- 1
    }
    else  {

        w1 <- sqrt(as.numeric((ab_dot_prod - mvals[2] * w2^2)/mvals[1]))
        bias0 <- bias1 <- 0
        for(i in 1:times) {
            g <- w1*mvecs[, 1] + w2 * mvecs[, 2] + sqrt(1-w1^2 - w2^2) * N %*% rustiefel(p-2, 1)
            d <- Re(X %*% g)

            c1 <- as.numeric(t(xb/sqrt(sum(xb^2))) %*% d/sqrt(sum(d^2)))
            c2 <- sqrt(1-c1^2)
            normal_samples <- rnorm(10000, sd=1)
            e_d <- sapply(d, function(di) {
                mean(invlogit(escale*(c1*di + c2*normal_samples), a=0, b=1))
            })
            
            ## res <- glm(T ~ d - 1, family="binomial")
            ## e_d <- predict(res, type="response")
            
            bias0 <- bias0 + sum(1 / (1-e_d[T==0]) * Y[T==0]) / sum(1 / (1-e_d[T==0]))
            bias1 <- bias1 + sum(1 / e_d[T==1] * Y[T==1]) / sum(1 / e_d[T==1])

            if(is.na(bias0) | is.na(bias1))
                browser()
            
        }
    }
    list(bias0=bias0/times, bias1=bias1/times)
}

invlogit <- function(x, a, b) {
    exp(a+x*b) / (1+exp(a+x*b))
}

