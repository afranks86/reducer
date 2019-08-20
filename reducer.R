estimate_d <- function(Y0, X, T, arm=0, angle=pi/4){
  beta_hat <- coef(glm(T ~ X - 1))
  alpha_hat <- coef(lm(Y0[T==arm] ~ X[T==arm, ] - 1))
  
  ## normalize to get d(X) surface
  alpha_hat_norm <- sqrt(sum(alpha_hat^2))
  beta_hat_norm <- sqrt(sum(beta_hat^2))
  alpha_hat_normalized <- alpha_hat / alpha_hat_norm
  beta_hat_normalized <- beta_hat / beta_hat_norm
  
  # null space of alpha_hat
  A <- diag(p) - alpha_hat_normalized %*% t(alpha_hat_normalized)
  # null space of beta_hat
  B <- diag(p) - beta_hat_normalized %*% t(beta_hat_normalized)
  
  AnotB  <- eigen(A %*% (diag(p) - B))$vectors[, 1]
  BnotA  <- eigen(B %*% (diag(p) - A))$vectors[, 1]
  
  gen_g <- function(angle, svec, Bside){
    if(Bside){
      z <- cbind(BnotA, eigen(A %*% (diag(p) - BnotA %*% t(BnotA)))$vectors[, svec]) %*%
        c(cos(angle), sin(angle))
      
      w <- eigen(A %*% (diag(p) - z %*% t(z)))$vectors[, 1:(p-2)]
    } else {
      z <- cbind(AnotB, eigen(B %*% (diag(p) - AnotB %*% t(AnotB)))$vectors[, svec]) %*%
        c(cos(angle), sin(angle))
      
      w <- eigen(B %*% (diag(p) - z %*% t(z)))$vectors[, 1:(p-2)]
    }
    
    u <- eigen(z %*% t(z) + w %*% t(w))$vectors[, 1:(p-1)]
    
    M <- u %*% t(u)
    
    
    ## gamma is now the eigenvector of I - M (step 4)
    eig <- eigen(diag(p) - M)
    
    Re(eig$vectors[, 1])
  }
  
  d <- X %*% gen_g(angle, 1, FALSE)
  print(cor(X %*% alpha_hat, d))
  print(cor(X %*% beta_hat, d))
  d
}

one_sided <- function(Y, X, T, arm, angle=pi/4, include_m=FALSE){
  d <- estimate_d(Y, X, T, arm, angle)
  if(!include_m || angle==0){
    matrix(d, nc=1)
  } else {
    m <- estimate_d(Y, X, T, arm, 0)
    cbind(d, m)
  }
}

two_sided <- function(Y, X, T, angle=pi/4, include_m=FALSE){
  #d1 <- estimate_d(Y, X, T, 1, angle)
  #d0 <- estimate_d(Y, X, T, 0, angle)
  #if(!include_m || angle==0){
  #  cbind(d1, d0)
  #} else {
  #  m1 <- estimate_d(Y, X, T, 1, 0)
  #  m0 <- estimate_d(Y, X, T, 0, 0)
  #  cbind(d1, d0, m1, m0)
  #}
  cbind(one_sided(Y, X, T, 0, angle, include_m),
        one_sided(Y, X, T, 1, angle, include_m))
}