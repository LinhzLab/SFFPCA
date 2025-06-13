X_minus_AH <- function(X, A, H) { 
  X_minus_AH <- array(0, c(n, p, T) )
  for (t in 1 : T) X_minus_AH[, , t] <- X[, , t] - H[, , t] %*% t(A)
  return(X_minus_AH)
}

AB_array_multi <- function(A, B) {
  AB <- matrix(0, length(A[1, , 1]), length(B[1, , 1]) )
  for (t in 1 : T) AB <- AB + t(A[, , t]) %*% B[, , t]
  return(AB)
}


Alpha_ite <- function(X, B, H, HH_inv, Bs, BBs_inv) {
  X_minus_BH <- X_minus_AH(X, B, H)
  XH_BH <- AB_array_multi(X_minus_BH, H) 
  Alpha <- HH_inv %*% t(XH_BH) %*% Bs %*% BBs_inv
  return(Alpha)
}

Theta_score_ite <- function(Wl, K) {
  Wl_svd <- svd( t(Wl) %*% Wl )
  Theta_q_l <- t( Wl_svd$u[, 1 : K] ) 
  score_n_l <- Wl %*% t(Theta_q_l)
  pro_l <- sum( (Wl_svd$d)[1 : K] ) / sum( (Wl_svd$d) )
  return(list(Theta_q_l, score_n_l, pro_l) )
}

library("fda")
library("MASS")
library("Matrix")
generation <- function(qq) {
  q <- qq
  K <- 2
  n <- 100
  p <- 500 
  tau_Bt <- 5
  T <- 100
  
  set.seed(1)
  xyz <- matrix(0, p + 100, 3)
  for (k in 1 : 3) xyz[ , k] <- sort(runif(p + 100, 0, 1), decreasing = FALSE) 
  tau_Bs_k <- 4
  break_num_s <- 4
  tau_Bs <- tau_Bs_k^3 
  Bs1 <- bsplineS(x = xyz[, 1],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s) 
  Bs2 <- bsplineS(x = xyz[, 2],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s) 
  Bs3 <- bsplineS(x = xyz[, 3],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s) 
  Bs_all <- matrix(t( KhatriRao( KhatriRao(t(Bs1), t(Bs2)), t(Bs3) ) ) , p + 100, tau_Bs ) * sqrt(p/4)
  
  
  set.seed(1)
  test_index <- floor(seq(1,600,length.out = 100))
  train_index <- (1:600)[-test_index]
  xyz_test <- xyz[test_index, ]
  xyz <- xyz[train_index, ]
  Bs <- Bs_all[train_index, ]
  BBs_inv <- solve( t(Bs) %*% Bs + 0.1 * diag(tau_Bs)) 
  Weight_s <- matrix(0, p, p)
  for (j in 1 : p) Weight_s[j, -j] <- apply(xyz[-j, ], 1, function(x) return(1/norm(x - xyz[j, ], type="2")) )
  Weight_s <- t(Weight_s)
  
  
  set.seed(1)
  Alpha_true <-  mvrnorm(q, rep(0, tau_Bs), diag(tau_Bs)) 
  Alpha_true <- t((svd(Alpha_true))$v)
  f_true <- Bs %*% t(Alpha_true)
  B_true <- matrix(0, p, q)
  B_true[1 : 100, ] <- t(matrix(c(1,1)/10, q, 100)) * sqrt(p/4)
  f_true <- Bs %*% t(Alpha_true)
  
  
  set.seed(1)
  time <- runif(T, 0, 1)
  phi_true <- array(0, c(q, K, T))
  for (l in 1 : q)  {
    for (i in 1 : (K/2)) {
      phi_true[l, 2 * i - 1, ] <- sin(2 * i * pi * time) * sqrt(2)
      phi_true[l, 2 * i, ] <- cos(2 * i * pi * time) * sqrt(2)
    }
  }
  tau_Bt <- 5
  break_num_t <- 5 
  Bt <- bsplineS(time,breaks = seq(0, 1, length = break_num_t), norder = tau_Bt + 2 - break_num_t)
  Bt <- (svd( Bt %*% t(Bt) ))$u[, 1 : tau_Bt] * sqrt(T)
  for (k in 1 : tau_Bt) Bt[, k] <- sign(Bt[1, k]) * Bt[, k]
  BBt_inv <- solve( t(Bt) %*% Bt )
  
  
  set.seed(1)
  score_var <- rep(1, K * q)
  for (j in 1 : q) score_var[ (j-1) * K + 1 ] <- 3 * (q + 1 - j)
  score_var <- diag(score_var)
  
  score_true <- mvrnorm(n = n, rep(0, K * q), score_var)
  for (k in 1 : (K*q) ) score_true[, k] <- sign(score_true[1, k]) * score_true[, k]
  score_n_true <- array(0, c(n, q, K) )
  for (i in 1 : n) score_n_true[i, ,] <- t( matrix(score_true[i, ], K, q) )
  H_true <- array(0, c(n, q, T))
  for (l in 1 : q) H_true[, l, ] <-  score_n_true[, l, ] %*% phi_true[l, , ]
  
  HBf_true <- X <- array(0, c(n, p, T))
  for (t in 1 : T) HBf_true[, , t] <- H_true[, ,t] %*% t(B_true+f_true)
  return(HBf_true)
}

p <- 500
n <- 100
T <- 100
e_var <- diag(1, p, p)
Q2 <- Q3 <- Q4 <- Q5 <- Q6 <- rep(0, 100)
Q2P <- Q3P <- Q4P <- Q5P <- Q6P <- rep(0, 100)
# q = 2 
HBf_true <- generation(2)
X <- array(0, c(n, p, T))
for (seed in 1 : 100) {
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  eigen_random <- rep(0, p)
  for (kk in 1 : 50) {
    set.seed(kk)
    X1 <- array( rnorm(n*p*T, 0, sd(X)), c(n,p,T)  )
    XX1 <- AB_array_multi(X1, X1)
    XX_svd1 <- svd(XX1)
    eigen_random <- eigen_random + XX_svd1$d / 50
  }
  q <- max(which( (XX_svd$d-eigen_random)>0  ))
  Q2P[seed] <- sum(XX_svd$d[1:q]) / sum(XX_svd$d)
  Q2[seed] <- q 
  #print(seed)
}  


# q = 3 
HBf_true <- generation(3)
X <- array(0, c(n, p, T))
for (seed in 1 : 200) {
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  eigen_random <- rep(0, p)
  for (kk in 1 : 50) {
    set.seed(kk)
    X1 <- array( rnorm(n*p*T, 0, sd(X)), c(n,p,T)  )
    XX1 <- AB_array_multi(X1, X1)
    XX_svd1 <- svd(XX1)
    eigen_random <- eigen_random + XX_svd1$d / 50
  }
  q <- max(which( (XX_svd$d-eigen_random)>0  ))
  Q3P[seed] <- sum(XX_svd$d[1:q]) / sum(XX_svd$d)
  Q3[seed] <- q 
  #print(seed)
}  


# q = 4 
HBf_true <- generation(4)
X <- array(0, c(n, p, T))
for (seed in 1 : 200) {
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  eigen_random <- rep(0, p)
  for (kk in 1 : 50) {
    set.seed(kk)
    X1 <- array( rnorm(n*p*T, 0, sd(X)), c(n,p,T)  )
    XX1 <- AB_array_multi(X1, X1)
    XX_svd1 <- svd(XX1)
    eigen_random <- eigen_random + XX_svd1$d / 50
  }
  q <- max(which( (XX_svd$d-eigen_random)>0  ))
  Q4P[seed] <- sum(XX_svd$d[1:q]) / sum(XX_svd$d)
  Q4[seed] <- q 
  #print(seed)
}  


# q = 5 
HBf_true <- generation(5)
X <- array(0, c(n, p, T))
for (seed in 1 : 200) {
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  eigen_random <- rep(0, p)
  for (kk in 1 : 50) {
    set.seed(kk)
    X1 <- array( rnorm(n*p*T, 0, sd(X)), c(n,p,T)  )
    XX1 <- AB_array_multi(X1, X1)
    XX_svd1 <- svd(XX1)
    eigen_random <- eigen_random + XX_svd1$d / 50
  }
  q <- max(which( (XX_svd$d-eigen_random)>0  ))
  Q5P[seed] <- sum(XX_svd$d[1:q]) / sum(XX_svd$d)
  Q5[seed] <- q 
  #print(seed)
}  


# q = 6 
HBf_true <- generation(6)
X <- array(0, c(n, p, T))
for (seed in 1 : 200) {
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  eigen_random <- rep(0, p)
  for (kk in 1 : 50) {
    set.seed(kk)
    X1 <- array( rnorm(n*p*T, mean(X), sd(X)), c(n,p,T)  )
    XX1 <- AB_array_multi(X1, X1)
    XX_svd1 <- svd(XX1)
    eigen_random <- XX_svd1$d / 50
  }
  q <- max(which( (XX_svd$d-eigen_random)>0  ))
  Q6P[seed] <- sum(XX_svd$d[1:q]) / sum(XX_svd$d)
  Q6[seed] <- q 
  #print(seed)
}  


print(c("Q=2",mean(Q2)))
print(c("Q=3",mean(Q3)))
print(c("Q=4",mean(Q4)))
print(c("Q=5",mean(Q5)))
print(c("Q=6",mean(Q6)))