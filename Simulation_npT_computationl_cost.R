library("fda")
library("MASS")
library("Matrix")

SFFPCA <- function(n1,p1,T1,Kstep,seed) {
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
  
  q <- 2
  K <- 2
  n <- n1
  p <- p1 
  tau_Bt <- 5
  T <- T1
  
  set.seed(1)
  xyz <- matrix(0, p, 3)
  for (k in 1 : 3) xyz[ , k] <- sort(runif(p, 0, 1), decreasing = FALSE) 
  tau_Bs_k <- 4
  break_num_s <- 4
  tau_Bs <- tau_Bs_k^3 
  Bs1 <- bsplineS(x = xyz[, 1],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s) 
  Bs2 <- bsplineS(x = xyz[, 2],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s) 
  Bs3 <- bsplineS(x = xyz[, 3],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s) 
  Bs <- matrix(t( KhatriRao( KhatriRao(t(Bs1), t(Bs2)), t(Bs3) ) ) , p, tau_Bs ) * sqrt(p/4)
  BBs_inv <- solve( t(Bs) %*% Bs + 0.1 * diag(tau_Bs)) 
  Weight_s <- matrix(0, p, p)
  for (j in 1 : p) Weight_s[j, -j] <- apply(xyz[-j, ], 1, function(x) return(1/norm(x - xyz[j, ], type="2")) )
  Weight_s <- t(Weight_s)
  
  set.seed(1)
  Alpha_true <-  mvrnorm(q, rep(0, tau_Bs), diag(tau_Bs)) 
  Alpha_true <- t((svd(Alpha_true))$v)
  f_true <- Bs %*% t(Alpha_true)
  B_true <- matrix(0, p, q)
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
  score_var <- diag(c(6,1,3,1))
  score_true <- mvrnorm(n = n, rep(0, K * q), score_var)
  for (k in 1 : (K*q) ) score_true[, k] <- sign(score_true[1, k]) * score_true[, k]
  score_n_true <- array(0, c(n, q, K) )
  for (i in 1 : n) score_n_true[i, ,] <- t( matrix(score_true[i, ], K, q) )
  H_true <- array(0, c(n, q, T))
  for (l in 1 : q) H_true[, l, ] <-  score_n_true[, l, ] %*% phi_true[l, , ]
  
  
  HBf_true <- X <- array(0, c(n, p, T))
  e_var <- diag(3, p, p)
  for (t in 1 : T) HBf_true[, , t] <- H_true[, ,t] %*% t(B_true+f_true)
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  Bf0 <- XX_svd$u[, 1 : q] * sqrt(p) 
  H0 <- array(0, c(n, q, T))
  for (i in 1 : n) H0[i, , ] <- t(Bf0) %*% X[i, , ] / p 
  Alpha0 <- t( BBs_inv %*% t(Bs) %*% Bf0 )
  f0 <- Bs %*% t(Alpha0)
  B0 <- Bf0 - f0
  
  
  W <- array(0, c(n, tau_Bt, q) )
  phi <- array(0,c(q, K, T) )
  Theta_q <- array(0, c(K, tau_Bt, q) )
  score_n <- array(0, c(n, q, K) )
  B <- B0
  f <- Bf0 - B
  H <- H0
  lambda1 <- 0.01 
  lambda2 <- 0.3
  D <- kronecker(diag(p) * (p) - matrix(rep(1, p*p) ,p ,p), diag(q) )
  E <- rep(0, p*q)
  mu <- 0.001
  B_minus_B <-  array(0, c(p, p, q))
  for (j in 1 : p ) B_minus_B[, j, ]  <- t(apply(B, 1, function(x) x - B[j, ] ))  
  Gam <- B_minus_B
  nu <- array(0, c(p, p, q)) 
  loss <- 0
  
  aa <- system.time({
    for (kstep in 1 : Kstep) {
    loss1 <- loss
    wtH <- matrix(0, p*q, p*q)
    for (j in 1 : p) wtH[ ((j-1)*q+1) : (j*q) , ((j-1)*q+1) : (j*q) ] <- AB_array_multi(H, H) 
    X_minus_fH <- X_minus_AH(X, f, H)
    XH_fH <- AB_array_multi(X_minus_fH, H) 
    wtH_X_minus_fH  <- c(t(XH_fH))
    HH_inv <- solve( AB_array_multi(H, H) +0.01*diag(q))
    
    for (lll in 1 : 1) {
      for (j in 1 : p) E[ ((j-1)*q+1)  :  (j*q) ] <- colSums(nu[j, -j,] - nu[-j, j, ] - Gam[j, -j, ] + Gam[-j, j, ]) 
      wtB <- solve(2 * wtH + mu * D+0.01*diag(p*q)) %*% (2 * wtH_X_minus_fH - E )
      B <- t(matrix(wtB, q, p))
      B <- t(apply(B, 1, function(x) ifelse( (norm(x,type = "2")-lambda1) > 0, return((norm(x,type="2")-lambda1) * x / norm(x,type="2")), return(rep(0, q)) )) )
      for (j in 1 : p ) B_minus_B[, j, ] <- t(apply(B, 1, function(x) x - B[j, ] ))  
      Gam <- nu / mu + B_minus_B
      for (j in 1 : p) {
        for (j1 in 1 : p) { 
          if (norm(c(Gam[j, j1, ]),type = "2") - lambda2*Weight_s[j, j1] > 0) { Gam[j, j1, ] <- (norm(Gam[j, j1, ],type = "2")-lambda2*Weight_s[j, j1]) * Gam[j, j1, ] / norm(Gam[j, j1, ],type = "2") }
          else {Gam[j, j1, ] <- rep(0,q)}
        }
      }
      nu <- nu + mu * (B_minus_B - Gam)/10
    }
    Alpha <- Alpha_ite(X, B, H, HH_inv, Bs, BBs_inv)     
    f <- Bs %*% t(Alpha)
    for (t in 1 : T) H[, , t] = X[, , t] %*% (B+f) %*% solve( t(B+f) %*% (B+f) +0.01*diag(q))
    for (i in 1 : n) W[i, , ] =  BBt_inv %*% t(Bt) %*% t(H[i, , ])
    for (l in 1 : q) {
      Theta_score_result_l <- Theta_score_ite(W[, , l], K)
      Theta_q[, , l] <- Theta_score_result_l[[1]]     
      score_n[, l, ] <- Theta_score_result_l[[2]]     
      phi[l, , ] <- t( Bt %*% t( Theta_q[, , l] ))
      phi[l, , ] <- t( Bt %*% t( Theta_q[, , l] ))
      H[, l, ] <-  score_n[, l, ] %*% phi[l, , ]
    }
    f <- Bs %*% t(Alpha)
    X_predict <- array(0, c(n, p, T) )
    for (t in 1 : T) X_predict[, , t] <- H[, ,t] %*% t(B+f) 
    loss <- mean((X_predict - X)^2) 
    if (abs(loss-loss1)< 1e-4) {break}
  }
    })
  return(c(aa))
}



TIME_100_100_100 <- TIME_100_200_100 <- TIME_100_500_100 <- TIME_200_100_100 <- TIME_200_200_100 <- TIME_200_500_100 <- TIME_500_100_100 <- TIME_500_200_100 <- TIME_500_500_100 <- matrix(0,10,5)
TIME_100_100_200 <- TIME_100_200_200 <- TIME_100_500_200 <- TIME_200_100_200 <- TIME_200_200_200 <- TIME_200_500_200 <- TIME_500_100_200 <- TIME_500_200_200 <- TIME_500_500_200 <- matrix(0,10,5)
TIME_100_100_500 <- TIME_100_200_500 <- TIME_100_500_500 <- TIME_200_100_500 <- TIME_200_200_500 <- TIME_200_500_500 <- TIME_500_100_500 <- TIME_500_200_500 <- TIME_500_500_500 <- matrix(0,10,5)

for (seed in 1 : 10) {
  set.seed(seed)
  TIME_100_100_100[seed,] <- SFFPCA(100, 100, 100,5,seed)
  #print( c(seed, 1))
  TIME_100_200_100[seed,] <- SFFPCA(100, 200, 100,5,seed)
  #print( c(seed, 2))
  TIME_100_500_100[seed,] <- SFFPCA(100, 500, 100,5,seed)
  #print( c(seed, 3))
  TIME_200_100_100[seed,] <- SFFPCA(200, 100, 100,5,seed)
  #print( c(seed, 4))
  TIME_200_200_100[seed,] <- SFFPCA(200, 200, 100,5,seed)
  #print( c(seed, 5))
  TIME_200_500_100[seed,] <- SFFPCA(200, 500, 100,5,seed)
  #print( c(seed, 6))
  TIME_500_100_100[seed,] <- SFFPCA(500, 100, 100,5,seed)
  #print( c(seed, 7))
  TIME_500_200_100[seed,] <- SFFPCA(500, 200, 100,5,seed)
  #print( c(seed, 8))
  TIME_500_500_100[seed,] <- SFFPCA(500, 500, 100,5,seed)
  #print( c(seed, 9))
  
  TIME_100_100_200[seed,] <- SFFPCA(100, 100, 200,5,seed)
  #print( c(seed, 11))
  TIME_100_200_200[seed,] <- SFFPCA(100, 200, 200,5,seed)
  #print( c(seed, 12))
  TIME_100_500_200[seed,] <- SFFPCA(100, 500, 200,5,seed)
  #print( c(seed, 13))
  TIME_200_100_200[seed,] <- SFFPCA(200, 100, 200,5,seed)
  #print( c(seed, 14))
  TIME_200_200_200[seed,] <- SFFPCA(200, 200, 200,5,seed)
  print( c(seed, 15))
  TIME_200_500_200[seed,] <- SFFPCA(200, 500, 200,5,seed)
  #print( c(seed, 16))
  TIME_500_100_200[seed,] <- SFFPCA(500, 100, 200,5,seed)
  #print( c(seed, 17))
  TIME_500_200_200[seed,] <- SFFPCA(500, 200, 200,5,seed)
  #print( c(seed, 18))
  TIME_500_500_200[seed,] <- SFFPCA(500, 500, 200,5,seed)
  #print( c(seed, 19))

  TIME_100_100_500[seed,] <- SFFPCA(100, 100, 500,5,seed)
  #print( c(seed, 21))
  TIME_100_200_500[seed,] <- SFFPCA(100, 200, 500,5,seed)
  #print( c(seed, 22))
  TIME_100_500_500[seed,] <- SFFPCA(100, 500, 500,5,seed)
  #print( c(seed, 23))
  TIME_200_100_500[seed,] <- SFFPCA(200, 100, 500,5,seed)
  #print( c(seed, 24))
  TIME_200_200_500[seed,] <- SFFPCA(200, 200, 500,5,seed)
  #print( c(seed, 25))
  TIME_200_500_500[seed,] <- SFFPCA(200, 500, 500,5,seed)
  #print( c(seed, 26))
  TIME_500_100_500[seed,] <- SFFPCA(500, 100, 500,5,seed)
  #print( c(seed, 27))
  TIME_500_200_500[seed,] <- SFFPCA(500, 200, 500,5,seed)
  #print( c(seed, 28))
  TIME_500_500_500[seed,] <- SFFPCA(500, 500, 500,5,seed)
  #print( c(seed, 29))
}

mean(TIME_100_100_100[,1])
mean(TIME_100_200_100[,1])
mean(TIME_100_500_100[,1])
mean(TIME_200_100_100[,1])
mean(TIME_200_200_100[,1])
mean(TIME_200_500_100[,1])
mean(TIME_500_100_100[,1])
mean(TIME_500_200_100[,1])
mean(TIME_500_500_100[,1])


mean(TIME_100_100_200[,1])
mean(TIME_100_200_200[,1])
mean(TIME_100_500_200[,1])
mean(TIME_200_100_200[,1])
mean(TIME_200_200_200[,1])
mean(TIME_200_500_200[,1])
mean(TIME_500_100_200[,1])
mean(TIME_500_200_200[,1])
mean(TIME_500_500_200[,1])

mean(TIME_100_100_500[,1])
mean(TIME_100_200_500[,1])
mean(TIME_100_500_500[,1])
mean(TIME_200_100_500[,1])
mean(TIME_200_200_500[,1])
mean(TIME_200_500_500[,1])
mean(TIME_500_100_500[,1])
mean(TIME_500_200_500[,1])
mean(TIME_500_500_500[,1])