#current_path <- getwd()
#parent_path <- dirname(current_path)
#setwd(parent_path)
library("MASS")
library("fda")
library("Matrix")

#load("SimulationX1_mu.RData")
###  Some numbers which need to be choosed
n <- 200
p <- 100
T <- 20


K0 <- 3
K1 <- K2 <- 4


set.seed(1)
time <- sort( runif(T, 0, 1) )
xyz <- matrix(0, p, 3)
for (k in 1 : 3) xyz[ , k] <- runif(p, 0, 1) 
HBf1_true <- array(0, c(n, p, T))
HBf2_true <- array(0, c(n, p, T))


score_1 <- matrix(0, n, K0) 
score_1[, 1] <- rnorm(n, 0, 12)
score_1[, 2] <- rnorm(n, 0, 7)
score_1[, 3] <- rnorm(n, 0, 4)

psi <- array(0, c(p, T, K0))
for (t in 1 : T) {
  for (j in 1 : p) {
    psi[j, t, 1] <- sqrt(2)/45 * sin( pi*(sum(xyz[j, ])/15) ) * cos(pi * t /18)
    psi[j, t, 2] <- sqrt(2)/45 * cos( pi*(sum(xyz[j, ])/15) ) * cos(pi * t /18)
    psi[j, t, 3] <- 1/(30*sqrt(3)) *  cos(pi * t /12)
  }
}

for (t in 1 : T) for (j in 1 : p) for (i in 1 : n) HBf1_true[i, j, t] <- score_1[i, ] %*% psi[j, t, ]


score_2 <- mvrnorm(n, rep(0, K1 * K2), diag( ((K1*K2) :1)/4 ) )
phis <- matrix(0, p, K1)
phit <- matrix(0, T, K2)

for (i in 1 : (K1/2) ) {
  phis[ ,2 * i - 1] <- sin(2 * i * pi * apply(xyz, 1, function(x) norm(x, type = "2")/sqrt(3)) ) /10 
  phis[ ,2 * i] <- cos(2 * i * pi * apply(xyz, 1, function(x) norm(x, type = "2")/sqrt(3)) ) /10
}


for (i in 1 : (K2/2) ) {
  phit[ ,2 * i - 1] <- sin(2 * i * pi * time) * sqrt(2)
  phit[ ,2 * i] <- cos(2 * i * pi * time) * sqrt(2)
}

for (t in 1 : T) for (j in 1 : p) for (i in 1 : n) HBf2_true[i, j, t] <- sum(matrix(score_2[i, ], K1, K2) * ( phis[j, ] %*% t(phit[t, ]) ) )

HBf_true <- HBf1_true + HBf2_true
e_var <- diag(0.2, p, p)
n <- n_test <- 100

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


tau_Bs_k <- 4
break_num_s <- 3
tau_Bs <- tau_Bs_k^3 
Bs1 <- bsplineS(x = xyz[, 1],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs2 <- bsplineS(x = xyz[, 2],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs3 <- bsplineS(x = xyz[, 3],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs <- matrix(t( KhatriRao( KhatriRao(t(Bs1),t(Bs2)), t(Bs3) ) ) , p, tau_Bs )
Bs <- (svd( Bs %*% t(Bs) ))$u[, 1 : tau_Bs] * sqrt(p)
for (k in 1 : tau_Bs) Bs[, k] <- sign(Bs[1, k]) * Bs[, k]
BBs_inv <- solve( t(Bs) %*% Bs ) # the last term in Alpha

tau_Bt <- 5
break_num_t <- 5 
Bt <- bsplineS(time,breaks = seq(0, 1, length = break_num_t), norder = tau_Bt + 2 - break_num_t)
Bt <- (svd( Bt %*% t(Bt) ))$u[, 1 : tau_Bt] * sqrt(T)
for (k in 1 : tau_Bt) Bt[, k] <- sign(Bt[1, k]) * Bt[, k]
BBt_inv <- solve( t(Bt) %*% Bt ) # the first term in W

tau_Bs_k <- 4
break_num_s <- 3
tau_Bs <- tau_Bs_k^3 
Bs1 <- bsplineS(x = xyz[, 1],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs2 <- bsplineS(x = xyz[, 2],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs3 <- bsplineS(x = xyz[, 3],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs <- matrix(t( KhatriRao( KhatriRao(t(Bs1),t(Bs2)), t(Bs3) ) ) , p, tau_Bs )
Bs <- (svd( Bs %*% t(Bs) ))$u[, 1 : tau_Bs] * sqrt(p)
for (k in 1 : tau_Bs) Bs[, k] <- sign(Bs[1, k]) * Bs[, k]
BBs_inv <- solve( t(Bs) %*% Bs ) # the last term in Alpha

tau_Bt <- 5
break_num_t <- 5 
Bt <- bsplineS(time,breaks = seq(0, 1, length = break_num_t), norder = tau_Bt + 2 - break_num_t)
Bt <- (svd( Bt %*% t(Bt) ))$u[, 1 : tau_Bt] * sqrt(T)
for (k in 1 : tau_Bt) Bt[, k] <- sign(Bt[1, k]) * Bt[, k]
BBt_inv <- solve( t(Bt) %*% Bt ) # the first term in W

PE <- rep(0,100)
q_select <- rep(0, 100)
n <- n_test <- 100
q <- 4
K <- 5
step_max <- 2
for (seed in 1 : 100) {
  
  set.seed(seed)
  X_all <- HBf_true
  for (t in 1 : T) X_all[, ,t] <- HBf_true[, ,t] + matrix(mvrnorm(n = 200, rep(0, p), e_var), 200, p, byrow = T) 
  index_train <- sort( sample(1 : 200, n, replace = FALSE), decreasing = FALSE)
  index_test <- (1 : 200)[-index_train]
  X <- X_all[index_train, , ]
  X_test <- X_all[index_test, , ]
  
  X_SFFPCA <- X
  X_test_SFFPCA <- X_test
  X_mean <- mu_SFFPCA
  for (i in 1 : n) X_SFFPCA[i, ,] <- X[i, , ] - X_mean 
  for (i in 1 : n_test) X_test_SFFPCA[i, ,] <- X_test[i, , ] - X_mean
  
  
  XX <- AB_array_multi(X_SFFPCA, X_SFFPCA)
  XX_svd <- svd(XX)
  Bf0 <- XX_svd$u[, 1 : q] * sqrt(p)     
  Alpha0 <- t( BBs_inv %*% t(Bs) %*% Bf0 )
  f0 <- Bs %*% t(Alpha0)
  B0 <- Bf0 - f0
  H0 <- array(0, c(n, q, T))
  for (i in 1 : n) H0[i, , ] <- t(Bf0) %*% X_SFFPCA[i, , ] / p 
  
  B <- B0
  f <- f0
  H <- H0
  lambda1 <- 1
  lambda2 <- 3
  Weight_s <- matrix(0, p, p)
  for (j in 1 : p) Weight_s[j, -j] <- apply(xyz[-j, ], 1, function(x) return(1/norm(x - xyz[j, ], type="2")) )
  Weight_s <- t(Weight_s)
  D <- kronecker(diag(p) * (p) - matrix(rep(1, p*p) ,p ,p), diag(q) )
  E <- rep(0, p*q)
  mu <- 0.001
  B_minus_B <-  array(0, c(p, p, q))
  for (j in 1 : p ) B_minus_B[, j, ]  <- t(apply(B, 1, function(x) x - B[j, ] ))  
  Gam <- B_minus_B
  nu <- array(0, c(p, p, q)) 
  W <- array(0, c(n, tau_Bt, q) )
  phi <- array(0,c(q, K, T) )
  Theta_q <- array(0, c(K, tau_Bt, q) )
  score_n <- array(0, c(n, q, K) )
  Pro <- matrix(0, q, step_max)
  for (k in 1 : step_max) {
    
    wtH <- matrix(0, p*q, p*q)
    for (j in 1 : p) wtH[ ((j-1)*q+1) : (j*q) , ((j-1)*q+1) : (j*q) ] <- AB_array_multi(H, H) 
    X_minus_fH <- X_minus_AH(X_SFFPCA, f, H)
    XH_fH <- AB_array_multi(X_minus_fH, H) 
    wtH_X_minus_fH  <- c(t(XH_fH))
    HH_inv <- solve( AB_array_multi(H, H) )
    for (l in 1 : 1) {
      for (j in 1 : p) E[ ((j-1)*q+1)  :  (j*q) ] <- colSums(nu[j, -j,] - nu[-j, j, ] - Gam[j, -j, ] + Gam[-j, j, ]) 
      wtB <- solve(2 * wtH + mu * D/10) %*% (2 * wtH_X_minus_fH - E/10 )
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
    
    Alpha <- Alpha_ite(X_SFFPCA, B, H, HH_inv, Bs, BBs_inv)    
    f <- Bs %*% t(Alpha)
    for (t in 1 : T) H[, , t] = X_SFFPCA[, , t] %*% (B+f) %*% solve( t(B+f) %*% (B+f) )
    for (i in 1 : n) W[i, , ] =  BBt_inv %*% t(Bt) %*% t(H[i, , ])
    for (l in 1 : q) {
      Theta_score_result_l <- Theta_score_ite(W[, , l], K)
      Theta_q[, , l] <- Theta_score_result_l[[1]]     
      score_n[, l, ] <- Theta_score_result_l[[2]]     
      phi[l, , ] <- t( Bt %*% t( Theta_q[, , l] ))
      H[, l, ] <-  score_n[, l, ] %*% phi[l, , ]
    }
  }
  
  H_test <- array(0, c(n_test, q, T) )
  W_test <- array(0, c(n_test, tau_Bt, q) )
  score_n_test <- array(0, c(n_test, q, K) )
  X_predict <- array(0, c(n_test, p, T) )
  for (t in 1 : T) H_test[, , t] = X_test_SFFPCA[, , t] %*% (B+f) %*% solve( t(B+f) %*% (B+f) )
  for (i in 1 : n_test) W_test[i, , ] =  BBt_inv %*% t(Bt) %*% t(H_test[i, , ])
  for (l in 1 : q) {
    score_n_test[, l, ] <- W_test[, , l] %*% t(Theta_q[, , l]) 
    H_test[, l, ] <-  score_n_test[, l, ] %*% phi[l, , ]
  }
  for (t in 1 : T) X_predict[, , t] <- H_test[, ,t] %*% t(B+f) 
  for (i in 1 : n_test) X_predict[i, ,] <- X_predict[i, , ] + X_mean
  pe <- mean(sort((X_predict - X_test)^2, decreasing = FALSE)[1 : (0.95 * length(X_predict))] ) / var(X_test) 
  PE[seed] <- pe
}

print(c("PE",mean(PE),sd(PE)))
