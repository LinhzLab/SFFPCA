#current_path <- getwd()
#setwd(current_path)
#load("SimulationX1_mu.RData")

library("MASS")
library("fda")
library("Matrix")

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
BBs_inv <- solve( t(Bs) %*% Bs ) 

tau_Bt <- 5
break_num_t <- 5 
Bt <- bsplineS(time,breaks = seq(0, 1, length = break_num_t), norder = tau_Bt + 2 - break_num_t)
Bt <- (svd( Bt %*% t(Bt) ))$u[, 1 : tau_Bt] * sqrt(T)
for (k in 1 : tau_Bt) Bt[, k] <- sign(Bt[1, k]) * Bt[, k]
BBt_inv <- solve( t(Bt) %*% Bt ) 

tau_Bs_k <- 4
break_num_s <- 3
tau_Bs <- tau_Bs_k^3 
Bs1 <- bsplineS(x = xyz[, 1],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs2 <- bsplineS(x = xyz[, 2],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs3 <- bsplineS(x = xyz[, 3],seq(0, 1, length = break_num_s), tau_Bs_k + 2 - break_num_s)
Bs <- matrix(t( KhatriRao( KhatriRao(t(Bs1),t(Bs2)), t(Bs3) ) ) , p, tau_Bs )
Bs <- (svd( Bs %*% t(Bs) ))$u[, 1 : tau_Bs] * sqrt(p)
for (k in 1 : tau_Bs) Bs[, k] <- sign(Bs[1, k]) * Bs[, k]
BBs_inv <- solve( t(Bs) %*% Bs ) 

tau_Bt <- 5
break_num_t <- 5 
Bt <- bsplineS(time,breaks = seq(0, 1, length = break_num_t), norder = tau_Bt + 2 - break_num_t)
Bt <- (svd( Bt %*% t(Bt) ))$u[, 1 : tau_Bt] * sqrt(T)
for (k in 1 : tau_Bt) Bt[, k] <- sign(Bt[1, k]) * Bt[, k]
BBt_inv <- solve( t(Bt) %*% Bt ) 

q <- 4
K <- 5
step_max <- 10
X_SFFPCA <- X
X_mean <- mu_SFFPCA
for (i in 1 : n) X_SFFPCA[i, ,] <- X[i, , ] - X_mean 
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
lambda1 <- 10
lambda2 <- 0.1
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
  HH_inv <- solve( AB_array_multi(H, H) + 0.001*diag(q) )
  for (l in 1 : 1) {
    for (j in 1 : p) E[ ((j-1)*q+1)  :  (j*q) ] <- colSums(nu[j, -j,] - nu[-j, j, ] - Gam[j, -j, ] + Gam[-j, j, ]) 
    wtB <- solve(2 * wtH + mu * D/10 + 0.01*diag(p*q) ) %*% (2 * wtH_X_minus_fH - E/10 )
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
xi_SFFPCA <- matrix(0, n, K * q) 
for (i in 1 : n) xi_SFFPCA[i, ] <- c(score_n[i, , ])




library("glmnet")
g <- rep(0, n)
for (i in 1 : n) {
  for (j in 1 : (p-1)) g[i] <- g[i] + sum( (sin(time) + cos(time) ) * (X[i, j, ] * X[i, j+1, ]) / T)
  for (j in 1 : p )  g[i] <- g[i] + sum( ( sin(time) ) * ( X[i, j, ] ) / T)
  for (j in 1 : p )  g[i] <- g[i] + sum( ( cos(time) ) * ( X[i, j, ]^2 * sign(X[i, j, ]) ) / T)
  for (j in 1 : p )  g[i] <- g[i] + sum( ( cos(time) ) * ( X[i, j, ]^3 ) / T)
  
}

PE_SFFPCA <- rep(0, 200)

for (seed in 1 : 200) {
  
  set.seed(seed)
  Y <- g + rnorm(n, 0, 2)
  data_SFFPCA <- cbind(Y, xi_SFFPCA)
  
  train_index <- sample(1 : n, 0.5 * n, replace = FALSE)
  test_index <- c(1 : n)[-train_index]
  data_SFFPCA_train <- data_SFFPCA[train_index, ]
  data_SFFPCA_test <- data_SFFPCA[test_index, ]

  data_SFFPCA <- cbind(Y, xi_SFFPCA)
  data_SFFPCA_train <- data_SFFPCA[train_index, ]
  data_SFFPCA_test <-  data_SFFPCA[test_index, ]
  cv_ridge_SFFPCA <- cv.glmnet(data_SFFPCA_train[, -1], Y[train_index], alpha = 0,intercept = TRUE)
  best_lambda_ridge_SFFPCA <- cv_ridge_SFFPCA$lambda.min 
  ridge_SFFPCA <- glmnet(data_SFFPCA_train[, -1], Y[train_index], alpha = 0, intercept = TRUE, lambda = cv_ridge_SFFPCA$lambda.min  )
  pe_SFFPCA <- mean( sort((predict(ridge_SFFPCA, data_SFFPCA_test[, -1]) - Y[test_index]) ^ 2, decreasing = FALSE)[1:(0.95*length(test_index))] ) / var(Y[test_index])
  
  PE_SFFPCA[seed] <- pe_SFFPCA
}
print(c("PE_SFFPCA",mean(PE_SFFPCA), sd(PE_SFFPCA)))
