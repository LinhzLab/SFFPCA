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
q <- 2
K <- 2
p <- 500 
tau_Bt <- 5
T <- 100

set.seed(143)
xyz <- matrix(0, p + 100, 3)
for (k in 1 : 3) xyz[ , k] <- sort(runif(p + 100, 0, 1), decreasing = FALSE) #xyz[ , k] <- seq(0, 1, length.out = p)
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


set.seed(70)
Alpha_true <-  mvrnorm(q, rep(0, tau_Bs), diag(tau_Bs)) 
Alpha_true <- t((svd(Alpha_true))$v)
f_true <- Bs %*% t(Alpha_true)
B_true <- matrix(0, p, q)
B_true[1 : 100, ] <- t(matrix(c(1,1)/10, q, 100)) * sqrt(p/4)
B_true <- B_true[, c(1,2)]
Alpha_true <- Alpha_true[c(1,2), ]
f_true <- Bs %*% t(Alpha_true)


set.seed(2)
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


n <- 400
set.seed(3)
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


BF_result <- FF_result <- PP_result <- SS_result <- NNMI <- rep(0, 47)


for (seed in  1 : 47) {
  set.seed(seed)
  for (t in 1 : T) X[, ,t] <- HBf_true[, ,t]  + mvrnorm(n = n, rep(0, p), e_var)
  XX <- AB_array_multi(X, X)
  XX_svd <- svd(XX)
  
  Bf0 <- XX_svd$u[, 1 : q] * sqrt(p) 
  for (k in 1 : q) {
    index_k <- which( abs((B_true+f_true)[, k]) == max( abs((B_true+f_true)[, k]) ) )
    Bf0[, k] <- Bf0[, k] * sign(Bf0[index_k, k] * (B_true+f_true)[index_k, k])
  }
  
  
  set.seed(seed)
  B0 <-  B_true + mvrnorm(p, rep(0, q), 0.3 * diag(q) )
  f0 <- Bf0 - B0
  
  W <- array(0, c(n, tau_Bt, q) )
  phi <- array(0,c(q, K, T) )
  Theta_q <- array(0, c(K, tau_Bt, q) )
  score_n <- array(0, c(n, q, K) )
  B <- B0
  f <- f0
  H <- array(0, c(n, q, T))
  for (t in 1 : T) H[, , t] = X[, , t] %*% (B+f) %*% solve( t(B+f) %*% (B+f) )
  
  lambda1 <- 0.12
  lambda2 <- 0.5
  D <- kronecker(diag(p) * (p) - matrix(rep(1, p*p) ,p ,p), diag(q) )
  E <- rep(0, p*q)
  mu <- 0.001
  B_minus_B <-  array(0, c(p, p, q))
  for (j in 1 : p ) B_minus_B[, j, ]  <- t(apply(B, 1, function(x) x - B[j, ] ))  
  Gam <- B_minus_B
  nu <- array(0, c(p, p, q)) 
  
  for (kstep in 1 : 3) {
    wtH <- matrix(0, p*q, p*q)
    for (j in 1 : p) wtH[ ((j-1)*q+1) : (j*q) , ((j-1)*q+1) : (j*q) ] <- AB_array_multi(H, H) 
    X_minus_fH <- X_minus_AH(X, f, H)
    XH_fH <- AB_array_multi(X_minus_fH, H) 
    wtH_X_minus_fH  <- c(t(XH_fH))
    HH_inv <- solve( AB_array_multi(H, H) )
    
    for (lll in 1 : 1) {
      for (j in 1 : p) E[ ((j-1)*q+1)  :  (j*q) ] <- colSums(nu[j, -j,] - nu[-j, j, ] - Gam[j, -j, ] + Gam[-j, j, ]) 
      wtB <- solve(2 * wtH + mu * D) %*% (2 * wtH_X_minus_fH - E )
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

    for (t in 1 : T) H[, , t] = X[, , t] %*% (B+f) %*% solve( t(B+f) %*% (B+f) )
    for (i in 1 : n) W[i, , ] =  BBt_inv %*% t(Bt) %*% t(H[i, , ])
    for (l in 1 : q) {
      Theta_score_result_l <- Theta_score_ite(W[, , l], K)
      Theta_q[, , l] <- Theta_score_result_l[[1]]     
      score_n[, l, ] <- Theta_score_result_l[[2]]     
      for (kk in 1 : k) if( ( score_n[1, l, kk] * score_n_true[1, l, kk] ) < 0) {score_n[, l, kk] <- -score_n[, l, kk]}
      phi[l, , ] <- t( Bt %*% t( Theta_q[, , l] ))
      for (kk in 1 : K) if( ( phi[l, kk, 1] * phi_true[l, kk, 1] ) < 0 )  {Theta_q[kk, , l] <- -Theta_q[kk, , l]}
      phi[l, , ] <- t( Bt %*% t( Theta_q[, , l] ))
      H[, l, ] <-  score_n[, l, ] %*% phi[l, , ]
    }
    Alpha <- Alpha_ite(X, B, H, HH_inv, Bs, BBs_inv)     
    f <- Bs %*% t(Alpha)
    for (k in 1 : q) {
      index_k <- which( abs((f_true)[, k]) == max( abs((f_true)[, k]) ) )
      Alpha[k, ] <- Alpha[k, ] * sign( (f_true)[index_k, k] * (f)[index_k, k] )
    }
    f <- Bs %*% t(Alpha)
  }
  Alpha2 = Alpha
  B2 = B
  f2 = f 
  phi2 = phi
  score_n2 = score_n
  Theta_q2 = Theta_q
  H2 = H
  
  
  
  
  BF_result[seed] <- mean( (B2+f2-B_true-f_true)^2 )
  FF_result[seed] <- mean( (f2-f_true)^2 )
  PP_result[seed] <- mean( (phi2-phi_true)^2 )
  SS_result[seed] <- mean((score_n2-score_n_true)^2 )
  
  BB_select_true <- c( rep(1,100), rep(2,400) )
  BB_select <-  rep(0, p)
  BB_select[which(B[,1] > 0.5 )] <- 1
  BB_select[which(B[,1] < 0.2 )] <- 2
  
  CD_11 <- intersect(which(BB_select==1), which(BB_select_true==1))
  CD_12 <- intersect(which(BB_select==1), which(BB_select_true==2))
  CD_21 <- intersect(which(BB_select==2), which(BB_select_true==1))
  CD_22 <- intersect(which(BB_select==2), which(BB_select_true==2))
  CD_31 <- intersect(which(BB_select==0), which(BB_select_true==1))
  CD_32 <- intersect(which(BB_select==0), which(BB_select_true==2))
  
  CD_11 <- length(CD_11)
  CD_12 <- length(CD_12)
  CD_21 <- length(CD_21)
  CD_22 <- length(CD_22)
  CD_31 <- length(CD_31)
  CD_32 <- length(CD_32)
  
  I11 <- CD_11/p * log( p * CD_11 / length(which(BB_select==1)) / length(which(BB_select_true==1)) )
  I12 <- CD_12/p * log( p * CD_12 / length(which(BB_select==1)) / length(which(BB_select_true==2)) )
  I21 <- CD_21/p * log( p * CD_21 / length(which(BB_select==2)) / length(which(BB_select_true==1)) )
  I22 <- CD_22/p * log( p * CD_22 / length(which(BB_select==2)) / length(which(BB_select_true==2)) )
  I31 <- CD_31/p * log( p * CD_31 / length(which(BB_select==0)) / length(which(BB_select_true==1)) )
  I32 <- CD_32/p * log( p * CD_32 / length(which(BB_select==0)) / length(which(BB_select_true==2)) )
  
  if( CD_11 == 0  ) I11 <- 0
  if( CD_12 == 0  ) I12 <- 0
  if( CD_21 == 0  ) I21 <- 0
  if( CD_22 == 0  ) I22 <- 0
  if( CD_31 == 0  ) I31 <- 0
  if( CD_32 == 0  ) I32 <- 0
  
  ICD <- I11 + I12 + I21 + I22 + I31 + I32
  
  HC1 <-  -length( which(BB_select==1)) / p  * log( length( which(BB_select==1)) / p  )
  HC2 <-  -length( which(BB_select==2)) / p  * log( length( which(BB_select==2)) / p  )
  HC3 <-  -length( which(BB_select==0)) / p  * log( length( which(BB_select==0)) / p  )
  HD1 <-  -length(which(BB_select_true==1))     / p  * log( length(which(BB_select_true==1))     / p  )
  HD2 <-  -length(which(BB_select_true==2))     / p  * log( length(which(BB_select_true==2))     / p  )
  
  if (length( which(BB_select==1)) == 0) HC1 <- 0 
  if (length( which(BB_select==2)) == 0) HC2 <- 0 
  if (length( which(BB_select==0)) == 0) HC3 <- 0 
  
  HC <- HC1 + HC2  + HC3
  HD <- HD1 + HD2
  NMI <- 2 * ICD / (HD + HC)
  NNMI[seed] <- NMI
  
}

print(c(mean(NNMI), sd(NNMI)))
print(c(mean(FF_result),sd(FF_result)))
print(c(mean(PP_result),sd(PP_result)))
print(c(mean(SS_result),sd(SS_result)))


NMI_400 <- NNMI
FF_400 <- FF_result
PP_400 <- PP_result
SS_400 <- SS_result