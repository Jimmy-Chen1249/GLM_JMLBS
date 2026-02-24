# Fixed tau = 0.15; rho.0 = 9 * tau / 2 = 0.675
Gdata_FGM <- function(n, r, CR = c(4, 6), rho.0)
{
  set.seed(r)
  #*********Survival data first*********************
  #FGM Copula
  U1 = runif(n); V1 = runif(n)
  A1 = 1 + rho.0 * (1 - 2 * U1)
  B1 = sqrt(A1 ^ 2 - 4 * (A1 - 1) * V1)
  U2 = 2 * V1 / (B1 + A1)
  T.star = cbind(U1, U2)
  
  W <- as.matrix(runif(n, 0, 1))
  b <- rnorm(n, 0, D.0)
  
  B11 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B12 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B13 = runif(n, 0, 1.75)
  B21 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B22 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B23 = runif(n, 0, 1.75)
  ##Solve the original data
  
  Temp.1 = (log(1 + k2 * (-log(T.star)[, 1]) / (k1 * exp(
    W %*% gamma.011 + beta.01 * b + gamma.012 * B11))) / k2)
  Temp.1[which(is.na(Temp.1) == 1)] = 0
  T.orig.1 = Temp.1 * (Temp.1 <= B13)
  
  Temp.2 = (log((k2 * (-log(T.star)[, 1]) / (k1 * exp(W %*% gamma.011 + beta.01 * b)) - exp(k2 * B13) * (exp(
    B11 * gamma.012) - exp(B12 * gamma.012)) + exp(B11 * gamma.012)) / exp(B12 * gamma.012)) / k2)
  Temp.2[which(is.na(Temp.2) == 1)] = 0
  T.orig.1 = T.orig.1 + Temp.2 * (Temp.2 > B13)
  
  
  Temp.1 = (log(1 + k2 * (-log(T.star)[, 2]) / (k1 * exp(
    W %*% gamma.021 + beta.02 * b + gamma.022 * B21))) / k2)
  Temp.1[which(is.na(Temp.1) == 1)] = 0
  T.orig.2 = Temp.1 * (Temp.1 <= B23)
  
  Temp.2 = (log((k2 * (-log(T.star)[, 2]) / (k1 * exp(W %*% gamma.021 + beta.02 * b)) - exp(k2 * B23) * (exp(
    B21 * gamma.022) - exp(B22 * gamma.022)) + exp(B21 * gamma.022)) / exp(B22 * gamma.022)) / k2)
  Temp.2[which(is.na(Temp.2) == 1)] = 0
  T.orig.2 = T.orig.2 + Temp.2 * (Temp.2 > B23)
  
  T.orig.1[which(T.orig.1 == 'NaN')] = 10
  T.orig.2[which(T.orig.2 == 'NaN')] = 10
  T.orig <- matrix(c(T.orig.1, T.orig.2), n, 2)
  
  ##Censored survival data
  C1 <- (runif(n, 0, CR[1]))
  C2 <- (runif(n, 0, CR[2]))
  C <- (matrix(c(C1, C2), n, 2))
  
  T.cen <- matrix(0, n, 2)
  delta <- matrix(as.numeric(T.orig <= C), n, 2)
  T.cen[T.orig[, 1] <= C1, 1] <- T.orig[, 1][T.orig[, 1] <= C1]
  T.cen[T.orig[, 1] >= C1, 1] <- C1[T.orig[, 1] > C1]
  T.cen[T.orig[, 2] <= C2, 2] <- T.orig[, 2][T.orig[, 2] <= C2]
  T.cen[T.orig[, 2] > C2, 2] <- C2[T.orig[, 2] > C2]
  T.max <- apply(T.cen, 1, max)
  
  t_cond = seq(0, 2, length = 20)
  #***************Longitudinal data second*******************************
  t <- lapply(1 : n, function(i) t_cond[t_cond < T.max[i]])
  ni = unlist(lapply(1 : n, function(i) length(t[[i]])))
  X3 <- (rbinom(n, 1, 0.5)); X4 <- (runif(n, 0, 2))
  X <- cbind(1, unlist(t), rep(X3, ni), rep(X4, ni))
  
  
  e <- lapply(1 : n, function(i) rnorm(ni[i], 0, sigma.0e))
  y <- X %*% alpha.0 + rep(b, ni) + unlist(e)
  
  ##organize the longitudinal data and survival data into a table
  obse.no <- unlist(lapply(1 : n, function(i) 1 : ni[i]))
  id = rep(1 : n, ni)
  
  Mydata <- list(id = id, ni = ni, obse.no = obse.no, X = X, b.long = rep(
    b, ni), rho = rep(1, length(id)), y = y, T.cen = T.cen, W = W, b.surv = b, 
    delta = delta, X_al = cbind(X3, X4), t = t, B1 = cbind(B11, B12, B13), B2 = cbind(B21, B22, B23))
  return(Mydata)
  
}

Cen_FGM <- function(n, r, MR, rho.0)
{
  CR = seq(1, 3, length = 20)
  CRR = cbind(rep(CR, each = length(CR)), rep(CR, length(CR)))
  mr = matrix(0, nrow(CRR), 2)
  for(i in 1 : nrow(CRR))
  {
    dat = Gdata_FGM(n, r, CR = CRR[i, ], rho.0) 
    mr[i, ] = 1 - apply(dat$delta, 2, mean)
  }
  i.opt = which.min(apply(abs(mr - MR), 1, mean))
  return(Gdata_FGM(n, r, CR = c(CRR[i.opt, ]), rho.0))
}