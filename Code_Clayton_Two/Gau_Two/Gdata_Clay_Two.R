# Fixed tau = 0.2; rho.0 = (1 / tau - 1) ^{-1} = 0.5
# Gdata_Clay_Two <- function(n, r, CR = c(4, 6), rho.0)
# {
#   set.seed(r)
#   #*********Survival data first*********************
#   #Clay Copula
#   U1 = runif(n); V1 = runif(n)
#   U2 = (U1 ^ (-rho.0) * V1 ^ (-rho.0/(rho.0+1)) - U1 ^ (-rho.0) + 1) ^ (-1/rho.0)
#   T.star = cbind(U1, U2)
#   
#   W <- as.matrix(runif(n, 0, 1))
#   b <- mvrnorm(n, rep(0, ncol(Sigma.b0)), Sigma.b0)
#   
#   B11 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B12 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B13 = runif(n, 0, 1.75)
#   B21 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B22 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B23 = runif(n, 0, 1.75)
#   
#   ##Solve the original data
#   B11.t = exp(W %*% gamma.011 + beta.01 * b[, 1] + gamma.012 * B11)
#   B12.t = exp(W %*% gamma.011 + beta.01 * b[, 1] + gamma.012 * B12)
#   
#   Temp.1 = log(1 + (k2 + b[, 2] * beta.01) * (-log(1 - T.star)[, 1]) / (k1 * B11.t)) / (k2 + b[, 2] * beta.01)
#   Temp.1[which(is.na(Temp.1) == 1)] = 0
#   T.orig.1 = Temp.1 * (Temp.1 <= B13)
#   
#   Temp.2 = log(((k2 + b[, 2] * beta.01) * (-log(1 - T.star)[, 1]) + k1 * B11.t + exp(
#     k2 * B13 + b[, 2] * beta.01 * B13) * (k1 * B12.t - k1 * B11.t)) / (k1 * B12.t)) / (k2 + b[, 2] * beta.01)
#   Temp.2[which(is.na(Temp.2) == 1)] = 0
#   T.orig.1 = T.orig.1 + Temp.2 * (Temp.2 > B13)
#   
#   
#   B21.t = exp(W %*% gamma.021 + beta.02 * b[, 1] + gamma.022 * B21)
#   B22.t = exp(W %*% gamma.021 + beta.02 * b[, 1] + gamma.022 * B22)
#   
#   Temp.1 = log(1 + (k2 + b[, 2] * beta.02) * (-log(1 - T.star)[, 2]) / (k1 * B21.t)) / (k2 + b[, 2] * beta.02)
#   Temp.1[which(is.na(Temp.1) == 1)] = 0
#   T.orig.2 = Temp.1 * (Temp.1 <= B23)
#   
#   Temp.2 = log(((k2 + b[, 2] * beta.02) * (-log(1 - T.star)[, 2]) + k1 * B21.t + exp(
#     k2 * B23 + b[, 2] * beta.02 * B23) * (k1 * B22.t - k1 * B21.t)) / (k1 * B22.t)) / (k2 + b[, 2] * beta.02)
#   Temp.2[which(is.na(Temp.2) == 1)] = 0
#   T.orig.2 = T.orig.2 + Temp.2 * (Temp.2 > B23)
#   
#   
#   T.orig.1[which(T.orig.1 == 'NaN')] = mean(T.orig.1, na.rm = TRUE)
#   T.orig.2[which(T.orig.2 == 'NaN')] = mean(T.orig.2, na.rm = TRUE)
#   T.orig <- matrix(c(T.orig.1, T.orig.2), n, 2)
#   
#   ##Censored survival data
#   C1 <- (runif(n, 0, CR[1]))
#   C2 <- (runif(n, 0, CR[2]))
#   C <- (matrix(c(C1, C2), n, 2))
#   
#   T.cen <- matrix(0, n, 2)
#   delta <- matrix(as.numeric(T.orig <= C), n, 2)
#   T.cen[T.orig[, 1] <= C1, 1] <- T.orig[, 1][T.orig[, 1] <= C1]
#   T.cen[T.orig[, 1] >= C1, 1] <- C1[T.orig[, 1] > C1]
#   T.cen[T.orig[, 2] <= C2, 2] <- T.orig[, 2][T.orig[, 2] <= C2]
#   T.cen[T.orig[, 2] > C2, 2] <- C2[T.orig[, 2] > C2]
#   T.max <- apply(T.cen, 1, max)
#   
#   t_cond = seq(0, 2, length = 20)
#   #***************Longitudinal data second*******************************
#   t <- lapply(1 : n, function(i) t_cond[t_cond < T.max[i]])
#   ni = unlist(lapply(1 : n, function(i) length(t[[i]])))
#   X3 <- (rbinom(n, 1, 0.5)); X4 <- (runif(n, 0, 2))
#   X <- cbind(1, unlist(t), rep(X3, ni), rep(X4, ni))
#   
#   
#   e <- lapply(1 : n, function(i) rnorm(ni[i], 0, sigma.0e))
#   y <- X %*% alpha.0 + rep(b[, 1], ni) + unlist(t) * rep(b[, 2], ni) + unlist(e)
#   
#   ##organize the longitudinal data and survival data into a table
#   obse.no <- unlist(lapply(1 : n, function(i) 1 : ni[i]))
#   id = rep(1 : n, ni)
#   
#   Mydata <- list(id = id, ni = ni, obse.no = obse.no, X = X, b.long1 = rep(
#     b[, 1], ni), b.long2 = rep(b[, 2], ni), rho = rep(1, length(id)), y = y,
#     T.cen = T.cen, W = W, b.surv = b, delta = delta, X_al = cbind(X3, X4), t = t, B1 = cbind(B11, B12, B13), B2 = cbind(B21, B22, B23))
#   return(Mydata)
#   
# }

Gdata_Clay_Two <- function(n, r, CR = c(4, 6), rho.0)
{
  set.seed(r)
  
  #********* Survival data first *********************
  # Clayton Copula: (U1, U2) ~ C_rho
  U1 <- runif(n); V1 <- runif(n)
  U2 <- (U1 ^ (-rho.0) * V1 ^ (-rho.0/(rho.0+1)) - U1 ^ (-rho.0) + 1) ^ (-1/rho.0)
  T.star <- cbind(U1, U2)
  
  # h_k = -log(U_k)  (IMPORTANT: use -log(U), NOT -log(1-U))
  h1 <- -log(T.star[, 1])
  h2 <- -log(T.star[, 2])
  
  W <- as.matrix(runif(n, 0, 1))
  b <- mvrnorm(n, rep(0, ncol(Sigma.b0)), Sigma.b0)
  
  B11 <- sample(c(-1,1), n, prob = c(0.5, 0.5), replace = TRUE)
  B12 <- sample(c(-1,1), n, prob = c(0.5, 0.5), replace = TRUE)
  B13 <- runif(n, 0, 1.75)
  
  B21 <- sample(c(-1,1), n, prob = c(0.5, 0.5), replace = TRUE)
  B22 <- sample(c(-1,1), n, prob = c(0.5, 0.5), replace = TRUE)
  B23 <- runif(n, 0, 1.75)
  
  ## -------- Solve original T1 --------
  B11.t <- exp(W %*% gamma.011 + beta.01 * b[, 1] + gamma.012 * B11)
  B12.t <- exp(W %*% gamma.011 + beta.01 * b[, 1] + gamma.012 * B12)
  
  d1 <- (k2 + b[, 2] * beta.01)
  
  Temp.1 <- log(1 + d1 * h1 / (k1 * B11.t)) / d1
  Temp.1[is.na(Temp.1) | is.nan(Temp.1) | is.infinite(Temp.1)] <- 0
  T.orig.1 <- Temp.1 * (Temp.1 <= B13)
  
  Temp.2 <- log((d1 * h1 + k1 * B11.t +
                   exp(k2 * B13 + b[, 2] * beta.01 * B13) * (k1 * B12.t - k1 * B11.t)) /
                  (k1 * B12.t)) / d1
  Temp.2[is.na(Temp.2) | is.nan(Temp.2) | is.infinite(Temp.2)] <- 0
  T.orig.1 <- T.orig.1 + Temp.2 * (Temp.2 > B13)
  
  ## -------- Solve original T2 --------
  B21.t <- exp(W %*% gamma.021 + beta.02 * b[, 1] + gamma.022 * B21)
  B22.t <- exp(W %*% gamma.021 + beta.02 * b[, 1] + gamma.022 * B22)
  
  d2 <- (k2 + b[, 2] * beta.02)
  
  Temp.1 <- log(1 + d2 * h2 / (k1 * B21.t)) / d2
  Temp.1[is.na(Temp.1) | is.nan(Temp.1) | is.infinite(Temp.1)] <- 0
  T.orig.2 <- Temp.1 * (Temp.1 <= B23)
  
  Temp.2 <- log((d2 * h2 + k1 * B21.t +
                   exp(k2 * B23 + b[, 2] * beta.02 * B23) * (k1 * B22.t - k1 * B21.t)) /
                  (k1 * B22.t)) / d2
  Temp.2[is.na(Temp.2) | is.nan(Temp.2) | is.infinite(Temp.2)] <- 0
  T.orig.2 <- T.orig.2 + Temp.2 * (Temp.2 > B23)
  
  ## final safety (avoid non-finite)
  bad1 <- is.na(T.orig.1) | is.nan(T.orig.1) | is.infinite(T.orig.1)
  bad2 <- is.na(T.orig.2) | is.nan(T.orig.2) | is.infinite(T.orig.2)
  if (any(bad1)) T.orig.1[bad1] <- mean(T.orig.1[!bad1], na.rm = TRUE)
  if (any(bad2)) T.orig.2[bad2] <- mean(T.orig.2[!bad2], na.rm = TRUE)
  
  T.orig <- cbind(T.orig.1, T.orig.2)
  
  ## -------- Censoring --------
  C1 <- runif(n, 0, CR[1])
  C2 <- runif(n, 0, CR[2])
  
  delta <- matrix(as.numeric(T.orig <= cbind(C1, C2)), n, 2)
  T.cen <- cbind(pmin(T.orig.1, C1), pmin(T.orig.2, C2))
  T.max <- apply(T.cen, 1, max)
  
  #*************** Longitudinal data second *******************************
  t_cond <- seq(0, 2, length = 20)
  t <- lapply(1:n, function(i) t_cond[t_cond < T.max[i]])
  ni <- vapply(t, length, integer(1))
  
  X3 <- rbinom(n, 1, 0.5)
  X4 <- runif(n, 0, 2)
  
  X <- cbind(1, unlist(t), rep(X3, ni), rep(X4, ni))
  e <- lapply(1:n, function(i) rnorm(ni[i], 0, sigma.0e))
  
  y <- X %*% alpha.0 + rep(b[, 1], ni) + unlist(t) * rep(b[, 2], ni) + unlist(e)
  
  obse.no <- unlist(lapply(1:n, function(i) seq_len(ni[i])))
  id <- rep(1:n, ni)
  
  Mydata <- list(
    id = id, ni = ni, obse.no = obse.no, X = X,
    b.long1 = rep(b[, 1], ni), b.long2 = rep(b[, 2], ni),
    rho = rep(1, length(id)), y = y,
    T.cen = T.cen, W = W, b.surv = b, delta = delta,
    X_al = cbind(X3, X4), t = t,
    B1 = cbind(B11, B12, B13), B2 = cbind(B21, B22, B23)
  )
  return(Mydata)
}

Cen_Clay_Two <- function(n, r, MR, rho.0)
{
  CR = seq(1, 3, length = 20)
  CRR = cbind(rep(CR, each = length(CR)), rep(CR, length(CR)))
  mr = matrix(0, nrow(CRR), 2)
  for(i in 1 : nrow(CRR))
  {
    dat = Gdata_Clay_Two(n, r, CR = CRR[i, ], rho.0) 
    mr[i, ] = 1 - apply(dat$delta, 2, mean)
  }
  i.opt = which.min(apply(abs(mr - MR), 1, mean))
  return(Gdata_Clay_Two(n, r, CR = c(CRR[i.opt, ]), rho.0))
}