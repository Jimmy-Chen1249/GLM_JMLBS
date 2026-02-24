#******************************************************************************************
Gdata = function(n, r, CR = c(4, 6))
{
  set.seed(r)
  
  #*********Survival data first**********************
  T.star <- mvrnorm(n, c(0, 0), matrix(c(1, rho.0, rho.0, 1), 2, 2))
  W <- as.matrix(runif(n, -2, 2))
  b <- mvrnorm(n, rep(0, ncol(Sigma.b0)), Sigma.b0)
  
  B11 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B12 = sample(c(
    -1,1), n, prob = c(0.5, 0.5), replace = T); B13 = runif(n, 0, 1.75)
  B21 = sample(c(-1,1), n, prob = c(0.5, 0.5), replace = T); B22 = sample(c(
    -1,1), n, prob = c(0.5, 0.5), replace = T); B23 = runif(n, 0, 1.75)
  
  ##Solve the original data
  B11.t = exp(W %*% gamma.011 + beta.01 * b[, 1] + gamma.012 * B11)
  B12.t = exp(W %*% gamma.011 + beta.01 * b[, 1] + gamma.012 * B12)
  
  Temp.1 = log(1 + (k2 + b[, 2] * beta.01) * (-log(1 - pnorm(T.star))[, 1]) / (k1 * B11.t)) / (k2 + b[, 2] * beta.01)
  Temp.1[which(is.na(Temp.1) == 1)] = 0
  T.orig.1 = Temp.1 * (Temp.1 <= B13)
  
  Temp.2 = log(((k2 + b[, 2] * beta.01) * (-log(1 - pnorm(T.star))[, 1]) + k1 * B11.t + exp(
    k2 * B13 + b[, 2] * beta.01 * B13) * (k1 * B12.t - k1 * B11.t)) / (k1 * B12.t)) / (k2 + b[, 2] * beta.01)
  Temp.2[which(is.na(Temp.2) == 1)] = 0
  T.orig.1 = T.orig.1 + Temp.2 * (Temp.2 > B13)
  
  
  B21.t = exp(W %*% gamma.021 + beta.02 * b[, 1] + gamma.022 * B21)
  B22.t = exp(W %*% gamma.021 + beta.02 * b[, 1] + gamma.022 * B22)
  
  Temp.1 = log(1 + (k2 + b[, 2] * beta.02) * (-log(1 - pnorm(T.star))[, 2]) / (k1 * B21.t)) / (k2 + b[, 2] * beta.02)
  Temp.1[which(is.na(Temp.1) == 1)] = 0
  T.orig.2 = Temp.1 * (Temp.1 <= B23)
  
  Temp.2 = log(((k2 + b[, 2] * beta.02) * (-log(1 - pnorm(T.star))[, 2]) + k1 * B21.t + exp(
    k2 * B23 + b[, 2] * beta.02 * B23) * (k1 * B22.t - k1 * B21.t)) / (k1 * B22.t)) / (k2 + b[, 2] * beta.02)
  Temp.2[which(is.na(Temp.2) == 1)] = 0
  T.orig.2 = T.orig.2 + Temp.2 * (Temp.2 > B23)
  
  
  T.orig.1[which(T.orig.1 == 'NaN')] = mean(T.orig.1, na.rm = TRUE)
  T.orig.2[which(T.orig.2 == 'NaN')] = mean(T.orig.2, na.rm = TRUE)
  T.orig <- matrix(c(T.orig.1, T.orig.2), n, 2)
  
  ##Censored survival data
  C1 = (runif(n, 0, CR[1]))
  C2 = (runif(n, 0, CR[2]))
  C = (matrix(c(C1, C2), n, 2))
  
  T.cen = matrix(0, n, 2)
  delta = matrix(as.numeric(T.orig <= C), n, 2)
  T.cen[T.orig[, 1] <= C1, 1] = T.orig[, 1][T.orig[, 1] <= C1]
  T.cen[T.orig[, 1] >= C1, 1] = C1[T.orig[, 1] > C1]
  T.cen[T.orig[, 2] <= C2, 2] = T.orig[, 2][T.orig[, 2] <= C2]
  T.cen[T.orig[, 2] > C2, 2] = C2[T.orig[, 2] > C2]
  T.max = apply(T.cen, 1, max)
  T.min = runif(length(T.max), 0, min(T.cen))
  
  #***************Longitudinal data second*******************************
  ni = (sample(2 : mtimes, n, replace = T))
  t = lapply(1 : n, function(i) seq(T.min[i], T.max[i], length = ni[i]))
  X3 = (rbinom(n, 1, 0.5)); X4 = (runif(n, 0, 2))
  X = cbind(1, unlist(t), rep(X3, ni), rep(X4, ni))

  y = rpois(length(unlist(t)), exp(X %*% alpha.0 + rep(b[, 1], ni) + unlist(t) * rep(b[, 2], ni)))
  
    
  ##organize the longitudinal data and survival data into a table
  obse.no = unlist(lapply(1 : n, function(i) 1 : ni[i]))
  id = rep(1 : n, ni)
  
  Mydata = list(id = id, ni = ni, obse.no = obse.no, X = X, b.long1 = rep(
    b[, 1], ni), b.long2 = rep(b[, 2], ni), rho = rep(1, length(id)), y = y, 
    T.cen = T.cen, W = W, b.surv = b, delta = delta, X_al = cbind(X3, X4), t = t, B1 = cbind(B11, B12, B13), B2 = cbind(B21, B22, B23))
  return(Mydata)
  
}

##To find C
Cen = function(R = 500, n = 100)
{
  cen_r = matrix(0, R, 2)
  for (r in 1 : R) 
  {
    dat = Gdata(n, r ^ 2, 1 : n, CR = c(2, 1.8))
    #print(max(dat$y))
    cen_r[r, ] = apply(dat$delta, 2, mean)
  }
  cr = 1 - apply(cen_r, 2, mean)
  cr
}
# Cen()

Cen_d <- function(n, r, MR)
{
  CR = seq(1, 3, length = 20)
  CRR = cbind(rep(CR, each = length(CR)), rep(CR, length(CR)))
  mr = matrix(0, nrow(CRR), 2)
  for(i in 1 : nrow(CRR))
  {
    dat = Gdata(n, r, CR = CRR[i, ]) 
    mr[i, ] = 1 - apply(dat$delta, 2, mean)
  }
  i.opt = which.min(apply(abs(mr - MR), 1, mean))
  return(Gdata(n, r, CR = c(CRR[i.opt, ])))
}
