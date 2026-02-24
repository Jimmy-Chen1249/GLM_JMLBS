#***************This program is to get finial estimate***************************
mjoint_Gau = function(Mixeffect, data, dat, Sur1, Sur2, formSurv1, formSurv2, nMC, iter_max)
{
  #*********************Get initial value*****************************
  #****************Longitudinal first*******************
  # lfit = lmer(formula = Mixeffect, data = data)
  lfit = nlme::lme(fixed = y ~ time + X1 + X2, random = ~1|id,
                   data = data, method = "ML",
                   control = nlme::lmeControl(opt = "optim"))
  
  # Longitudinal outcomes
  y = by(data$y, data$id, as.vector)
  X = data.frame(id = data$id, X.Intercept = data$intercept, time = data$time,
                 X1 = data$X1, X2 = data$X2)
  X = by(X, X$id, function(u){as.matrix(u[, -1])})
  
  Z.lon = data.frame(id = data$id, X.Intercept = data$intercept)
  Z.lon = by(Z.lon, Z.lon$id, function(u) {as.matrix(u[, -1])})
  Z.lont = lapply(Z.lon, t)
  
  XtX = lapply(X, crossprod)
  XtX.inv = solve(Reduce("+", XtX))
  Xty = mapply(function(x, y) {
    crossprod(x, y)
  },
  x = X, y = y,
  SIMPLIFY = FALSE)
  XtZ = mapply(function(x, z) {
    crossprod(x, z)
  },
  x = X, z = Z.lon,
  SIMPLIFY = FALSE)
  l = list(ni = dat$ni, y = y, X = X, Z.lon = Z.lon, Z.lont = Z.lont, 
           XtX.inv = XtX.inv, Xty = Xty, XtZ = XtZ)
  
  #*******************Survival Second*********************************
  #Estimate the paramters in Cox model
  jitter_surtime <- function(surtime, eps = 1e-4) {
    dup_idx <- duplicated(surtime) | duplicated(surtime, fromLast = TRUE)
    if (any(dup_idx)) {
      surtime[dup_idx] <- surtime[dup_idx] + rnorm(sum(dup_idx), mean = 0, sd = eps)
    }
    surtime <- pmax(surtime, .Machine$double.eps)
    return(surtime)
  }
  #***************For S1
  Sur1$surtime = jitter_surtime(Sur1$surtime)
  ##Prapare for further estimating
  sfit.1 = survival::coxph(formSurv1, data = Sur1, x = TRUE)
  q.1 = ncol(sfit.1$x)
  sfit.start.1 = survival::survfit(sfit.1)
  tj.1 = sfit.start.1$time[sfit.start.1$n.event > 0]
  nev.1 = sfit.start.1$n.event[sfit.start.1$n.event > 0]
  nev.uniq.1 = length(tj.1)
  Sur1.1 = data.frame(Sur1$id, sfit.1$x, sfit.1$y[, 1], sfit.1$y[, 2])
  
  
  colnames(Sur1.1)[c(1, (q.1 + 2):(q.1 + 3))] = c("id", "T", "delta")
  
  Sur1.1$tj.ind = sapply(1:n, function(i) {
    sum(tj.1 <= Sur1.1$T[i])
  })
  Sur1.1.list = by(Sur1.1, Sur1.1$id, list)
  W.1 = by(Sur1.1, Sur1.1$id, function(u) {
    unlist(u[, 2:(q.1 + 1)])
  })
  
  maxT.1 = which.max(Sur1$surtime * Sur1$cens)
  t.1 = list(W.1 = W.1, Sur1 = Sur1, Sur1.1 = Sur1.1, Sur1.1.list = Sur1.1.list,
             q.1 = q.1, nev.1 = nev.1, nev.uniq.1 = nev.uniq.1, tj.1 = tj.1, maxT.1 = maxT.1)
  
  Zdat.sur.1 = data.frame(
    "id" = rep(unique(Sur1.1$id), pmax(Sur1.1$tj.ind, 1)),
    "time" = unlist(sapply(1:n, function(i) {
      tj.1[1:max(Sur1.1$tj.ind[i], 1)]
    },
    simplify = FALSE)))
  names(Zdat.sur.1)[1] = "id"
  names(Zdat.sur.1)[2] = "time"
  
  Z.sur.1 = data.frame("id" = Zdat.sur.1[, "id"], X.Intercept = 1)
  Z.sur.1 = by(Z.sur.1, Z.sur.1$id, function(u) {
    as.matrix(u[, -1])})
  Zt.sur.1 = lapply(Z.sur.1, t)
  
  
  Zdat.1 = data.frame("id" = Zdat.sur.1[, "id"], Z.1 = unlist(lapply(1 : n, function(i) {
    temp.1 = (Zdat.sur.1[which(Zdat.sur.1[, "id"] == i), ])[, 2];
    temp.1})))
  
  Z.1 = by(Zdat.1, Zdat.1$id, function(u) {
    B1 = Sur1$B1[u[1,1]]; B2 = Sur1$B2[u[1,1]]; B3 = Sur1$B3[u[1,1]]
    as.matrix(B1 * (u[, -1] <= B3) + B2 * (u[, -1] > B3))})
  Zt.1 = lapply(Z.1, t)
  IW.1 = by(Sur1.1, Sur1.1$id, function(u) {
    do.call("cbind", lapply(1:1, function(i) diag(max(u$tj.ind, 1))))
  })
  
  s.1 = list(Zdat.sur.1 = Zdat.sur.1, Z.sur.1 = Z.sur.1, Zt.sur.1  =Zt.sur.1,
             Zdat.1 = Zdat.1, Z.1 = Z.1, Zt.1 = Zt.1, IW.1 = IW.1)
  
  su1.para = initsSurv_new(data, Sur1, lfit, Sur1.1, "id", "time",  1)
  
  rm(W.1, Sur1, Sur1.1, Sur1.1.list, q.1, nev.1, nev.uniq.1, tj.1,
     Zdat.sur.1, Z.sur.1, Zt.sur.1, Zdat.1, Z.1, IW.1)
  
  #***************For S2
  Sur2$surtime = jitter_surtime(Sur2$surtime)
  ##Prapare for further estimating
  sfit.2 = survival::coxph(formSurv2, data = Sur2, x = TRUE)
  q.2 = ncol(sfit.2$x)
  sfit.start.2 = survival::survfit(sfit.2)
  tj.2 = sfit.start.2$time[sfit.start.2$n.event > 0]
  nev.2 = sfit.start.2$n.event[sfit.start.2$n.event > 0]
  nev.uniq.2 = length(tj.2)
  Sur2.1 = data.frame(Sur2$id, sfit.2$x, sfit.2$y[, 1], sfit.2$y[, 2])
  
  
  colnames(Sur2.1)[c(1, (q.2 + 2):(q.2 + 3))] = c("id", "T", "delta")
  
  Sur2.1$tj.ind = sapply(1:n, function(i) {
    sum(tj.2 <= Sur2.1$T[i])
  })
  Sur2.1.list = by(Sur2.1, Sur2.1$id, list)
  W.2 = by(Sur2.1, Sur2.1$id, function(u) {
    unlist(u[, 2:(q.2 + 1)])
  })
  
  maxT.2 = which.max(Sur2$surtime * Sur2$cens)
  t.2 = list(W.2 = W.2, Sur2 = Sur2, Sur2.1 = Sur2.1, Sur2.1.list = Sur2.1.list,
             q.2 = q.2, nev.2 = nev.2, nev.uniq.2 = nev.uniq.2, tj.2 = tj.2, maxT.2 = maxT.2)
  
  Zdat.sur.2 = data.frame(
    "id" = rep(unique(Sur2.1$id), pmax(Sur2.1$tj.ind, 1)),
    "time" = unlist(sapply(1:n, function(i) {
      tj.2[1:max(Sur2.1$tj.ind[i], 1)]
    },
    simplify = FALSE)))
  names(Zdat.sur.2)[1] = "id"
  names(Zdat.sur.2)[2] = "time"
  
  Z.sur.2 = data.frame("id" = Zdat.sur.2[, "id"], X.Intercept = 1)
  Z.sur.2 = by(Z.sur.2, Z.sur.2$id, function(u) {
    as.matrix(u[, -1])})
  Zt.sur.2 = lapply(Z.sur.2, t)
  
  Zdat.2 = data.frame("id" = Zdat.sur.2[, "id"], Z.2 = unlist(lapply(1 : n, function(i) {
    temp.2 = (Zdat.sur.2[which(Zdat.sur.2[, "id"] == i), ])[, 2];
    temp.2})))
  Z.2 = by(Zdat.2, Zdat.2$id, function(u) {
    B1 = Sur2$B1[u[1,1]]; B2 = Sur2$B2[u[1,1]]; B3 = Sur2$B3[u[1,1]]
    as.matrix(B1 * (u[, -1] <= B3) + B2 * (u[, -1] > B3))})
  Zt.2 = lapply(Z.2, t)
  IW.2 = by(Sur2.1, Sur2.1$id, function(u) {
    do.call("cbind", lapply(1:1, function(i) diag(max(u$tj.ind, 1))))
  })
  
  s.2 = list(Zdat.sur.2 = Zdat.sur.2, Z.sur.2 = Z.sur.2, Zt.sur.2  =Zt.sur.2,
             Zdat.2 = Zdat.2, Z.2 = Z.2, Zt.2 = Zt.2, IW.2 = IW.2)
  su2.para = initsSurv_new(data, Sur2, lfit, Sur2.1, "id", "time",  1)
  
  rm(W.2, Sur2, Sur2.1, Sur2.1.list, q.2, nev.2, nev.uniq.2, tj.2,
     Zdat.sur.2, Z.sur.2, Zt.sur.2, Zdat.2, Z.2, IW.2)
  
  #***********************************************************
  alpha0 = fixef(lfit)
  sigma20 = lfit$sigma^2
  D20 = getVarCov(lfit)
  gamma.10 = su1.para$gamma
  haz.10 = su1.para$haz
  gamma.20 = su2.para$gamma
  haz.20 = su2.para$haz
  rho = 0.2; rho.region = c(-0.5, 0.5);
  theta = list(alpha = alpha0, sigma2 = sigma20, D2 = D20, gamma.1 = gamma.10,
               haz.1 = haz.10, gamma.2 = gamma.20, haz.2 = haz.20, rho = rho)
  
  #Start estimation
  #*****************************************************************************
  #######iteration#########
  alpha.old = theta$alpha
  sigma2.old = theta$sigma2
  D2.old = theta$D2
  gamma.1.old = theta$gamma.1
  gamma.2.old = theta$gamma.2
  haz.1.old = theta$haz.1
  haz.2.old = theta$haz.2
  rho.old = theta$rho
  
  tolerance = 1
  iter = 0
  tol = 5e-3
  #iter_max = 500
  
  print(c(theta$alpha, theta$sigma2, theta$D2, theta$gamma.1, theta$gamma.2, theta$rho))
  while(tolerance>=tol & iter <= iter_max)
  {
    out = stepEM_Gau(theta = theta, l = l, t.1 = t.1, s.1 = s.1,
                     t.2 = t.2, s.2 = s.2, nMC = nMC, rho.region = rho.region)
    theta.new = out$theta
    alpha.new = theta.new$alpha
    sigma2.new = theta.new$sigma2
    D2.new = theta.new$D2
    gamma.1.new = theta.new$gamma.1
    gamma.2.new = theta.new$gamma.2
    haz.1.new = theta.new$haz.1
    haz.2.new = theta.new$haz.2
    rho.new = theta.new$rho
    
    abssub = c(abs(alpha.new - alpha.old), abs(sigma2.new - sigma2.old),
               abs(D2.new - D2.old), abs(gamma.1.new - gamma.1.old),
               abs(gamma.2.new - gamma.2.old), abs(rho.new - rho.old))
    rbssub = abssub / abs(c(alpha.old, sigma2.old, D2.old, gamma.1.old, gamma.2.old, rho.old))
    # print(c(max(abssub), max(rbssub)))
    tolerance = min(max(abssub), max(rbssub))#max(max(abssub), max(c(abssub, ee1)))
    # print(abssub)
    # print(c(alpha.new, sigma2.new, D2.new, gamma.1.new, gamma.2.new, rho.new))
    # print(c(gamma.1.new, gamma.2.new))
    iter = iter+1
    # print(tolerance)
    
    theta = theta.new
    alpha.old = theta$alpha
    sigma2.old = theta$sigma2
    D2.old = theta$D2
    gamma.1.old = theta$gamma.1
    gamma.2.old = theta$gamma.2
    haz.1.old = theta$haz.1
    haz.2.old = theta$haz.2
    rho.old = theta$rho
    if(tolerance < tol | iter > iter_max){message("EM algorithm has converged!\n")
      theta = theta.new
      break} else {
        nMC = min(nMC + floor(nMC / 3), 2000)
        # print(c(nMC, iter))
      }
  }
  return(Re = list(theta = theta, iter = iter))
}