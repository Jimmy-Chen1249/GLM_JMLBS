#***************This program is to get finial estimate***************************
mjoint = function(Mixeffect, data, dat, Sur1, Sur2, formSurv1, formSurv2, nMC, iter_max)
{
  start_time <- Sys.time()
  time_limit <- 20000
  #*********************Get initial value*****************************
  #****************Longitudinal first*******************
  lfit = glmer(formula = Mixeffect, data = data, family = poisson)
  
  # Longitudinal outcomes
  y = by(data$y, data$id, as.vector)
  X = data.frame(id = data$id, X.Intercept = data$intercept, time = data$time,
                 X1 = data$X1, X2 = data$X2)
  X = by(X, X$id, function(u){as.matrix(u[, -1])})
  
  Z.lon = data.frame(id = data$id, X.Intercept = data$intercept, time = data$time)
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
  
  #centering
  xcenter.1 <- NULL
  xcenter.1 <- apply(Sur1.1[2:(q.1 + 1)], 2, mean)
  Sur1.1[2:(q.1 + 1)] <- scale(Sur1.1[2:(q.1 + 1)], center = xcenter.1, scale = FALSE)
  
  Sur1.1$tj.ind = sapply(1:n, function(i) {
    sum(tj.1 <= Sur1.1$T[i])
  })
  Sur1.1.list = by(Sur1.1, Sur1.1$id, list)
  W.1 = by(Sur1.1, Sur1.1$id, function(u) {
    unlist(u[, 2:(q.1 + 1)])
  })
  maxT.1 = which.max(dat$T.cen[,1] * dat$delta[,1])
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
  #Zt.sur.1 = lapply(Z.sur.1, t)
  
  
  Zdat.1 = data.frame("id" = Zdat.sur.1[, "id"], Z.1 = unlist(lapply(1 : n, function(i) {
    temp.1 = (Zdat.sur.1[which(Zdat.sur.1[, "id"] == i), ])[, 2];
    temp.1})))
  Z.1 = by(Zdat.1, Zdat.1$id, function(u) {
    B1 = Sur1$B1[u[1,1]]; B2 = Sur1$B2[u[1,1]]; B3 = Sur1$B3[u[1,1]]
    as.matrix(B1 * (u[, -1] <= B3) + B2 * (u[, -1] > B3))})
  Zt.1 = lapply(Z.1, t)
  
  Z.t1 = by(Zdat.1, Zdat.1$id, function(u) {
    as.matrix(u[, -1])})
  #add the function of t into Z.sur.1
  Z.sur.1 = mapply(function(x, id) {
    y = Zdat.sur.1$time[which(Zdat.sur.1$id == id)]
    cbind(x, y)
  },
  x = Z.sur.1, id = 1 : n, SIMPLIFY = FALSE)
  Zt.sur.1 = lapply(Z.sur.1, t)
  
  IW.1 = by(Sur1.1, Sur1.1$id, function(u) {
    do.call("cbind", lapply(1:1, function(i) diag(max(u$tj.ind, 1))))
  })
  
  s.1 = list(Zdat.sur.1 = Zdat.sur.1, Z.sur.1 = Z.sur.1, Zt.sur.1  =Zt.sur.1,
             Zdat.1 = Zdat.1, Z.1 = Z.1, Zt.1 = Zt.1, IW.1 = IW.1)
  su1.para = initsSurv(data, Sur1, lfit, Sur1.1, "id", "time",  1)
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
  
  #centering
  xcenter.2 <- NULL
  xcenter.2 <- apply(Sur2.1[2:(q.2 + 1)], 2, mean)
  Sur2.1[2:(q.2 + 1)] <- scale(Sur2.1[2:(q.2 + 1)], center = xcenter.2, scale = FALSE)
  
  Sur2.1$tj.ind = sapply(1:n, function(i) {
    sum(tj.2 <= Sur2.1$T[i])
  })
  Sur2.1.list = by(Sur2.1, Sur2.1$id, list)
  W.2 = by(Sur2.1, Sur2.1$id, function(u) {
    unlist(u[, 2:(q.2 + 1)])
  })
  maxT.2 = which.max(dat$T.cen[,2] * dat$delta[,2])
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
  # Zt.sur.2 = lapply(Z.sur.2, t)
  
  
  Zdat.2 = data.frame("id" = Zdat.sur.2[, "id"], Z.2 = unlist(lapply(1 : n, function(i) {
    temp.2 = (Zdat.sur.2[which(Zdat.sur.2[, "id"] == i), ])[, 2];
    temp.2})))
  Z.2 = by(Zdat.2, Zdat.2$id, function(u) {
    B1 = Sur2$B1[u[1,1]]; B2 = Sur2$B2[u[1,1]]; B3 = Sur2$B3[u[1,1]]
    as.matrix(B1 * (u[, -1] <= B3) + B2 * (u[, -1] > B3))})
  Zt.2 = lapply(Z.2, t)
  
  Z.t2 = by(Zdat.2, Zdat.2$id, function(u) {
    as.matrix(u[, -1])})
  #add the function of t into Z.sur
  Z.sur.2 = mapply(function(x, id) {
    y = Zdat.sur.2$time[which(Zdat.sur.2$id == id)]
    cbind(x, y)
  },
  x = Z.sur.2, id = 1 : n, SIMPLIFY = FALSE)
  Zt.sur.2 = lapply(Z.sur.2, t)
  
  IW.2 = by(Sur2.1, Sur2.1$id, function(u) {
    do.call("cbind", lapply(1:1, function(i) diag(max(u$tj.ind, 1))))
  })
  
  s.2 = list(Zdat.sur.2 = Zdat.sur.2, Z.sur.2 = Z.sur.2, Zt.sur.2  =Zt.sur.2,
             Zdat.2 = Zdat.2, Z.2 = Z.2, Zt.2 = Zt.2, IW.2 = IW.2)
  su2.para = initsSurv(data, Sur2, lfit, Sur2.1, "id", "time",  1)
  rm(W.2, Sur2, Sur2.1, Sur2.1.list, q.2, nev.2, nev.uniq.2, tj.2,
     Zdat.sur.2, Z.sur.2, Zt.sur.2, Zdat.2, Z.2, IW.2)
  
  #***********************************************************
  Sigma.b0 = matrix(unlist(summary(lfit)$varcor), nrow = sqrt(length(unlist(summary(lfit)$varcor))))
  alpha0 <- summary(lfit)$coefficients[, 1]
  gamma.10 = su1.para$gamma #c(gamma.011, gamma.012, beta.01)# su1.para$gamma
  haz.10 = su1.para$haz
  gamma.20 = su2.para$gamma #c(gamma.021, gamma.022, beta.02)#su2.para$gamma
  haz.20 = su2.para$haz
  if(rho.0 < 0) {rho = -0.5; rho.region <- c(-1, 0)} else if(rho.0 == 0) {
    rho = 0; rho.region <- c(-0.5, 0.5)} else if(rho.0 > 0) {rho = 0.5; rho.region <- c(0, 1)}
  theta = list(alpha = alpha0, Sigma.b = Sigma.b0, gamma.1 = gamma.10, rho = rho,
               haz.1 = haz.10, gamma.2 = gamma.20, haz.2 = haz.20)
  
  #Start estimation
  #*****************************************************************************
  #######iteration#########
  alpha.old = theta$alpha
  Sigma.b.old = theta$Sigma.b
  gamma.1.old = theta$gamma.1
  gamma.2.old = theta$gamma.2
  haz.1.old = theta$haz.1
  haz.2.old = theta$haz.2
  rho.old = theta$rho
  
  tolerance = 1
  iter = 0
  tol = 5e-3; nMCmax = 500; iter1 = iter2 = 0
  while(tolerance>=tol & iter <= iter_max)
  {
    if(iter  > floor(iter_max * 0.3) && iter1 != 1)  {nMCmax = nMCmax* 2; iter1 =1}
    if(iter  > floor(iter_max * 0.5) && iter2 != 1) {nMCmax = nMCmax* 2; iter2 =1}
    
    # 想给的超时阈值（秒）
    timeout_sec <- 60
    t0 <- Sys.time()
    res <- tryCatch({
      # 给当前这次调用设置“软”超时；只对这一次调用生效（transient=TRUE）
      setTimeLimit(elapsed = timeout_sec, transient = TRUE)
      stepEM(theta = theta, l = l, t.1 = t.1, s.1 = s.1,
             t.2 = t.2, s.2 = s.2, nMC = nMC, rho.region = rho.region)
    }, error = function(e) {
      # 如果是“到达时间上限”的报错，返回 NULL 触发降级重跑；其他错误原样抛出
      msg <- conditionMessage(e)
      if (grepl("reached elapsed time limit", msg)) return(NULL)
      stop(e)
    })
    # 及时恢复时间限制（取消）
    setTimeLimit()
    
    # 如果超时（res 为 NULL），或虽完成但耗时仍 > timeout_sec，则用 nMC = 100 重跑
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(res) || elapsed > timeout_sec) {
      nMC = 100
      out <- stepEM(theta = theta, l = l, t.1 = t.1, s.1 = s.1,
                    t.2 = t.2, s.2 = s.2, nMC = nMC, rho.region = rho.region)
    } else {
      out <- res
    }
    
    # out = stepEM(theta = theta, l = l, t.1 = t.1, s.1 = s.1,
    #              t.2 = t.2, s.2 = s.2, nMC = nMC, rho.region = rho.region)
    theta.new = out$theta
    alpha.new = theta.new$alpha
    Sigma.b.new = theta.new$Sigma.b
    gamma.1.new = theta.new$gamma.1
    gamma.2.new = theta.new$gamma.2
    haz.1.new = theta.new$haz.1
    haz.2.new = theta.new$haz.2
    rho.new = theta.new$rho
    abssub = c(abs(alpha.new - alpha.old), abs(as.numeric(Sigma.b.new - Sigma.b.old)), 
               abs(gamma.1.new - gamma.1.old), abs(gamma.2.new - gamma.2.old),
               abs(rho.new - rho.old))
    rbssub = abssub / c(abs(alpha.new), abs(as.numeric(Sigma.b.new)), 
                        abs(gamma.1.new), abs(gamma.2.new ), abs(rho.new ))
    
    tolerance = min(max(abssub), max(rbssub))#max(max(abssub), max(c(abssub, ee1)))
    print(abssub)
    print(c(alpha.new, c(Sigma.b.new), gamma.1.new, gamma.2.new, rho.new))
    # print(rho.new)
    # print(list(alpha.new = theta.new$alpha,
    #            Sigma.b.new = theta.new$Sigma.b,
    #            gamma.1.new = theta.new$gamma.1,
    #            gamma.2.new = theta.new$gamma.2,
    #            rho.new = theta.new$rho))
    #print(c(alpha.new, D2.new, gamma.1.new, gamma.2.new, rho.new))
    iter = iter+1
    print(c(tolerance, nMC))
    
    theta = theta.new
    alpha.old = theta$alpha
    Sigma.b.old = theta$Sigma.b
    gamma.1.old = theta$gamma.1
    gamma.2.old = theta$gamma.2
    haz.1.old = theta$haz.1
    rho.old = theta$rho
    if(tolerance < tol | iter > iter_max){message("EM algorithm has converged!\n")
      theta = theta.new
      break} else {
        nMC = min(nMC + floor(nMC / 3), nMCmax)
        print(c(iter, nMCmax))
      }
    if (as.numeric(Sys.time() - start_time, units = "secs") > time_limit) {
      theta = c(NA, 1)
      break
    }
  }
  return(Re = list(theta = theta, iter = iter))
}