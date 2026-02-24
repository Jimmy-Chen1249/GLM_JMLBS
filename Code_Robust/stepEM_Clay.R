stepEM_Clay = function(theta, l, t.1, s.1, t.2, s.2, nMC, rho.region)
{
  
  # Input parameter estimates
  alpha = theta$alpha
  D2 = theta$D2
  sigma2 = theta$sigma2
  gamma.1 = theta$gamma.1
  gamma.2 = theta$gamma.2
  haz.1 = theta$haz.1
  haz.2 = theta$haz.2
  rho = theta$rho
  
  # Longitudinal
  ni = l$ni
  X = l$X
  y = l$y
  XtX.inv = l$XtX.inv
  Xty = l$Xty
  XtZ = l$XtZ
  Z.lon = l$Z.lon
  Z.lont = l$Z.lont
  
  # Time-to-event data
  #For T1
  W.1 = t.1$W.1
  Sur1 = t.1$Sur1
  Sur1.1 = t.1$Sur1.1
  Sur1.1.list = t.1$Sur1.1.list
  q.1 = t.1$q.1
  p.1 = 1
  r.1 = 1
  nev.1 = t.1$nev.1
  nev.uniq.1 = t.1$nev.uniq.1
  tj.1 = t.1$tj.1
  Zdat.sur.1 = s.1$Zdat.sur.1
  Z.sur.1 = s.1$Z.sur.1
  Zt.sur.1 = s.1$Zt.sur.1
  Zdat.1 = s.1$Zdat.1
  Z.1 = s.1$Z.1
  Zt.1 = s.1$Zt.1
  IW.1 = s.1$IW.1
  maxT.1 = t.1$maxT.1
  
  #For T2
  W.2 = t.2$W.2
  Sur2 = t.2$Sur2
  Sur2.1 = t.2$Sur2.1
  Sur2.1.list = t.2$Sur2.1.list
  q.2 = t.2$q.2
  p.2 = 1
  r.2 = 1
  nev.2 = t.2$nev.2
  nev.uniq.2 = t.2$nev.uniq.2
  tj.2 = t.2$tj.2
  Zdat.sur.2 = s.2$Zdat.sur.2
  Z.sur.2 = s.2$Z.sur.2
  Zt.sur.2 = s.2$Zt.sur.2
  Zdat.2 = s.2$Zdat.2
  Z.2 = s.2$Z.2
  Zt.2 = s.2$Zt.2
  IW.2 = s.2$IW.2
  maxT.2 = t.2$maxT.2
  
  
  
  
  # smapling size
  nMC = nMC
  
  t0 = Sys.time()
  
  #*****************************************************
  # Monte Carlo set-up
  #********************** Sampling ***********************
  # Sigma_i (error covariance matrix; diagonal matrix)
  Sigmai = lapply(ni, function(i) {
    diag(x = rep(sigma2, i), ncol = sum(i))
  })
  
  # Inverse-Sigma_i (error precision matrix; diagonal matrix)
  Sigmai.inv = lapply(ni, function(i) {
    diag(x = rep(1 / sigma2, i), ncol = sum(i))
  })
  
  {
    Dinv <- solve(D2)
    Ai <- mapply(FUN = function(zt, s, z) {
      solve((zt %*% s %*% z) + Dinv)
    },
    z = Z.lon, zt = Z.lont, s = Sigmai.inv,
    SIMPLIFY = FALSE)
    
    # MVN mean vector for [y | b]
    Mi <- mapply(function(a, z, s, yi, xi) {
      as.vector(a %*% (z %*% s %*% (yi - xi %*% alpha)))
    },
    a = Ai, z = Z.lont, s = Sigmai.inv, yi = y, xi = X,
    SIMPLIFY = FALSE)
  }
  
  
  # Monte Carlo sample of [b | y]
  Zq = randtoolbox::sobol(nMC, dim = 1, normal = TRUE, scrambling = 1)
  bi.y = mapply(function(m, a) {
    C = chol(a)
    matrix(rep(m, nMC), nrow = nMC, byrow = TRUE) + (Zq %*% C)
  },
  m = Mi, a = Ai,
  SIMPLIFY = FALSE)
  names(bi.y) = names(Ai)
  #**********************Sampling End*****************************************************
  
  # Calculation of W.1^T gamma_11 in Cox model
  Wtgam.1 = mapply(function(w) {
    as.numeric(w %*% gamma.1[1:q.1])
  }, w = W.1, SIMPLIFY = FALSE)
  
  #exp{Z^*{11}(t) * gamma_{12}}
  expZ.1 = mapply(function(z) {
    exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)]))
  }, z = Z.1, SIMPLIFY = FALSE)
  
  # Expanded gamma_y (repeated for each random effect term)
  gam.1.scale = diag(rep(gamma.1[-(1 : (q.1+p.1))], r.1), ncol = sum(r.1))
  #exp{Z_sur(t)^Tb * beta}
  IZ.sur.1 = mapply(function(x, y) {
    t(x %*% y)
  }, x = IW.1, y = Z.sur.1, SIMPLIFY = FALSE)
  # subjects who are censored before first failure time do not contribute anything
  # -> this information is captured through expRhob
  expZ.surb.1 = expWArma(IZ.sur.1, bi.y, gam.1.scale, Sur1.1.list)
  
  # Calculation of W.2^T gamma_21 in Cox model
  Wtgam.2 = mapply(function(w) {
    as.numeric(w %*% gamma.2[1:q.2])
  }, w = W.2, SIMPLIFY = FALSE)
  
  #exp{Z^*{21}(t) * gamma_{22}}
  expZ.2 = mapply(function(z) {
    exp(as.numeric(z %*% gamma.2[(q.2+1):(q.2+p.2)]))
  }, z = Z.2, SIMPLIFY = FALSE)
  
  # Expanded gamma_y (repeated for each random effect term)
  gam.2.scale = diag(rep(gamma.2[-(1 : (q.2+p.2))], r.2), ncol = sum(r.2))
  #exp{Z_sur(t)^Tb * beta}
  IZ.sur.2 = mapply(function(x, y) {
    t(x %*% y)
  }, x = IW.2, y = Z.sur.2, SIMPLIFY = FALSE)
  # subjects who are censored before first failure time do not contribute anything
  # -> this information is captured through expRhob
  expZ.surb.2 = expWArma(IZ.sur.2, bi.y, gam.2.scale, Sur2.1.list)
  
  #**********************Pdf for Survival and longitudianl*************************************
  logf_sur = mapply(function(w.1, z.1, zsur.1, h.1, w.2, z.2, zsur.2, h.2) {
    H.1 = as.vector(t(t(zsur.1) * z.1) %*% haz.1[1:ncol(zsur.1)]) * exp(w.1) # cummulative hazard
    H.2 = as.vector(t(t(zsur.2) * z.2) %*% haz.2[1:ncol(zsur.2)]) * exp(w.2) # cummulative hazard
    T.1 = qnorm(1 - exp(-H.1))
    T.2 = qnorm(1 - exp(-H.2))
    if (h.1$delta == 1 && h.2$delta == 1) { # event
      (log(haz.1[ncol(zsur.1)]) + w.1 + log(z.1[ncol(zsur.1)]) + log(zsur.1[, ncol(zsur.1)])) - H.1 + (log(
        haz.2[ncol(zsur.2)]) + w.2 + log(z.2[ncol(zsur.2)]) + log(zsur.2[, ncol(zsur.2)])) - H.2 + log((rho + 1) * (exp(-H.1) * exp(
          -H.2)) ^ (-(rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^ (-(2 * rho + 1) / rho))
    }
    else if (h.1$delta == 1 && h.2$delta == 0) { # event
      (log(haz.1[ncol(zsur.1)]) + w.1 + log(z.1[ncol(zsur.1)]) + log(zsur.1[, ncol(zsur.1)])) - H.1 + log(exp(-H.1) ^ (-(
        rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-(rho + 1) / rho))
    }
    else if (h.1$delta == 0 && h.2$delta == 1) { # event
      (log(haz.2[ncol(zsur.2)]) + w.2 + log(z.2[ncol(zsur.2)]) + log(zsur.2[, ncol(zsur.2)])) - H.2 + log(exp(-H.2) ^ (-(
        rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-(rho + 1) / rho))
    }
    else if (h.1$delta == 0 && h.2$delta == 0) { # event
      log((exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-1/rho))
    }
  },
  w.1 = Wtgam.1, z.1 = expZ.1, zsur.1 = expZ.surb.1, h.1 = Sur1.1.list,
  w.2 = Wtgam.2, z.2 = expZ.2, zsur.2 = expZ.surb.2, h.2 = Sur2.1.list,
  SIMPLIFY = FALSE)
  
  f_sur = lapply(logf_sur, exp) # f(T_1, T_2, delta_1, delta_2 | b)
  f_pdf = f_sur
  
  
  # Expectation denominator
  den = lapply(f_pdf, mean)
  denkk = length(which(unlist(den) == 0))
  for(denk in which(unlist(den) == 0))
  {
    f_pdf[[denk]] = rep(1, length(f_pdf[[1]]))
    den[[denk]] = 1
  }
  
  # f(T_1, T_2, delta_1, delta_2 | b) / den
  pb.yt = mapply(function(f, d) {
    if(d == 1 | is.nan(d) | is.na(d) | is.infinite(d)) {rep(1, length(f))} else {f / d}
  },
  f = f_pdf, d = den, SIMPLIFY = FALSE)
  
  
  # logf_sur: list，每个元素长度 = nMC（对应该受试者的 MC 样本）
  # 目标 pb.yt: list，每个元素长度 = nMC，且均值为 1（等价于除以 mean）
  
  # pb.yt <- lapply(logf_sur, function(lw) {
  #   lw <- as.numeric(lw)
  #   
  #   # 如果全是 NA/Inf，退化成均匀权重（至少不崩）
  #   ok <- is.finite(lw)
  #   if (!any(ok)) return(rep(1, length(lw)))
  #   
  #   # 用 log-sum-exp 稳定化：w = exp(lw - max)
  #   m <- max(lw[ok])
  #   w <- rep(0, length(lw))
  #   w[ok] <- exp(lw[ok] - m)
  #   
  #   # 归一化为 mean=1（等价于除以 mean(f)）
  #   denom <- mean(w[ok])
  #   if (!is.finite(denom) || denom <= 0) return(rep(1, length(lw)))
  #   
  #   w / denom
  # })
  
  t1 = Sys.time()
  
  #********************************************************************************
  # E-step starts here
  #********************************************************************************
  # E[b]
  Eb = mapply(function(b, pb) {
    colMeans(b * pb)
  },
  b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  
  # E[bb^T]
  EbbT = mapply(function(b, pb) {
    crossprod(b, (b * pb)) / nrow(b)
  },
  b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  
  # exp{W1 %*% gamma_11 + Z^*_11(t) * gamma_12 + Z^*_12(t)b * beta1}(T_1)
  expvstargam.1 = mapply(function(w.1, z.1, zsur.1) {
    exp(w.1) * t(t(zsur.1) * z.1)
  },
  w.1 = Wtgam.1, z.1 = expZ.1, zsur.1 = expZ.surb.1, SIMPLIFY = FALSE)
  
  #exp{W2 %*% gamma_21 + Z^*_21(t) * gamma_22 + Z^*_22(t)b * beta2}(T_2)
  expvstargam.2 = mapply(function(w.2, z.2, zsur.2) {
    exp(w.2) * t(t(zsur.2) * z.2)
  },
  w.2 = Wtgam.2, z.2 = expZ.2, zsur.2 = expZ.surb.2, SIMPLIFY = FALSE)
  
  
  beta1 <- gamma.1[-(1:(q.1 + p.1))]   # 共享随机效应系数（标量）
  
  haz.hat.1 <- update_haz_breslow(
    tj = tj.1,
    Sur = Sur1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ = expZ.1,
    Zsur_list = Z.sur.1,
    pb_y = pb.yt,
    b_y = bi.y,
    beta_shared = beta1
  )

  beta2 <- gamma.2[-(1:(q.2 + p.2))]   # 共享随机效应系数（标量）
  haz.hat.2 <- update_haz_breslow(
    tj = tj.2,
    Sur = Sur2,
    Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2,
    expZ = expZ.2,
    Zsur_list = Z.sur.2,
    pb_y = pb.yt,
    b_y = bi.y,
    beta_shared = beta2
  )
  
 
 
  
  gDelta.1 = gammaUpdate_fixDelta(bi.y, Zt.sur.1, expvstargam.1, pb.yt, haz.hat.1,
                         W.1, Zt.1, Sur1.1.list, 1, q.1, p.1, nev.uniq.1, nev.1)$gDelta
  
  gDelta.2 = gammaUpdate_fixDelta(bi.y, Zt.sur.2, expvstargam.2, pb.yt, haz.hat.2,
                         W.2, Zt.2, Sur2.1.list, 1, q.2, p.2, nev.uniq.2, nev.2)$gDelta
  
  
  
  t2 = Sys.time()
  
  #*****************************************************
  # M-step starts here
  #*****************************************************
  
  # D2
  D2.new = Reduce("+", EbbT) / n
  rownames(D2.new) = colnames(D2.new) = rownames(D2)
  
  #-----------------------------------------------------
  
  # alpha
  rr = mapply(function(x1, x2, b) {
    x1 - (x2 %*% b)
  },
  x1 = Xty, x2 = XtZ, b = Eb)
  rr.sum = rowSums(rr)
  
  alpha.new = as.vector(XtX.inv %*% rr.sum)
  names(alpha.new) = names(alpha)
  
  #---------------------------------------------------------------------------  
  #sigmae^2
  SSq = mapply(function(yi, xi, zi, b, b2) {
    residFixed = (yi - xi %*% alpha.new)
    t(residFixed) %*% (residFixed - 2*(zi %*% b)) + sum(diag(crossprod(zi) %*% b2))
  },
  yi = y, xi = X, zi = Z.lon, b = Eb, b2 = EbbT)
  sigma2.new = sum(SSq) / sum(dat$ni)
  
  #-----------------------------------------------------
  
  # gamma
  gamma.1.new = gamma.1 + as.vector(gDelta.1)
  gamma.2.new = gamma.2 + as.vector(gDelta.2)
  
  #---------------------------------------lam and rho----------------------------------------------
  gamma.1 = gamma.1.new; gamma.2 = gamma.2.new
  # W_{ik}^T %*% gamma_{k1}
  Wtgam.1 = mapply(function(w) {
    as.numeric(w %*% gamma.1[1:q.1])
  }, w = W.1, SIMPLIFY = FALSE)

  #exp{Z^*{11}(t) * gamma_{12}}
  expZ.1 = mapply(function(z) {
    exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)]))
  }, z = Z.1, SIMPLIFY = FALSE)

  #exp{Z.sur^*{k1}(t) * bi * beta_k}
  expZurb.1 = mapply(function(b, pb){
    mean(exp(b * gamma.1[-(1 : (q.1+p.1))]) * pb)
  }, b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)

  # W_{ik}^T %*% gamma_{k1} in Cox model
  Wtgam.2 = mapply(function(w) {
    as.numeric(w %*% gamma.2[1:q.2])
  }, w = W.2, SIMPLIFY = FALSE)

  #exp{Z^*{21}(t) * gamma_{22}}
  expZ.2 = mapply(function(z) {
    exp(as.numeric(z %*% gamma.2[(q.2+1):(q.2+p.2)]))
  }, z = Z.2, SIMPLIFY = FALSE)

  #exp{Z.sur^*{k1}(t) * bi * beta_k}
  expZurb.2 = mapply(function(b, pb){
    mean(exp(b * gamma.2[-(1 : (q.2+p.2))]) * pb)
  }, b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)


  gam.1.scale = diag(rep(gamma.1[-(1 : (q.1+p.1))], r.1), ncol = sum(r.1))
  gam.2.scale = diag(rep(gamma.2[-(1 : (q.2+p.2))], r.2), ncol = sum(r.2))
 
  beta1 <- gamma.1[-(1:(q.1 + p.1))]  # 共享随机效应系数（标量）
  
  haz.1.new <- update_haz_breslow(
    tj = tj.1,
    Sur = Sur1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ = expZ.1,
    Zsur_list = Z.sur.1,
    pb_y = pb.yt,
    b_y = bi.y,
    beta_shared = beta1
  )
  
  beta2 <- gamma.2[-(1:(q.2 + p.2))]
  
  haz.2.new <- update_haz_breslow(
    tj = tj.2,
    Sur = Sur2,
    Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2,
    expZ = expZ.2,
    Zsur_list = Z.sur.2,
    pb_y = pb.yt,
    b_y = bi.y,
    beta_shared = beta2
  )
  
  #************************** rho ************************************
  #For T1
  Z.1.new = mapply(function(zi) {exp(zi * gamma.1[(q.1+1):(q.1+p.1)])
  }, zi = Z.1)
  # haz.1.new = lambdaUpdate(bi.y, IW.1, Z.sur.1, pb.yt, W.1, Z.1.new, gam.1.scale,
  #                          gamma.1, q.1, nev.1, Sur1.1.list)
  # haz.1.new = as.vector(haz.1.new)
  
  #For T2
  Z.2.new = mapply(function(zi) {exp(zi * gamma.2[(q.2+1):(q.2+p.2)])
  }, zi = Z.2)
  # haz.2.new = lambdaUpdate(bi.y, IW.2, Z.sur.2, pb.yt, W.2, Z.2.new, gam.2.scale,
  #                          gamma.2, q.2, nev.2, Sur2.1.list)
  # haz.2.new = as.vector(haz.2.new)
  
  haz.1 = haz.1.new; haz.2 = haz.2.new;
  
  d_11 = which((Sur1.1$delta + Sur2.1$delta) == 2)
  d_10 = which((Sur1.1$delta - Sur2.1$delta) == 1)
  d_01 = which((Sur1.1$delta - Sur2.1$delta) == -1)
  d_00 = which((Sur1.1$delta + Sur2.1$delta) == 0)
  
  
  log_rho = function(rho)
  {
    k = mapply(function(w.1, z.1, zsur.1, h.1, w.2, z.2, zsur.2, h.2) {
      H.1 = as.vector((zsur.1 * z.1) %*% haz.1[1:length(z.1)]) * exp(w.1)  # cummulative hazard
      H.2 = as.vector((zsur.2 * z.2) %*% haz.2[1:length(z.2)]) * exp(w.2) # cummulative hazard
      if (h.1$delta == 1 && h.2$delta == 1) { # event
        (log(haz.1[length(z.1)]) + w.1 + log(z.1[length(z.1)]) + log(zsur.1)) - H.1 + (log(
          haz.2[length(z.2)]) + w.2 + log(z.2[length(z.2)]) + log(zsur.2)) - H.2 + log((rho + 1) * (exp(-H.1) * exp(
            -H.2)) ^ (-(rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^ (-(2 * rho + 1) / rho))
      }
      else if (h.1$delta == 1 && h.2$delta == 0) { # event
        (log(haz.1[length(z.1)]) + w.1 + log(z.1[length(z.1)]) + log(zsur.1)) - H.1 + log(exp(-H.1) ^ (-(
          rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-(rho + 1) / rho))
      }
      else if (h.1$delta == 0 && h.2$delta == 1) { # event
        (log(haz.2[length(z.2)]) + w.2 + log(z.2[length(z.2)]) + log(zsur.2)) - H.2 + log(exp(-H.2) ^ (-(
          rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-(rho + 1) / rho))
      }
      else if (h.1$delta == 0 && h.2$delta == 0) { # event
        log((exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-1/rho))
      }
    },
    w.1 = Wtgam.1, z.1 = expZ.1, zsur.1 = expZurb.1, h.1 = Sur1.1.list,
    w.2 = Wtgam.2, z.2 = expZ.2, zsur.2 = expZurb.2, h.2 = Sur2.1.list,
    SIMPLIFY = FALSE)
    return(-sum(unlist(k)))
  }
  rho.region = c(1e-4, 30)
  x.inv <- try(optimize(log_rho, rho.region)$minimum, silent=TRUE)
  if ('try-error' %in% class(x.inv)) {rho.new <- 10} else {rho.new <- optimize(log_rho, rho.region)$minimum}
  
  
 
  theta.new = list("D2" = D2.new, "alpha" = alpha.new, "sigma2" = sigma2.new,
                   "haz.1" = haz.1.new, "gamma.1" = gamma.1.new,
                   "haz.2" = haz.2.new, "gamma.2" = gamma.2.new, rho = rho.new)
  out = list(theta = theta.new)
  return(out)
}
