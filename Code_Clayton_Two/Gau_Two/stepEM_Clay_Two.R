stepEM_Clay_Two_old = function(theta, l, t.1, s.1, t.2, s.2, nMC, rho.region)
{
  
  # Input parameter estimates
  alpha = theta$alpha
  Sigma.b = theta$Sigma.b
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
    Dinv <- solve(Sigma.b)
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
  Zq = randtoolbox::sobol(nMC, dim = ncol(Sigma.b), normal = TRUE, scrambling = 1)
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
  gam.1.scale = diag(rep(gamma.1[-(1 : (q.1+p.1))], r.1), ncol(Sigma.b), ncol(Sigma.b))
  
  #exp{Z_sur(t)^Tb * beta}
  # IZ.sur.1 = mapply(function(x, y, z) {
  #   rbind(t(x %*% y), t(z))
  # }, x = IW.1, y = Z.sur.1, z = Z.1, SIMPLIFY = FALSE)
  
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
  gam.2.scale = diag(rep(gamma.2[-(1 : (q.2+p.2))], r.2), 2, 2)
  #exp{Z_sur(t)^Tb * beta}
  # IZ.sur.2 = mapply(function(x, y, z) {
  #   rbind(t(x %*% y), t(z))
  # }, x = IW.2, y = Z.sur.2, z = Z.2, SIMPLIFY = FALSE)
  
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
  
  
  time_sur_list_1 <- split(Zdat.sur.1$time, Zdat.sur.1$id)
  
  beta1 <- gamma.1[-(1:(q.1 + p.1))]  # 标量
  haz.hat.1 <- update_haz_breslow_one_t(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ = expZ.1,
    time_sur_list = time_sur_list_1,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta1
  )
  
  time_sur_list_2 <- split(Zdat.sur.2$time, Zdat.sur.2$id)
  
  beta2 <- gamma.2[-(1:(q.2 + p.2))]  # 标量
  haz.hat.2 <- update_haz_breslow_one_t(
    tj = tj.2,
    Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2,
    expZ = expZ.2,
    time_sur_list = time_sur_list_2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  
  # haz.hat.1 = unlist(lapply(1 : nrow(rSur1.1), function(i) 1 / sum(exp(unlist(Wtgam.1))[rSur1.1$id[i:nev.uniq.1]] * (
  #   Ebgamma1[, i])[rSur1.1$id[i:nev.uniq.1]] * rep(expZ.1[[maxT.1]][i], (nev.uniq.1- i + 1)))))
  # 
  # haz.hat.2 = unlist(lapply(1 : nrow(rSur2.1), function(i) 1 / sum(exp(unlist(Wtgam.2))[rSur2.1$id[i:nev.uniq.2]] * (
  #   Ebgamma2[, i])[rSur2.1$id[i:nev.uniq.2]] * rep(expZ.2[[maxT.2]][i], (nev.uniq.2- i + 1)))))
  
  
  gDelta.1 = gammaUpdate(bi.y, Zt.sur.1, expvstargam.1, pb.yt, haz.hat.1,
                         W.1, Zt.1, Sur1.1.list, 1, q.1, p.1, nev.uniq.1, nev.1)$gDelta
  gDelta.2 = gammaUpdate(bi.y, Zt.sur.2, expvstargam.2, pb.yt, haz.hat.2,
                         W.2, Zt.2, Sur2.1.list, 1, q.2, p.2, nev.uniq.2, nev.2)$gDelta
  
  t2 = Sys.time()
  
  #*****************************************************
  # M-step starts here
  #*****************************************************
  
  # D2
  Sigma.b.new = Reduce("+", EbbT) / n
  rownames(Sigma.b.new) = colnames(Sigma.b.new) = rownames(Sigma.b)
  
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
  #For haz1&haz2 in T1 T2
  gamma.1 = gamma.1.new; gamma.2 = gamma.2.new
  
  # Calculation of W.1^T gamma_11 in Cox model
  Wtgam.1 = mapply(function(w) {
    as.numeric(w %*% gamma.1[1:q.1])
  }, w = W.1, SIMPLIFY = FALSE)
  
  #exp{Z^*{11}(t) * gamma_{12}}
  expZ.1 = mapply(function(z) {
    exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)]))
  }, z = Z.1, SIMPLIFY = FALSE)
  
  #exp{Z.sur^*{k1}(t) * bi * beta_k}
  expZurb.1 = mapply(function(b, pb, zi){
    apply(exp(b %*% t(zi) * gamma.1[-(1 : (q.1+p.1))]) * pb, 2, mean)
  }, b = bi.y, pb = pb.yt, zi = Z.sur.1,  SIMPLIFY = FALSE) 
  
  # Expanded gamma_y (repeated for each random effect term)
  gam.1.scale = diag(rep(gamma.1[-(1 : (q.1+p.1))], r.1), ncol(Sigma.b), ncol(Sigma.b))
  
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
  
  #exp{Z.sur^*{k1}(t) * bi * beta_k}
  expZurb.2 = mapply(function(b, pb, zi){
    apply(exp(b %*% t(zi) * gamma.2[-(1 : (q.2+p.2))]) * pb, 2, mean)
  }, b = bi.y, pb = pb.yt, zi = Z.sur.2, SIMPLIFY = FALSE) 
  
  # Expanded gamma_y (repeated for each random effect term)
  gam.2.scale = diag(rep(gamma.2[-(1 : (q.2+p.2))], r.2), 2, 2)
  
  IZ.sur.2 = mapply(function(x, y) {
    t(x %*% y)
  }, x = IW.2, y = Z.sur.2, SIMPLIFY = FALSE)
  # subjects who are censored before first failure time do not contribute anything
  # -> this information is captured through expRhob
  expZ.surb.2 = expWArma(IZ.sur.2, bi.y, gam.2.scale, Sur2.1.list)
  
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
  
  
  ## --- T1 baseline hazard update (Z_sur = (1, t)) ---
  time_sur_list_1 <- split(Zdat.sur.1$time, Zdat.sur.1$id)
  
  beta1 <- gamma.1[-(1:(q.1 + p.1))]  # 标量
  haz.1.new <- update_haz_breslow_one_t(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ = expZ.1,
    time_sur_list = time_sur_list_1,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta1
  )
  
  ## --- T2 baseline hazard update (Z_sur = (1, t)) ---
  time_sur_list_2 <- split(Zdat.sur.2$time, Zdat.sur.2$id)
  
  beta2 <- gamma.2[-(1:(q.2 + p.2))]  # 标量
  haz.2.new <- update_haz_breslow_one_t(
    tj = tj.2,
    Sur1.1 = Sur2.1,      # 注意：这里传 Sur2.1
    Wtgam = Wtgam.2,
    expZ = expZ.2,
    time_sur_list = time_sur_list_2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  
  #************************** rho ************************************
  #For T1
  Z.1.new = mapply(function(zi) {exp(zi * gamma.1[(q.1+1):(q.1+p.1)])
  }, zi = Z.1)

  
  #For T2
  Z.2.new = mapply(function(zi) {exp(zi * gamma.2[(q.2+1):(q.2+p.2)])
  }, zi = Z.2)

  haz.1 = haz.1.new; haz.2 = haz.2.new; 
  
  d_11 = which((Sur1.1$delta + Sur2.1$delta) == 2)
  d_10 = which((Sur1.1$delta - Sur2.1$delta) == 1)
  d_01 = which((Sur1.1$delta - Sur2.1$delta) == -1)
  d_00 = which((Sur1.1$delta + Sur2.1$delta) == 0)
  
  # log_rho = function(rho)
  # {
  #   k = mapply(function(w.1, z.1, zsur.1, h.1, w.2, z.2, zsur.2, h.2) {
  #     H.1 = as.vector((zsur.1 * z.1) %*% haz.1[1:length(z.1)]) * exp(w.1)  # cummulative hazard
  #     H.2 = as.vector((zsur.2 * z.2) %*% haz.2[1:length(z.2)]) * exp(w.2) # cummulative hazard
  #     if (h.1$delta == 1 && h.2$delta == 1) { # event
  #       (log(haz.1[length(z.1)]) + w.1 + log(z.1[length(z.1)]) + log(zsur.1)) - H.1 + (log(
  #         haz.2[length(z.2)]) + w.2 + log(z.2[length(z.2)]) + log(zsur.2)) - H.2 + log((rho + 1) * (exp(-H.1) * exp(
  #           -H.2)) ^ (-(rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^ (-(2 * rho + 1) / rho))
  #     }
  #     else if (h.1$delta == 1 && h.2$delta == 0) { # event
  #       (log(haz.1[length(z.1)]) + w.1 + log(z.1[length(z.1)]) + log(zsur.1)) - H.1 + log(exp(-H.1) ^ (-(
  #         rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-(rho + 1) / rho))
  #     }
  #     else if (h.1$delta == 0 && h.2$delta == 1) { # event
  #       (log(haz.2[length(z.2)]) + w.2 + log(z.2[length(z.2)]) + log(zsur.2)) - H.2 + log(exp(-H.2) ^ (-(
  #         rho + 1)) * (exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-(rho + 1) / rho))
  #     }
  #     else if (h.1$delta == 0 && h.2$delta == 0) { # event
  #       log((exp(-H.1) ^ (-rho) + exp(-H.2) ^ (-rho) - 1) ^  (-1/rho))
  #     }
  #   },
  #   w.1 = Wtgam.1, z.1 = expZ.1, zsur.1 = expZurb.1, h.1 = Sur1.1.list,
  #   w.2 = Wtgam.2, z.2 = expZ.2, zsur.2 = expZurb.2, h.2 = Sur2.1.list,
  #   SIMPLIFY = FALSE)
  #   return(-sum(unlist(k)))
  # }
  # rho.region = c(1e-4, 30)
  # x.inv <- try(optimize(log_rho, rho.region)$minimum, silent=TRUE)
  # if ('try-error' %in% class(x.inv)) {rho.new <- 10} else {rho.new <- optimize(log_rho, rho.region)$minimum}
  
  
  log_rho <- function(rho) {
    Q_rho(rho,
          Sur1.1 = Sur1.1, Sur2.1 = Sur2.1,
          tj.1 = tj.1, tj.2 = tj.2,
          haz.1 = haz.1, haz.2 = haz.2,
          Wtgam.1 = Wtgam.1, Wtgam.2 = Wtgam.2,
          expZ.1 = expZ.1, expZ.2 = expZ.2,
          time_sur_list_1 = time_sur_list_1,
          time_sur_list_2 = time_sur_list_2,
          b_y = bi.y, pb_y = pb.yt,
          beta1 = beta1, beta2 = beta2)
  }

  rho.region <- c(1e-4, 30)
  rho.new <- tryCatch(optimize(log_rho, rho.region)$minimum, error = function(e) 10)
  
  theta.new = list("Sigma.b" = Sigma.b.new, "alpha" = alpha.new, "sigma2" = sigma2.new,
                   "haz.1" = haz.1.new, "gamma.1" = gamma.1.new,
                   "haz.2" = haz.2.new, "gamma.2" = gamma.2.new, rho = rho.new)
  out = list(theta = theta.new)
  return(out)
}


stepEM_Clay_Two <- function(theta, l, t.1, s.1, t.2, s.2, nMC, rho.region)
{
  # =========================
  # unpack theta
  # =========================
  alpha  <- theta$alpha
  Sigma.b <- theta$Sigma.b
  sigma2 <- theta$sigma2
  gamma.1 <- theta$gamma.1
  gamma.2 <- theta$gamma.2
  haz.1 <- theta$haz.1
  haz.2 <- theta$haz.2
  rho   <- theta$rho
  
  # =========================
  # unpack longitudinal blocks
  # =========================
  ni <- l$ni
  X  <- l$X
  y  <- l$y
  XtX.inv <- l$XtX.inv
  Xty <- l$Xty
  XtZ <- l$XtZ
  Z.lon  <- l$Z.lon
  Z.lont <- l$Z.lont
  
  # =========================
  # unpack survival blocks
  # =========================
  # ---- T1
  W.1 <- t.1$W.1
  Sur1 <- t.1$Sur1
  Sur1.1 <- t.1$Sur1.1
  Sur1.1.list <- t.1$Sur1.1.list
  q.1 <- t.1$q.1
  p.1 <- 1L
  r.1 <- 1L
  nev.1 <- t.1$nev.1
  nev.uniq.1 <- t.1$nev.uniq.1
  tj.1 <- t.1$tj.1
  
  Zdat.sur.1 <- s.1$Zdat.sur.1
  Z.sur.1 <- s.1$Z.sur.1
  Zt.sur.1 <- s.1$Zt.sur.1
  Z.1 <- s.1$Z.1
  Zt.1 <- s.1$Zt.1
  IW.1 <- s.1$IW.1
  
  # ---- T2
  W.2 <- t.2$W.2
  Sur2 <- t.2$Sur2
  Sur2.1 <- t.2$Sur2.1
  Sur2.1.list <- t.2$Sur2.1.list
  q.2 <- t.2$q.2
  p.2 <- 1L
  r.2 <- 1L
  nev.2 <- t.2$nev.2
  nev.uniq.2 <- t.2$nev.uniq.2
  tj.2 <- t.2$tj.2
  
  Zdat.sur.2 <- s.2$Zdat.sur.2
  Z.sur.2 <- s.2$Z.sur.2
  Zt.sur.2 <- s.2$Zt.sur.2
  Z.2 <- s.2$Z.2
  Zt.2 <- s.2$Zt.2
  IW.2 <- s.2$IW.2
  
  # =========================
  # helper: robust gammaUpdate
  # =========================
  gammaUpdate_safe <- function(...) {
    out <- try(gammaUpdate(...), silent = TRUE)
    if (inherits(out, "try-error")) {
      # fallback: no update (prevents blowing up EM)
      return(list(gDelta = rep(0, q.1 + p.1 + 1), ok = FALSE))
    }
    list(gDelta = as.vector(out$gDelta), ok = TRUE)
  }
  
  # =========================
  # Monte Carlo: b | y sampling
  # =========================
  Sigmai.inv <- lapply(ni, function(ii) diag(x = rep(1 / sigma2, ii), ncol = ii))
  
  Dinv <- solve(Sigma.b)
  Ai <- mapply(
    FUN = function(zt, s, z) solve((zt %*% s %*% z) + Dinv),
    z = Z.lon, zt = Z.lont, s = Sigmai.inv,
    SIMPLIFY = FALSE
  )
  
  Mi <- mapply(
    function(a, zt, s, yi, xi) as.vector(a %*% (zt %*% s %*% (yi - xi %*% alpha))),
    a = Ai, zt = Z.lont, s = Sigmai.inv, yi = y, xi = X,
    SIMPLIFY = FALSE
  )
  
  Zq <- randtoolbox::sobol(nMC, dim = ncol(Sigma.b), normal = TRUE, scrambling = 1)
  bi.y <- mapply(
    function(m, a) {
      C <- chol(a)
      matrix(rep(m, nMC), nrow = nMC, byrow = TRUE) + (Zq %*% C)
    },
    m = Mi, a = Ai,
    SIMPLIFY = FALSE
  )
  names(bi.y) <- names(Ai)
  
  # =========================
  # build time lists for (1,t) once
  # =========================
  time_sur_list_1 <- split(Zdat.sur.1$time, Zdat.sur.1$id)
  time_sur_list_2 <- split(Zdat.sur.2$time, Zdat.sur.2$id)
  
  # =========================
  # Cox pieces under current gamma
  # =========================
  Wtgam.1 <- mapply(function(w) as.numeric(w %*% gamma.1[1:q.1]), w = W.1, SIMPLIFY = FALSE)
  expZ.1  <- mapply(function(z) exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)])), z = Z.1, SIMPLIFY = FALSE)
  
  Wtgam.2 <- mapply(function(w) as.numeric(w %*% gamma.2[1:q.2]), w = W.2, SIMPLIFY = FALSE)
  expZ.2  <- mapply(function(z) exp(as.numeric(z %*% gamma.2[(q.2+1):(q.2+p.2)])), z = Z.2, SIMPLIFY = FALSE)
  
  # shared beta (scalar)
  beta1 <- as.numeric(gamma.1[-(1:(q.1 + p.1))])
  beta2 <- as.numeric(gamma.2[-(1:(q.2 + p.2))])
  
  # expWArma inputs
  gam.1.scale <- diag(rep(beta1, r.1), ncol(Sigma.b), ncol(Sigma.b))
  gam.2.scale <- diag(rep(beta2, r.2), ncol(Sigma.b), ncol(Sigma.b))
  
  IZ.sur.1 <- mapply(function(x, y) t(x %*% y), x = IW.1, y = Z.sur.1, SIMPLIFY = FALSE)
  IZ.sur.2 <- mapply(function(x, y) t(x %*% y), x = IW.2, y = Z.sur.2, SIMPLIFY = FALSE)
  
  expZ.surb.1 <- expWArma(IZ.sur.1, bi.y, gam.1.scale, Sur1.1.list)
  expZ.surb.2 <- expWArma(IZ.sur.2, bi.y, gam.2.scale, Sur2.1.list)
  
  # =========================
  # likelihood weights pb_y
  # =========================
  logf_sur <- mapply(function(w1, z1, zsur1, h1, w2, z2, zsur2, h2) {
    
    H1 <- as.vector(t(t(zsur1) * z1) %*% haz.1[1:ncol(zsur1)]) * exp(w1)
    H2 <- as.vector(t(t(zsur2) * z2) %*% haz.2[1:ncol(zsur2)]) * exp(w2)
    
    if (h1$delta == 1 && h2$delta == 1) {
      (log(haz.1[ncol(zsur1)]) + w1 + log(z1[ncol(zsur1)]) + log(zsur1[, ncol(zsur1)])) - H1 +
        (log(haz.2[ncol(zsur2)]) + w2 + log(z2[ncol(zsur2)]) + log(zsur2[, ncol(zsur2)])) - H2 +
        log((rho + 1) * (exp(-H1) * exp(-H2)) ^ (-(rho + 1)) *
              (exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1) ^ (-(2 * rho + 1) / rho))
    } else if (h1$delta == 1 && h2$delta == 0) {
      (log(haz.1[ncol(zsur1)]) + w1 + log(z1[ncol(zsur1)]) + log(zsur1[, ncol(zsur1)])) - H1 +
        log(exp(-H1)^(-(rho + 1)) * (exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1) ^ (-(rho + 1) / rho))
    } else if (h1$delta == 0 && h2$delta == 1) {
      (log(haz.2[ncol(zsur2)]) + w2 + log(z2[ncol(zsur2)]) + log(zsur2[, ncol(zsur2)])) - H2 +
        log(exp(-H2)^(-(rho + 1)) * (exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1) ^ (-(rho + 1) / rho))
    } else {
      log((exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1)^(-1 / rho))
    }
    
  },
  w1 = Wtgam.1, z1 = expZ.1, zsur1 = expZ.surb.1, h1 = Sur1.1.list,
  w2 = Wtgam.2, z2 = expZ.2, zsur2 = expZ.surb.2, h2 = Sur2.1.list,
  SIMPLIFY = FALSE)
  
  f_pdf <- lapply(logf_sur, exp)
  
  den <- lapply(f_pdf, mean)
  for (k in which(unlist(den) == 0)) {
    f_pdf[[k]] <- rep(1, length(f_pdf[[1]]))
    den[[k]] <- 1
  }
  
  pb.yt <- mapply(function(f, d) {
    if (d == 1 || is.nan(d) || is.na(d) || is.infinite(d)) rep(1, length(f)) else f / d
  }, f = f_pdf, d = den, SIMPLIFY = FALSE)
  
  # =========================
  # E-step: Eb, EbbT
  # =========================
  Eb <- mapply(function(b, pb) colMeans(b * pb), b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  EbbT <- mapply(function(b, pb) crossprod(b, (b * pb)) / nrow(b), b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  
  # pieces for gammaUpdate
  expvstargam.1 <- mapply(function(w1, z1, zsur1) exp(w1) * t(t(zsur1) * z1),
                          w1 = Wtgam.1, z1 = expZ.1, zsur1 = expZ.surb.1, SIMPLIFY = FALSE)
  expvstargam.2 <- mapply(function(w2, z2, zsur2) exp(w2) * t(t(zsur2) * z2),
                          w2 = Wtgam.2, z2 = expZ.2, zsur2 = expZ.surb.2, SIMPLIFY = FALSE)
  
  # =========================
  # baseline hazard updates (Breslow) with (1,t)
  # =========================
  haz.hat.1 <- update_haz_breslow_one_t(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ  = expZ.1,
    time_sur_list = time_sur_list_1,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta1
  )
  
  haz.hat.2 <- update_haz_breslow_one_t(
    tj = tj.2,
    Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2,
    expZ  = expZ.2,
    time_sur_list = time_sur_list_2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  # =========================
  # gamma updates (safe, FIXED delta)
  # =========================
  
  P1 <- q.1 + p.1 + 1L   # q + qq + K   (here qq=p.1=1, K=1)
  g1 <- try(gammaUpdate_fixDelta(bi.y, Zt.sur.1, expvstargam.1, pb.yt, haz.hat.1,
                                 W.1, Zt.1, Sur1.1.list, 1, q.1, p.1, nev.uniq.1, nev.1),
            silent = TRUE)
  if (inherits(g1, "try-error") || is.null(g1$gDelta) || any(!is.finite(g1$gDelta))) {
    gDelta.1 <- rep(0, P1)
  } else {
    gDelta.1 <- as.vector(g1$gDelta)
    if (length(gDelta.1) != P1) gDelta.1 <- rep(0, P1)  # 强制防错
  }
  
  P2 <- q.2 + p.2 + 1L
  g2 <- try(gammaUpdate_fixDelta(bi.y, Zt.sur.2, expvstargam.2, pb.yt, haz.hat.2,
                                 W.2, Zt.2, Sur2.1.list, 1, q.2, p.2, nev.uniq.2, nev.2),
            silent = TRUE)
  if (inherits(g2, "try-error") || is.null(g2$gDelta) || any(!is.finite(g2$gDelta))) {
    gDelta.2 <- rep(0, P2)
  } else {
    gDelta.2 <- as.vector(g2$gDelta)
    if (length(gDelta.2) != P2) gDelta.2 <- rep(0, P2)
  }
  
  # =========================
  # M-step
  # =========================
  Sigma.b.new <- Reduce("+", EbbT) / nrow(Sur1.1)
  rownames(Sigma.b.new) <- colnames(Sigma.b.new) <- rownames(Sigma.b)
  
  rr <- mapply(function(x1, x2, b) x1 - (x2 %*% b), x1 = Xty, x2 = XtZ, b = Eb)
  rr.sum <- rowSums(rr)
  alpha.new <- as.vector(XtX.inv %*% rr.sum)
  names(alpha.new) <- names(alpha)
  
  SSq <- mapply(function(yi, xi, zi, b, b2) {
    residFixed <- (yi - xi %*% alpha.new)
    t(residFixed) %*% (residFixed - 2 * (zi %*% b)) + sum(diag(crossprod(zi) %*% b2))
  }, yi = y, xi = X, zi = Z.lon, b = Eb, b2 = EbbT)
  sigma2.new <- sum(SSq) / sum(unlist(ni))
  
  gamma.1.new <- gamma.1 + gDelta.1
  gamma.2.new <- gamma.2 + gDelta.2
  
  # =========================
  # update haz again under new gamma, then update rho under NEW haz/gamma
  # =========================
  gamma.1 <- gamma.1.new
  gamma.2 <- gamma.2.new
  
  Wtgam.1 <- mapply(function(w) as.numeric(w %*% gamma.1[1:q.1]), w = W.1, SIMPLIFY = FALSE)
  expZ.1  <- mapply(function(z) exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)])), z = Z.1, SIMPLIFY = FALSE)
  beta1   <- as.numeric(gamma.1[-(1:(q.1 + p.1))])
  
  Wtgam.2 <- mapply(function(w) as.numeric(w %*% gamma.2[1:q.2]), w = W.2, SIMPLIFY = FALSE)
  expZ.2  <- mapply(function(z) exp(as.numeric(z %*% gamma.2[(q.2+1):(q.2+p.2)])), z = Z.2, SIMPLIFY = FALSE)
  beta2   <- as.numeric(gamma.2[-(1:(q.2 + p.2))])
  
  haz.1.new <- update_haz_breslow_one_t(
    tj = tj.1, Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1, expZ = expZ.1,
    time_sur_list = time_sur_list_1,
    b_y = bi.y, pb_y = pb.yt, beta_shared = beta1
  )
  
  haz.2.new <- update_haz_breslow_one_t(
    tj = tj.2, Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2, expZ = expZ.2,
    time_sur_list = time_sur_list_2,
    b_y = bi.y, pb_y = pb.yt, beta_shared = beta2
  )
  
  # ----- rho update (your Q_rho already returns NEGATIVE sum) -----
  log_rho <- function(rho_tmp) {
    Q_rho(rho_tmp,
          Sur1.1 = Sur1.1, Sur2.1 = Sur2.1,
          tj.1 = tj.1, tj.2 = tj.2,
          haz.1 = haz.1.new, haz.2 = haz.2.new,
          Wtgam.1 = Wtgam.1, Wtgam.2 = Wtgam.2,
          expZ.1 = expZ.1, expZ.2 = expZ.2,
          time_sur_list_1 = time_sur_list_1,
          time_sur_list_2 = time_sur_list_2,
          b_y = bi.y, pb_y = pb.yt,
          beta1 = beta1, beta2 = beta2)
  }
  
  rho.new <- tryCatch(optimize(log_rho, rho.region)$minimum, error = function(e) 10)
  
  theta.new <- list(
    Sigma.b = Sigma.b.new,
    alpha   = alpha.new,
    sigma2  = sigma2.new,
    haz.1   = haz.1.new,
    gamma.1 = gamma.1.new,
    haz.2   = haz.2.new,
    gamma.2 = gamma.2.new,
    rho     = rho.new
  )
  
  list(theta = theta.new)
}
