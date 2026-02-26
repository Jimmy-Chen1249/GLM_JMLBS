update_haz_breslow_one_t <- function(
    tj,              # 全局事件时间向量
    Sur1.1,          # 每个subject一行: id, T, delta, ...
    Wtgam,           # list: W_i^T gamma_11（标量）
    expZ,            # list: exp(Z^*(t)^T gamma_12)，按该subject网格时间排列
    time_sur_list,   # list: 该subject网格时间点（递增）
    b_y,             # list: nMC x 2 的 b抽样
    pb_y,            # list: length nMC 的权重
    beta_shared,     # 标量 beta
    eps = 1e-12
) {
  J <- length(tj)
  haz <- numeric(J)
  
  # idx = max{m: time_i[m] <= uu}，用 findInterval（最快且稳）
  find_idx <- function(time_i, uu) {
    idx <- findInterval(uu, time_i)  # 0..length(time_i)
    if (idx < 1) idx <- 1            # 若 uu < 第一个网格点，取1（左连续近似）
    idx
  }
  
  for (j in seq_len(J)) {
    uu <- tj[j]
    
    # 分子：该时刻事件数（无ties=1；有ties=个数）
    dNj <- sum(Sur1.1$delta == 1 & Sur1.1$T == uu)
    
    # 风险集：T_i >= uu
    risk <- which(Sur1.1$T >= uu)
    if (length(risk) == 0) {
      haz[j] <- 0
      next
    }
    
    denom <- 0
    
    for (ii in risk) {
      id_i <- Sur1.1$id[ii]
      
      time_i <- time_sur_list[[id_i]]
      if (is.null(time_i) || length(time_i) == 0) next
      
      idx_i <- find_idx(time_i, uu)
      
      # 标量部分：exp(W^T gamma_11) * exp(Z*(uu)^T gamma_12)
      sc_ij <- exp(Wtgam[[id_i]]) * expZ[[id_i]][idx_i]
      
      # MC：b_i = (b0, b1)，线性项 = b0 + uu*b1
      b_i  <- as.matrix(b_y[[id_i]])      # nMC x 2
      if (ncol(b_i) < 2) stop("b_y[[i]] must be nMC x 2 (b0,b1) for Zsur=(1,t).")
      
      pb_i <- as.numeric(pb_y[[id_i]])    # length nMC
      
      lin <- b_i[, 1] + uu * b_i[, 2]     # nMC
      E_ij <- mean(exp(beta_shared * lin) * pb_i)
      
      denom <- denom + sc_ij * E_ij
    }
    
    denom <- max(denom, eps)
    haz[j] <- dNj / denom
  }
  
  haz
}


## ========= 用于 rho 的 Q 函数：posterior 加权 E[loglik | O] =========
## Z_sur(t) = (1, t)，b 是 2 维，beta 是标量： beta * b^T (1,t)

Q_rho <- function(rho,
                  Sur1.1, Sur2.1,
                  tj.1, tj.2,
                  haz.1, haz.2,
                  Wtgam.1, Wtgam.2,
                  expZ.1, expZ.2,
                  time_sur_list_1, time_sur_list_2,
                  b_y, pb_y,
                  beta1, beta2) {
  
  n <- nrow(Sur1.1)
  
  # Clayton survival copula density c(u,v) with u=S1, v=S2
  clayton_c <- function(u, v, rho) {
    # u,v in (0,1]
    # c(u,v) = (1+rho) (u v)^(-(1+rho)) (u^{-rho}+v^{-rho}-1)^(-(2+1/rho))
    (1 + rho) * (u * v)^(-(1 + rho)) * (u^(-rho) + v^(-rho) - 1)^(-(2 + 1/rho))
  }
  
  out_i <- numeric(n)
  
  for (i in 1:n) {
    
    # subject i 的事件信息
    T1 <- Sur1.1$T[i]; d1 <- Sur1.1$delta[i]
    T2 <- Sur2.1$T[i]; d2 <- Sur2.1$delta[i]
    
    # subject i 的网格时间（与你的 Zdat.sur.* 对齐）
    t1_i <- time_sur_list_1[[as.character(i)]]
    t2_i <- time_sur_list_2[[as.character(i)]]
    
    # 找到 Ti 对应的最后一个网格点 idx
    j1 <- max(which(t1_i <= T1)); if (!is.finite(j1)) j1 <- 1L
    j2 <- max(which(t2_i <= T2)); if (!is.finite(j2)) j2 <- 1L
    
    # baseline cov part
    eta1 <- exp(Wtgam.1[[i]])
    eta2 <- exp(Wtgam.2[[i]])
    
    # exp{Z^*(t) gamma_2}
    ez1 <- expZ.1[[i]]  # length = length(t1_i)
    ez2 <- expZ.2[[i]]  # length = length(t2_i)
    
    # MC samples
    bmat <- as.matrix(b_y[[i]])       # nMC x 2
    pb   <- as.numeric(pb_y[[i]])     # length nMC
    
    # 计算每个 MC 样本下的累计风险 Hk(Tk | b)
    # Hk = sum_{l<=j} haz[l] * eta * ez[l] * exp( beta * b^T (1,t_l) )
    one_t1 <- cbind(1, t1_i[1:j1])     # j1 x 2
    one_t2 <- cbind(1, t2_i[1:j2])     # j2 x 2
    
    # 对每个 MC：lin = beta * b %*% t(one_t) -> nMC x j
    lin1 <- beta1 * (bmat %*% t(one_t1))  # nMC x j1
    lin2 <- beta2 * (bmat %*% t(one_t2))  # nMC x j2
    
    H1m <- as.numeric( (exp(lin1) %*% (haz.1[1:j1] * (eta1 * ez1[1:j1])) ) )
    H2m <- as.numeric( (exp(lin2) %*% (haz.2[1:j2] * (eta2 * ez2[1:j2])) ) )
    
    # S = exp(-H)
    S1m <- exp(-H1m)
    S2m <- exp(-H2m)
    
    # 事件时刻 hazard 需要用 j1, j2 对应点的 exp(beta b^T (1,t))
    lin1_T <- beta1 * (bmat %*% c(1, t1_i[j1]))  # nMC
    lin2_T <- beta2 * (bmat %*% c(1, t2_i[j2]))  # nMC
    
    lam1_T <- haz.1[j1] * eta1 * ez1[j1] * exp(lin1_T)
    lam2_T <- haz.2[j2] * eta2 * ez2[j2] * exp(lin2_T)
    
    # 每个 MC 样本的 loglik（对应 δ 组合）
    # 用 survival-copula 的标准形式：log part = δ1 log λ1 + δ2 log λ2 + log C density + log S terms
    # 这里把 logS = -H 已经包含进 S 的 copula 输入和 δ=1 的 λ*S 结构中时要一致
    # 稳妥写法：
    # joint density for (T1,T2) with survival copula:
    # if δ1=δ2=1: log(λ1) + log(λ2) + log c(S1,S2) + log S1 + log S2
    # if δ1=1,δ2=0: log(λ1) + log c_{10}(S1,S2) ??? 复杂。
    # 你原公式用的是 “四种情况”的 closed form（与你之前一致），我们继续沿用：
    # 但要用 S1,S2 直接写：
    # δ11: log(λ1)+log(λ2)+log c(S1,S2)
    # δ10: log(λ1)+log(∂C/∂u)  (u=S1,v=S2)
    # δ01: log(λ2)+log(∂C/∂v)
    # δ00: log(C)
    # Clayton survival copula C(u,v) = (u^{-rho}+v^{-rho}-1)^(-1/rho)
    
    Cuv <- (S1m^(-rho) + S2m^(-rho) - 1)^(-1/rho)
    
    dC_du <- Cuv^(rho + 1) * S1m^(-rho - 1)   # ∂C/∂u for Clayton
    dC_dv <- Cuv^(rho + 1) * S2m^(-rho - 1)
    
    cden <- clayton_c(S1m, S2m, rho)
    
    loglik_m <- if (d1 == 1 && d2 == 1) {
      log(pmax(lam1_T, 1e-300)) + log(pmax(lam2_T, 1e-300)) + log(pmax(cden, 1e-300))
    } else if (d1 == 1 && d2 == 0) {
      log(pmax(lam1_T, 1e-300)) + log(pmax(dC_du, 1e-300))
    } else if (d1 == 0 && d2 == 1) {
      log(pmax(lam2_T, 1e-300)) + log(pmax(dC_dv, 1e-300))
    } else {
      log(pmax(Cuv, 1e-300))
    }
    
    # posterior 加权期望（注意 pb.yt 的均值是 1，你用 mean(loglik*pb) 就是 E[loglik|O]）
    out_i[i] <- mean(loglik_m * pb)
  }
  
  # optimize() 是最小化，所以返回负号
  -sum(out_i)
}