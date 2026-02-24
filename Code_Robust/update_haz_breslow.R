update_haz_breslow <- function(tj, Sur, Sur1.1, Wtgam, expZ, Zsur_list, pb_y, b_y, beta_shared) {
  # tj: 全局事件时间向量 (length = nev.uniq)
  # Sur1.1: 每个 subject 一行，包含 id, T, delta
  # Wtgam[[i]]: 标量 W_i^T gamma_1
  # expZ[[i]]: 向量 exp(Z^*(t) gamma_2)（按该 subject 的网格时间排列）
  # Zsur_list[[i]]: 矩阵/向量，给出 z_sur(t)（按该 subject 的网格时间排列）
  # pb_y[[i]]: nMC 权重向量
  # b_y[[i]]: nMC x 1 的 b 抽样
  # beta_shared: 标量共享随机效应系数
  
  n <- nrow(Sur1.1)
  haz <- numeric(length(tj))
  
  for (j in seq_along(tj)) {
    uu <- tj[j]
    
    # 分子：该时间点事件数（无 ties 时 =1）
    dNj <- sum(Sur1.1$delta == 1 & Sur1.1$T == uu)
    
    # 风险集
    risk <- which(Sur1.1$T >= uu)
    
    denom <- 0
    for (ii in risk) {
      # subject id（你这里 Sur1.1$id 就是 1:n）
      id_i <- Sur1.1$id[ii]
      
      # 取 subject i 自己网格的时间序列（必须和 expZ、Zsur_list 对齐）
      # 如果你没有单独存时间序列，就用 Zsur_list 的行号对应的 time 去匹配
      # 这里建议你在构造时顺便把 time 向量也存成 list，比如 time_sur[[i]]
      # 我这里假设 Zsur_list[[i]] 的行对应的时间就是 1:nrow(...)
      # ——更稳的做法是你传入 time_sur[[i]]，我下面用“索引兜底”写：
      
      m_i <- nrow(as.matrix(Zsur_list[[id_i]]))
      if (m_i == 0) next
      
      # 这里用“最接近且不超过 uu 的最后一个点”作为 idx
      # 如果你能提供该 subject 的网格时间向量 time_i，会更精确
      # 暂时用一个常见假设：expZ[[i]] 与 Zsur_list[[i]] 的行顺序一致，
      # 并且时间是递增的；你应当替换为真实 time_i：
      # idx_i <- max(which(time_i <= uu))
      
      # ---- 你应该用真实的 time_i：比如 time_i <- Zdat.sur.1$time[Zdat.sur.1$id==id_i]
      # 我这里按你已有对象推一个：若你能拿到 time_i，就替换下面两行
      # time_i <- ??? 
      # idx_i <- max(which(time_i <= uu))
      
      # 兜底：若你拿不到 time_i，就用 "min(j, m_i)" 会更安全但没那么严格
      idx_i <- min(j, m_i)
      
      zsur_ij <- as.numeric(as.matrix(Zsur_list[[id_i]])[idx_i, 1])
      ez_ij   <- expZ[[id_i]][idx_i]
      
      b_i  <- as.matrix(b_y[[id_i]])  # nMC x 1
      pb_i <- as.numeric(pb_y[[id_i]])
      
      # E[ exp(beta*b*zsur) | O_i ]，用 pb 权重
      E_ij <- mean( exp(b_i[,1] * beta_shared * zsur_ij) * pb_i )
      
      scalar_ij <- exp(Wtgam[[id_i]]) * ez_ij
      
      denom <- denom + scalar_ij * E_ij
    }
    
    haz[j] <- dNj / denom
  }
  
  haz
}