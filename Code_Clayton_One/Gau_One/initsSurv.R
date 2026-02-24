initsSurv <- function(data, Sur1, lfit, Sur1.1, id, timeVar, q) {

  b <- matrix(nrow = nrow(data), ncol = 1)
  XbetaZb <- predict(lfit, level = 1)
  Xbeta <- predict(lfit, level = 0)
  b[,1] <- XbetaZb - Xbeta
  b = data.frame(id = data$id, b = b)

  dataAG <- by(data, data[ , id], FUN = function(u) {
    id.col <- u[ , id]
    b1 = b$b[which((b$id == id.col[1]) > 0)]
    B1 = Sur1$B1[id.col[1]]; B2 = Sur1$B2[id.col[1]]; B3 = Sur1$B3[id.col[1]]
    T <- Sur1.1[Sur1.1$id == id.col[1], "T"] + 1e-06
    start <- u[ , timeVar]
    id.t = max(which((T > start) == 1))
    b1 = b1[1:id.t]
    start = u[1:id.t, timeVar]
    stop <- c(u[1:id.t, timeVar][-1], T)
    status <- rep(0, id.t)
    status[id.t] <- Sur1.1[Sur1.1$id == id.col[1], "delta"]
    if (q > 0) {
      X <- Sur1.1[Sur1.1$id == id.col[1], 2:(q + 1), drop = FALSE]
      X <- X[rep(1, id.t), ]
      if (q == 1) {
        X <- matrix(X, ncol = 1)
      }
      B = B1 * (start <= B3) + B2 * (start > B3)
      X = cbind(X, B)
      #Attention! We add the Z_ik1 = t = start
      colnames(X) <- c(names(Sur1.1)[2:(q + 1)], "time")
    }
    if (q > 0) {
      data.frame("id" = id.col[1:id.t], start, stop, status, X, "beta_1" = b1)
    } else {
      data.frame("id" = id.col[1:id.t], start, stop, status, "beta_1" = b1)
    }
  })
  dataAG <- do.call("rbind", dataAG)
  #dataAG <- cbind(dataAG, b)

  formK <- paste0("beta_", 1:1, collapse = " + ")
  if (q > 0) {
    formX <- paste0(c(names(Sur1.1)[2:(q + 1)], "time"), collapse = " + ")
    formS <- paste("Surv(start, stop, status) ~ ", formX, "+", formK)
  } else {
    formS <- paste("Surv(start, stop, status) ~", formK)
  }
  fitAG <- survival::coxph(as.formula(formS), data = dataAG)

  gamma <- coef(fitAG)
  haz <- survival::coxph.detail(fitAG)$hazard

  return(list("gamma" = gamma, "haz" = haz))

}


initsSurv_new <- function(data, Sur1, lfit, Sur1.1, id, timeVar, q, eps = 1e-6) {
  
  ## --------- shared RE proxy from longitudinal model ----------
  XbetaZb <- predict(lfit, level = 1)
  Xbeta   <- predict(lfit, level = 0)
  b <- data.frame(id = data$id, b = as.numeric(XbetaZb - Xbeta))
  
  ## --------- build AG data (counting process) ----------
  dataAG_list <- by(data, data[, id], FUN = function(u) {
    
    sid <- u[1, id]
    
    ## subject-specific b(t) values (one per longitudinal row)
    b1 <- b$b[b$id == sid]
    
    ## event/censor time and delta
    Ti_row <- Sur1.1[Sur1.1$id == sid, c("T", "delta")]
    if (nrow(Ti_row) != 1) return(NULL)
    T <- as.numeric(Ti_row$T) + eps
    delta_i <- as.integer(Ti_row$delta)
    
    start_all <- u[, timeVar]
    
    ## keep only longitudinal visits strictly before T
    idx <- which(T > start_all)
    if (length(idx) == 0) return(NULL)
    id.t <- max(idx)
    
    ## truncate to visits before T
    b1    <- b1[1:id.t]
    start <- u[1:id.t, timeVar]
    stop  <- c(u[1:id.t, timeVar][-1], T)
    
    status <- integer(id.t)
    status[id.t] <- delta_i
    
    ## time-dependent covariate B(t)
    B1 <- Sur1$B1[sid]; B2 <- Sur1$B2[sid]; B3 <- Sur1$B3[sid]
    B_t <- B1 * (start <= B3) + B2 * (start > B3)
    
    if (q > 0) {
      ## baseline covariates from Sur1.1: columns 2:(q+1)
      X <- Sur1.1[Sur1.1$id == sid, 2:(q + 1), drop = FALSE]
      X <- X[rep(1, id.t), , drop = FALSE]
      
      ## FIX: data.frame -> numeric matrix safely
      X <- as.matrix(X)
      storage.mode(X) <- "double"
      
      X <- cbind(X, B_t = B_t)
      
      ## set colnames robustly
      colnames(X)[1:q] <- names(Sur1.1)[2:(q + 1)]
      
      out <- data.frame(
        id     = rep(sid, id.t),
        start  = start,
        stop   = stop,
        status = status,
        X,
        beta_1 = b1
      )
    } else {
      out <- data.frame(
        id     = rep(sid, id.t),
        start  = start,
        stop   = stop,
        status = status,
        B_t    = B_t,
        beta_1 = b1
      )
    }
    
    out
  })
  
  ## filter NULLs (important!)
  dataAG_list <- Filter(Negate(is.null), as.list(dataAG_list))
  dataAG <- do.call("rbind", dataAG_list)
  if (is.null(dataAG) || nrow(dataAG) == 0) {
    stop("initsSurv_new: dataAG is empty (no subjects have visits before T).")
  }
  
  ## --------- build formulas ----------
  formK <- "beta_1"
  
  if (q > 0) {
    base_covs <- names(Sur1.1)[2:(q + 1)]  # e.g. "W"
    formX_full <- paste(c(base_covs, "B_t"), collapse = " + ")
    formS_full <- paste("Surv(start, stop, status) ~", formX_full, "+", formK)
    
    ## W-only (use the FIRST baseline covariate among base_covs)
    formS_Wonly <- paste("Surv(start, stop, status) ~", base_covs[1])
  } else {
    formS_full  <- paste("Surv(start, stop, status) ~ B_t +", formK)
    formS_Wonly <- NULL
  }
  
  ## --------- Step 1: W-only Cox for stable init ----------
  init_full <- NULL
  if (!is.null(formS_Wonly)) {
    fitW <- try(survival::coxph(as.formula(formS_Wonly), data = dataAG), silent = TRUE)
    if (!inherits(fitW, "try-error")) {
      gammaW0 <- as.numeric(coef(fitW))
      
      ## full model parameter order in dataAG:
      ## base_covs (q of them), B_t, beta_1  -> total q + 2
      init_full <- rep(0, q + 2)
      init_full[1] <- gammaW0
      # others remain 0
    }
  }
  
  ## --------- Step 2: Full Cox to get gamma + hazard jumps ----------
  if (!is.null(init_full)) {
    fitAG <- survival::coxph(as.formula(formS_full), data = dataAG, init = init_full)
  } else {
    fitAG <- survival::coxph(as.formula(formS_full), data = dataAG)
  }
  
  gamma <- coef(fitAG)
  haz   <- survival::coxph.detail(fitAG)$hazard
  
  return(list(gamma = gamma, haz = haz, dataAG = dataAG))
}