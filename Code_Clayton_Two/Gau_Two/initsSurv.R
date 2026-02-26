initsSurv <- function(data, Sur1, lfit, Sur1.1, id, timeVar, q) {
  b <- matrix(nrow = nrow(data), ncol = 1)
  XbetaZb <- predict(lfit, re.form = NULL)
  Xbeta <- predict(lfit, re.form = NA)
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

initsSurv_one_t <- function(data, Sur1, lfit, Sur1.1, id, timeVar, q, eps = 1e-6) {
  
  ## ---- 1) get (b0_i, b1_i) from lmer ----
  re <- lme4::ranef(lfit)[[id]]          # data.frame, rownames are subject ids
  re_id <- rownames(re)
  
  if (!("(Intercept)" %in% colnames(re))) {
    stop("initsSurv_one_t: ranef(lfit)[[id]] has no '(Intercept)'. Your lmer random effects may not match (1+time|id).")
  }
  if (!(timeVar %in% colnames(re))) {
    stop(paste0("initsSurv_one_t: ranef(lfit)[[id]] has no '", timeVar,
                "'. Your lmer random effects must include random slope for timeVar if Z_sur=(1,t)."))
  }
  
  b0_map <- setNames(re[,"(Intercept)"], re_id)
  b1_map <- setNames(re[, timeVar],       re_id)
  
  ## ---- 2) build AG data ----
  dataAG_list <- by(data, data[[id]], FUN = function(u) {
    
    sid <- as.character(u[[id]][1])
    
    Ti <- Sur1.1[Sur1.1$id == u[[id]][1], "T"] + eps
    di <- Sur1.1[Sur1.1$id == u[[id]][1], "delta"]
    
    start_all <- u[[timeVar]]
    keep <- which(start_all < Ti)
    if (length(keep) == 0) return(NULL)
    
    id.t <- max(keep)
    start <- start_all[1:id.t]
    stop  <- c(start_all[2:id.t], Ti)
    
    status <- rep(0, id.t)
    status[id.t] <- di
    
    ## baseline W (q columns), replicated across rows
    if (q > 0) {
      Wmat <- Sur1.1[Sur1.1$id == u[[id]][1], 2:(q + 1), drop = FALSE]
      Wmat <- as.matrix(Wmat)
      Wmat <- Wmat[rep(1, id.t), , drop = FALSE]
      colnames(Wmat) <- names(Sur1.1)[2:(q + 1)]
    }
    
    ## B_t (piecewise)
    B1 <- Sur1$B1[u[[id]][1]]; B2 <- Sur1$B2[u[[id]][1]]; B3 <- Sur1$B3[u[[id]][1]]
    B_t <- B1 * (start <= B3) + B2 * (start > B3)
    
    ## shared term shape: Z_sur(t)^T b_i = b0_i + t * b1_i  (t uses 'start' = left endpoint)
    if (!(sid %in% names(b0_map)) || !(sid %in% names(b1_map))) return(NULL)
    bshare <- as.numeric(b0_map[[sid]] + start * b1_map[[sid]])
    
    if (q > 0) {
      out <- data.frame(
        id = rep(u[[id]][1], id.t),
        start = start, stop = stop, status = status,
        Wmat,
        B_t = B_t,
        beta_shared = bshare
      )
    } else {
      out <- data.frame(
        id = rep(u[[id]][1], id.t),
        start = start, stop = stop, status = status,
        B_t = B_t,
        beta_shared = bshare
      )
    }
    
    out
  })
  
  dataAG <- do.call(rbind, dataAG_list)
  if (is.null(dataAG) || nrow(dataAG) == 0) stop("initsSurv_one_t: dataAG empty.")
  
  ## ---- 3) Cox init in two steps ----
  if (q > 0) {
    formW <- paste(names(Sur1.1)[2:(q + 1)], collapse = " + ")
    form_Wonly <- as.formula(paste0("Surv(start, stop, status) ~ ", formW))
    form_full  <- as.formula(paste0("Surv(start, stop, status) ~ ", formW, " + B_t + beta_shared"))
  } else {
    form_Wonly <- NULL
    form_full  <- as.formula("Surv(start, stop, status) ~ B_t + beta_shared")
  }
  
  init_full <- NULL
  if (!is.null(form_Wonly)) {
    fitW <- try(survival::coxph(form_Wonly, data = dataAG), silent = TRUE)
    if (!inherits(fitW, "try-error")) {
      gammaW0 <- as.numeric(coef(fitW))
      init_full <- rep(0, q + 2)  # q W's + B_t + beta_shared
      init_full[1:q] <- gammaW0
    }
  }
  
  fit <- survival::coxph(form_full, data = dataAG, init = init_full)
  gamma <- coef(fit)
  haz   <- survival::coxph.detail(fit)$hazard
  
  return(list(gamma = gamma, haz = haz))
}