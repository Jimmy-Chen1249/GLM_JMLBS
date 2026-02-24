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
