#################### Code_by Chan (refactored: type1/type2 + save tau only)

rm(list = ls())
setwd("~/Desktop/Chan-Ynu/Sci投稿与返修/CSDA/001_Data_CSDA_R2_2026-02-07/Clayton/Code_robust/")

## ---- packages ----
library(survival)
library(nlme)
library(lme4)
library(Rcpp)
library(MASS)
library(mvtnorm)
library(mnormt)
library(doSNOW)
library(parallel)
library(rootSolve)

## ---- always needed helpers ----
source("initsSurv.R")
source("update_haz_breslow.R")
sourceCpp("gammaUpdate.cpp")
sourceCpp("expW.cpp")
source("stepEM_Clay.R")
source("mjoint_Clay.R")
source("Gdata_Clay.R")
source("stepEM_FGM.R")
source("mjoint_FGM.R")
source("Gdata_FGM.R")
source("stepEM_Gau.R")
source("mjoint_Gau.R")
source("Gdata_Gau.R")

## =========================
## 1) Choose copula (type1) and method (type2)
## =========================
type1 <- "DF"   # DC: Clayton, DF: FGM, DG: Gaussian
type2 <- "MC"   # MC/MF/MG: match with type1 in your design (your original convention)
tau.0 <- 0.15

## ---- map type1 -> data generator + rho.0
if (type1 == "DC") {rho.0 <- 2 * (1/tau.0 - 1)^(-1) }     # Clayton: tau=rho/(rho+2)
if (type1 == "DF") {rho.0 <- 9 * tau.0 / 2 }              # FGM: tau=2rho/9
if (type1 == "DG") {rho.0 <- sin(0.5 * pi * tau.0) }      # Gau: tau=2/pi asin(rho)


## =========================
## 2) Set parameters
## =========================
k1 <- 0.1; k2 <- 5

alpha.0 <- matrix(c(0.5, -1, 1, 1), nrow = 4)
D.0 <- 1
sigma.0e <- 0.5

gamma.011 <- as.matrix(0.3); gamma.012 <- -0.3; beta.01 <- 0.4
gamma.021 <- as.matrix(0.4); gamma.022 <-  0.3; beta.02 <- 0.3

## true parameter vector: last component is tau (NOT rho)
Theta.true <- c(as.numeric(alpha.0), sigma.0e, D.0,
                as.numeric(gamma.011), gamma.012, beta.01,
                as.numeric(gamma.021), gamma.022, beta.02,
                tau.0)

## experiment setup
n <- 200
mtimes <- 20
nMC <- 100
MR <- 0.3
iter_max <- 50
B.boot <- 100
p <- length(Theta.true)

R <- 200
Theta <- NULL; SE <- NULL; CP <- NULL

fileTheta <- paste0(type1, type2, "_Theta_one_n_", n, "_MR_", MR, "_tau_", tau.0, "_R_", R, "_B_", B.boot, ".txt")
fileSE    <- paste0(type1, type2, "_SE_one_n_",    n, "_MR_", MR, "_tau_", tau.0, "_R_", R, "_B_", B.boot, ".txt")
fileCP    <- paste0(type1, type2, "_CP_one_n_",    n, "_MR_", MR, "_tau_", tau.0, "_R_", R, "_B_", B.boot, ".txt")

## =========================
## 3) Parallel setup
## =========================
t.start <- Sys.time()
cl <- makeCluster(getOption("cl.cores", 13))
registerDoSNOW(cl)

pb <- txtProgressBar(max = B.boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## =========================
## 4) Main loop over replications
## =========================
for (r in 1:1) {
  
  set.seed(r)
  
  ## generate data (copula determined by type1)
  if (type1 == "DC") { dat <- Cen_Clay(n, r, MR, rho.0) } 
  if (type1 == "DF") { dat <- Cen_FGM(n, r, MR, rho.0) } 
  if (type1 == "DG") { dat <- Cen_Gau(n, r, MR, rho.0) } 
  
  
  data <- data.frame(id = dat$id, y = dat$y, intercept = dat$X[, 1],
                     time = dat$X[, 2], X1 = dat$X[, 3], X2 = dat$X[, 4])
  
  Sur1 <- data.frame(id = 1:n,
                     surtime = dat$T.cen[, 1], cens = dat$delta[, 1],
                     W = dat$W,
                     B1 = dat$B1[, 1], B2 = dat$B1[, 2], B3 = dat$B1[, 3])
  
  Sur2 <- data.frame(id = 1:n,
                     surtime = dat$T.cen[, 2], cens = dat$delta[, 2],
                     W = dat$W,
                     B1 = dat$B2[, 1], B2 = dat$B2[, 2], B3 = dat$B2[, 3])
  
  Mixeffect <- y ~ time + X1 + X2 + (1 | id)
  formSurv1 <- Surv(surtime, cens) ~ W
  formSurv2 <- Surv(surtime, cens) ~ W
  
  ## fit joint model (estimator determined by type2)
  if (type2 == "MC") {Ree <- try(mjoint_Clay(Mixeffect, data = data, dat = dat,
                                            Sur1, Sur2, formSurv1, formSurv2,
                                            nMC = nMC, iter_max = iter_max),
                                 silent = TRUE)}
  if (type2 == "MF") {Ree <- try(mjoint_FGM(Mixeffect, data = data, dat = dat,
                                            Sur1, Sur2, formSurv1, formSurv2,
                                            nMC = nMC, iter_max = iter_max),
                                 silent = TRUE)}
  if (type2 == "MG") {Ree <- try(mjoint_Gau(Mixeffect, data = data, dat = dat,
                                            Sur1, Sur2, formSurv1, formSurv2,
                                            nMC = nMC, iter_max = iter_max),
                                 silent = TRUE)}
  
  
  if (!inherits(Ree, "try-error") && Ree$iter <= iter_max) {
    
    theta <- Ree$theta
    MR.r  <- 1 - mean(apply(dat$delta, 2, mean))
    
    ## ---- IMPORTANT: save tau, not rho ----
    if (type1 == "MC") {tau.hat <- theta$rho / (2 + theta$rho)}     # Clayton: tau=rho/(rho+2)
    if (type1 == "MF") {tau.hat <- 2 * theta$rho / 9}              # FGM: tau=2rho/9
    if (type1 == "MG") {tau.hat <- 2 * asin(theta$rho) / pi}      # Gau: tau=2/pi asin(rho)
    
    Theta.hatr <- c(theta$alpha,
                    sqrt(theta$sigma2),
                    sqrt(theta$D),
                    theta$gamma.1,
                    theta$gamma.2,
                    tau.hat,            # << store tau here
                    Ree$iter, MR.r, r)
    
    ## =========================
    ## Bootstrap (parallel)
    ## =========================
    #****************Bootstrap***************
    Re = foreach(iter = 1:B.boot,
                 .combine = "rbind",
                 .options.snow = opts,
                 .errorhandling = "pass",
                 .packages = c("survival", "nlme", "lme4", "Rcpp", "MASS",
                               "mvtnorm", "mnormt", "doSNOW", "parallel", "rootSolve"),
                 .noexport = c("expWArma", "gammaUpdate", "hazHat", "lambdaUpdate", "gammaUpdate_fixDelta"),
                 .export = c("dat", "data", "n", "nMC", "iter_max",
                             "Mixeffect", "formSurv1", "formSurv2", "p")) %dopar% {
                               
                               sourceCpp("gammaUpdate.cpp")
                               sourceCpp("expW.cpp")
                               set.seed(iter)
                               
                               boot_ids = sample(1:n, size = n, replace = TRUE)
                               
                               ni_b = dat$ni[boot_ids]
                               id_b = rep(1 : n, ni_b)
                               obse.no_b = unlist(lapply(1 : n, function(i) 1 : ni_b[i]))
                               t_b = lapply(1 : n, function(i) unlist(dat$t[boot_ids[[i]]]))
                               X_b = cbind(1, unlist(t_b), rep(dat$X_al[boot_ids, 1], ni_b), rep(dat$X_al[boot_ids, 2], ni_b))
                               b.long_b = rep(dat$b.surv[boot_ids], ni_b)
                               rho_b = rep(1, length(id_b))
                               y_b = unlist(lapply(1 : n, function(i) {
                                 id_temp = boot_ids[i]
                                 if(id_temp == 1){
                                   dat$y[1:dat$ni[id_temp]]
                                 } else {
                                   dat$y[sum(dat$ni[1:(id_temp-1)]) + (1:dat$ni[id_temp])]
                                 }
                               }))
                               T.cen_b = dat$T.cen[boot_ids, ]
                               W_b = as.matrix(dat$W[boot_ids, ])
                               b.surv_b = dat$b.surv[boot_ids]
                               delta_b = dat$delta[boot_ids, ]
                               B1_b = dat$B1[boot_ids, ]
                               B2_b = dat$B2[boot_ids, ]
                               
                               dat_b = list(id = id_b, ni = ni_b, obse.no = obse.no_b, X = X_b, b.long  = b.long_b, rho = rho_b,
                                            y = y_b, T.cen = T.cen_b, W = W_b, b.surv = b.surv_b,delta = delta_b,
                                            t = t_b, B1 = B1_b, B2 = B2_b)
                               
                               data_b = data.frame(id = dat_b$id, y = dat_b$y, intercept = dat_b$X[, 1],
                                                   time = dat_b$X[, 2], X1 = dat_b$X[, 3], X2 = dat_b$X[, 4])
                               
                               Sur1_b = data.frame(id = 1 : n, surtime = dat_b$T.cen[, 1], cens = dat_b$delta[, 1], W = dat_b$W, 
                                                   B1 = dat_b$B1[, 1], B2 = dat_b$B1[, 2], B3 = dat_b$B1[,3])
                               Sur2_b = data.frame(id = 1 : n, surtime = dat_b$T.cen[, 2], cens = dat_b$delta[, 2], W = dat_b$W, 
                                                   B1 = dat_b$B2[, 1], B2 = dat_b$B2[, 2], B3 = dat_b$B2[,3])
                               
                               Mixeffect = y ~ time + X1 + X2 + (1 | id)
                               formSurv1 = Surv(surtime, cens) ~ W
                               formSurv2 = Surv(surtime, cens) ~ W
                               
                               ## fit joint model (estimator determined by type2)
                               if (type2 == "MC") {fit_b <- try(mjoint_Clay(Mixeffect, data = data_b, dat = dat_b, Sur1_b, Sur2_b, formSurv1, formSurv2,
                                                                          nMC = nMC, iter_max = iter_max), silent = TRUE)}
                               if (type2 == "MF") {fit_b <- try(mjoint_FGM(Mixeffect, data = data_b, dat = dat_b, Sur1_b, Sur2_b, formSurv1, formSurv2,
                                                                          nMC = nMC, iter_max = iter_max), silent = TRUE)}
                               if (type2 == "MG") {fit_b <- try(mjoint_Gau(Mixeffect, data = data_b, dat = dat_b, Sur1_b, Sur2_b, formSurv1, formSurv2,
                                                                          nMC = nMC, iter_max = iter_max), silent = TRUE)}
                               
                               ## 默认：整串 NA，长度固定为 p，保证 rbind 不会崩
                               Theta.hatbr = rep(NA_real_, p)
                               
                               if (!inherits(fit_b, "try-error") && fit_b$iter <= iter_max) {
                                 theta_b = fit_b$theta
                                 MR.r_b  = 1 - mean(apply(dat_b$delta, 2, mean))
                                 ## ---- IMPORTANT: save tau, not rho ----
                                 if (type1 == "MC") {tau.hat_b <- theta_b$rho / (2 + theta_b$rho)}     # Clayton: tau=rho/(rho+2)
                                 if (type1 == "MF") {tau.hat_b <- 2 * theta_b$rho / 9}              # FGM: tau=2rho/9
                                 if (type1 == "MG") {tau.hat_b <- 2 * asin(theta_b$rho) / pi}      # Gau: tau=2/pi asin(rho)
                                 
                                 Theta.hatbr = c(theta_b$alpha, sqrt(theta_b$sigma2), sqrt(theta_b$D),
                                                 theta_b$gamma.1, theta_b$gamma.2,
                                                 tau.hat_b, fit_b$iter, MR.r_b, iter) # / (2 + theta_b$rho)
                               }
                               
                               Theta.hatbr
                             }
    
    ## ---- SE/CP ----
    SE.hatr <- apply(Re, 2, sd, na.rm = TRUE)
    CP.hatr <- (Theta.true >= Theta.hatr[1:p] - 1.96 * SE.hatr[1:p]) *
      (Theta.true <= Theta.hatr[1:p] + 1.96 * SE.hatr[1:p])
    
    Theta <- rbind(Theta, Theta.hatr)
    SE    <- rbind(SE, SE.hatr)
    CP    <- rbind(CP, CP.hatr)
    
    ## ---- running info ----
    CP_param <- CP[, 1:p, drop = FALSE]
    cat("\n====================\n")
    cat("Current replication r =", r, "\n")
    cat("Successful fits so far =", nrow(CP_param), "\n")
    cat("Running mean CP (first p params):\n")
    print(round(colMeans(CP_param, na.rm = TRUE), 4))
    cat("====================\n\n")
    
    ## ---- write results ----
    if (!file.exists(fileTheta)) {
      write.table(t(c(r, Theta.hatr)), file = fileTheta, row.names = FALSE, col.names = TRUE)
    } else {
      write.table(t(c(r, Theta.hatr)), file = fileTheta, append = TRUE, row.names = FALSE, col.names = FALSE)
    }
    
    if (!file.exists(fileSE)) {
      write.table(t(c(r, SE.hatr)), file = fileSE, row.names = FALSE, col.names = TRUE)
    } else {
      write.table(t(c(r, SE.hatr)), file = fileSE, append = TRUE, row.names = FALSE, col.names = FALSE)
    }
    
    if (!file.exists(fileCP)) {
      write.table(t(c(r, CP.hatr)), file = fileCP, row.names = FALSE, col.names = TRUE)
    } else {
      write.table(t(c(r, CP.hatr)), file = fileCP, append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

t.stop <- Sys.time()
print(t.stop - t.start)
stopCluster(cl)