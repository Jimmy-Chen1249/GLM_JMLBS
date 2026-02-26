rm(list = ls())
setwd("~/Desktop/Chan-Ynu/Sci投稿与返修/CSDA/001_Data_CSDA_R2_2026-02-07/Clayton/Code_two/Poi_two/")
library(survival)
library(nlme)
library(lme4)
library(Rcpp)
library(MASS)
library(mvtnorm)
library(mnormt)
library(rootSolve)
library(doSNOW)
library(parallel)
#library(ggplot2)
source("stepEM_Clay_Two_P.R")
source("mjoint_Clay_Two_P.R")
source("Gdata_Clay_Two_P.R")
source("update_haz_breslow_one_t.R")
source("initsSurv.R")
sourceCpp("expW.cpp")
sourceCpp("gammaUpdate.cpp")
sourceCpp("uSamplerPoisson_n.cpp")

#**********************************set parameters********************************************
#set parameters
k1 = 0.1; k2 = 5
n = 200; MR = 0.3; rho.0 = 1; R = 200
alpha.0 = matrix(c(0.5, -1, 1, 1), nrow = 4)
Sigma.b0 = matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2); sigma.0e = 0.5

gamma.011 = as.matrix(0.3); gamma.012 = -0.3; beta.01 = 0.4
gamma.021 = as.matrix(0.4); gamma.022 = 0.3; beta.02 = 0.3

Theta.true <- c(as.numeric(alpha.0), c(Sigma.b0), as.numeric(
  gamma.011), gamma.012, beta.01, as.numeric(gamma.021), gamma.022, beta.02, rho.0)

#*****************************************************************************
mtimes = 20   ##the maximum observed times
nMC = 100; iter_max = 50; B.boot = 100; p = length(Theta.true)
R_tol = R; r.start = 1; Theta = NULL; SE = NULL; CP = NULL

fileTheta = paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_R_", R, "_B_", B.boot, ".txt", sep = "")
fileSE = paste("SE_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_R_", R, "_B_", B.boot, ".txt", sep = "")
fileCP = paste("CP_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_R_", R, "_B_", B.boot, ".txt", sep = "")
#*************
t.start = Sys.time()
cl = makeCluster(getOption("cl.cores", 13))
envil = environment(Cen_Clay_Two_P)
registerDoSNOW(cl)
iterations = B.boot
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
#*************
r = 201
for(r in 176:200)
{
  #****************Estimate***************
  t1 = Sys.time()
  set.seed(r)
  dat = Cen_Clay_Two_P(n, r, MR, rho.0)
  data = data.frame(id = dat$id, y = dat$y, intercept = dat$X[, 1],
                    time = dat$X[, 2], X1 = dat$X[, 3], X2 = dat$X[, 4])
  Sur1 = data.frame(id = 1 : n, surtime = dat$T.cen[, 1], cens = dat$delta[, 1], W = dat$W, 
                    B1 = dat$B1[, 1], B2 = dat$B1[, 2], B3 = dat$B1[,3])
  Sur2 = data.frame(id = 1 : n, surtime = dat$T.cen[, 2], cens = dat$delta[, 2], W = dat$W, 
                    B1 = dat$B2[, 1], B2 = dat$B2[, 2], B3 = dat$B2[,3])
  Mixeffect = y ~ time + X1 + X2 + (1 + time | id)
  formSurv1 = Surv(surtime, cens) ~ W
  formSurv2 = Surv(surtime, cens) ~ W
  Ree = try(mjoint_Clay_Two_P(Mixeffect, data = data, dat = dat, Sur1, Sur2, formSurv1, formSurv2,
                   nMC = nMC, iter_max = iter_max), silent = TRUE)
  t2 = Sys.time()
  print(t2 - t1)
  if(!("try-error" %in% class(Ree)) && Ree$iter <= iter_max)
  {
    theta = Ree$theta
    MR.r = 1 - mean(apply(dat$delta, 2, mean))
    Theta.hatr = c(theta$alpha, c(theta$Sigma.b), theta$gamma.1, theta$gamma.2, 
                   theta$rho, Ree$iter, MR.r, r)
    
    
    #****************Bootstrap***************
    t1 = Sys.time()
    Re = foreach(iter = 1:B.boot,
                 .combine = "rbind",
                 .options.snow = opts,
                 .errorhandling = "pass",
                 .packages = c("survival", "nlme", "lme4", "Rcpp", "MASS",
                               "mvtnorm", "mnormt", "doSNOW", "parallel", "rootSolve"),
                 .noexport = c("expWArma", "gammaUpdate", "hazHat", "lambdaUpdate", "uSamplerPoissonCpp_n", "uSamplerPoissonCpp_n_acc", "gammaUpdate_fixDelta"),
                 .export = c("dat", "data", "n", "nMC", "iter_max",
                             "Mixeffect", "formSurv1", "formSurv2", "p")) %dopar% {
                               
                               sourceCpp("gammaUpdate.cpp")
                               sourceCpp("expW.cpp")
                               sourceCpp("uSamplerPoisson_n.cpp")
                               set.seed(iter)
                               
                               boot_ids = sample(1:n, size = n, replace = TRUE)
                               
                               ni_b = dat$ni[boot_ids]
                               id_b = rep(1 : n, ni_b)
                               obse.no_b = unlist(lapply(1 : n, function(i) 1 : ni_b[i]))
                               t_b = lapply(1 : n, function(i) unlist(dat$t[boot_ids[[i]]]))
                               X_b = cbind(1, unlist(t_b), rep(dat$X_al[boot_ids, 1], ni_b), rep(dat$X_al[boot_ids, 2], ni_b))
                               b.long_b = rep(dat$b.surv[boot_ids], ni_b)
                               rho_b = rep(1, length(id_b))
                               y_b = unlist(lapply(1 : n, function(i) {id_temp = boot_ids[i]; if(id_temp == 1){dat$y[1:dat$ni[id_temp]]}else{
                                 dat$y[sum(dat$ni[1:(id_temp-1)]) + (1:dat$ni[id_temp])]
                               }}))
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
                               #eps = 1e-4
                               Sur1_b = data.frame(id = 1 : n, surtime = dat_b$T.cen[, 1], cens = dat_b$delta[, 1], W = dat_b$W, 
                                                   B1 = dat_b$B1[, 1], B2 = dat_b$B1[, 2], B3 = dat_b$B1[,3])
                               Sur2_b = data.frame(id = 1 : n, surtime = dat_b$T.cen[, 2], cens = dat_b$delta[, 2], W = dat_b$W, 
                                                   B1 = dat_b$B2[, 1], B2 = dat_b$B2[, 2], B3 = dat_b$B2[,3])
                               
                               Mixeffect = y ~ time + X1 + X2 + (1 + time | id)
                               formSurv1 = Surv(surtime, cens) ~ W
                               formSurv2 = Surv(surtime, cens) ~ W
                               fit_b = try(mjoint_Clay_Two_P(Mixeffect, data = data_b, dat = dat_b, Sur1_b, Sur2_b, formSurv1, formSurv2,
                                                  nMC = nMC, iter_max = iter_max), silent = TRUE)
                               
                               ## 默认：整串 NA，长度固定为 p，保证 rbind 不会崩
                               Theta.hatbr = rep(NA_real_, p)
                               
                               if (!inherits(fit_b, "try-error") && fit_b$iter <= iter_max) {
                                 theta_b = fit_b$theta
                                 MR.r_b  = 1 - mean(apply(dat_b$delta, 2, mean))
                                 Theta.hatbr = c(theta_b$alpha, c(theta_b$Sigma.b),
                                                 theta_b$gamma.1, theta_b$gamma.2,
                                                 theta_b$rho, fit_b$iter, MR.r_b, iter)
                                 Theta.hatbr = unlist(Theta.hatbr)
                               }
                               
                               Theta.hatbr   # ★ 最后一行作为返回值，不用 return()
                               # write.csv(iter, paste("iter_", iter, ".csv", sep = ""))
                               # Theta.hatbr
                             }
    t2 = Sys.time()
    print(t2 - t1)
    SE.hatr = apply(Re, 2, sd, na.rm = TRUE)
    CP.hatr = (Theta.true >= Theta.hatr[1:p] - 1.96 * SE.hatr[1:p]) * (Theta.true <= Theta.hatr[1:p] + 1.96 * SE.hatr[1:p])
    
    Theta = rbind(Theta, Theta.hatr)
    SE = rbind(SE, SE.hatr)
    CP = rbind(CP, CP.hatr)
    
    # ================== 新增：输出当前前“成功次数”的平均 CP ==================
    CP_param <- CP[, 1:p, drop = FALSE]
    CP_zero_each <- apply(CP_param, 2, function(x) sum(x == 0, na.rm = TRUE))
    CP_mean_each <- colMeans(CP_param, na.rm = TRUE)
    
    cat("\n====================\n")
    cat("Current replication r =", r, "\n")
    cat("Successful fits so far =", nrow(CP_param), "\n")
    cat("Count of CP==0 (each param):\n")
    print(CP_zero_each)
    cat("Running mean CP (each param):\n")
    print(round(CP_mean_each, 4))
    cat("====================\n\n")
    # ======================================================================
    
    # #******************
    ## 写 Theta
    if (!file.exists(fileTheta)) {
      write.table(t(c(r, Theta.hatr)), file = fileTheta, row.names = FALSE, col.names = TRUE)
    } else {
      write.table(t(c(r, Theta.hatr)), file = fileTheta, append = TRUE, row.names = FALSE, col.names = FALSE)
    }

    ## 写 SE
    if (!file.exists(fileSE)) {
      write.table(t(c(r, SE.hatr)), file = fileSE, row.names = FALSE, col.names = TRUE)
    } else {
      write.table(t(c(r, SE.hatr)), file = fileSE, append = TRUE, row.names = FALSE, col.names = FALSE)
    }

    ## 写 CP
    if (!file.exists(fileCP)) {
      write.table(t(c(r, CP.hatr)), file = fileCP, row.names = FALSE, col.names = TRUE)
    } else {
      write.table(t(c(r, CP.hatr)), file = fileCP, append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

t.stop = Sys.time()
print(t.stop - t.start)
stopCluster(cl)