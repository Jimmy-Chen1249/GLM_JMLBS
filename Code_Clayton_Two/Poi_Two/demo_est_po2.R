rm(list = ls())
setwd("~/Desktop/Chan-Ynu/Sci投稿与返修/CSDA/001_Data_CSDA_R2_2026-02-07/Clayton/Code_Clayton_Two/Poi_Two/")
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
n = 200; MR = 0.3; rho.0 = 0.5; R = 50
alpha.0 = matrix(c(0.5, -1, 1, 1), nrow = 4)
Sigma.b0 = matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2);

gamma.011 = as.matrix(0.3); gamma.012 = -0.3; beta.01 = 0.4
gamma.021 = as.matrix(0.4); gamma.022 = 0.3; beta.02 = 0.3

Theta.true <- c(as.numeric(alpha.0), c(Sigma.b0), as.numeric(
  gamma.011), gamma.012, beta.01, as.numeric(gamma.021), gamma.022, beta.02, rho.0)

#*****************************************************************************
mtimes = 20   ##the maximum observed times
nMC = 100; iter_max = 50; p = length(Theta.true)
r.start = 1; Theta = NULL; SE = NULL; CP = NULL

#*************
t.start = Sys.time()
cl = makeCluster(getOption("cl.cores", 13))
envil = environment(Cen_Clay_Two_P)
registerDoSNOW(cl)
iterations = R
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
#*************
r = 1; index = NULL
for(r in 1:R)
{
  tempp = try(read.table(paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_court_", r, ".txt", sep = "")), silent = TRUE)
  if("try-error" %in% class(tempp))
  {
    index = c(index, r)
  }
}
Re = foreach(r = index,
             .combine = "rbind",
             .options.snow = opts,
             .errorhandling = "pass",
             .packages = c("survival", "nlme", "lme4", "Rcpp", "MASS",
                           "mvtnorm", "mnormt", "doSNOW", "parallel", "rootSolve"),
             .noexport = c("expWArma", "gammaUpdate", "hazHat", "lambdaUpdate", "uSamplerPoissonCpp_n", "uSamplerPoissonCpp_n_acc", "gammaUpdate_fixDelta"),
             .export = c("n", "nMC", "iter_max", "p")) %dopar% {
               
               sourceCpp("gammaUpdate.cpp")
               sourceCpp("expW.cpp")
               sourceCpp("uSamplerPoisson_n.cpp")
               
               #****************Estimate***************
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
               Theta.hatr = rep(NA_real_, p)
               if(!("try-error" %in% class(Ree)) && Ree$iter <= iter_max)
               {
                 theta = Ree$theta
                 MR.r = 1 - mean(apply(dat$delta, 2, mean))
                 Theta.hatr = c(theta$alpha, c(theta$Sigma.b), theta$gamma.1, theta$gamma.2, 
                                theta$rho, Ree$iter, MR.r, r)
                 Theta.hatr = unlist(Theta.hatr)
               }
               fileTheta = paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_court_", r, ".txt", sep = "")
               write.table(Theta.hatr, file = fileTheta, row.names = FALSE, col.names = FALSE)
               Theta.hatr
             }

t.stop = Sys.time()
print(t.stop - t.start)
stopCluster(cl)

index1 = NULL
for(r in 1:R)
{
  tempp = try(read.table(paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_court_", r, ".txt", sep = "")), silent = TRUE)
  if("try-error" %in% class(tempp))
  {
    index1 = c(index1, r)
  }
}

Theta.hat = NULL
if(is.null(index1))
{
  for(r in (1:R))
  {
    filename1 = paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_court_", r, ".txt", sep = "")
    tempp = unlist(read.table(filename1))
    Theta.hat = rbind(Theta.hat, tempp)
  }
}

fileTheta = paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_R_", R, ".txt", sep = "")
write.table(Theta.hat, file = fileTheta, row.names = FALSE, col.names = FALSE)

if(nrow(Theta.hat) == R)
{
  for(r in 1:R)
  {
    filename1 = paste("Theta_Poisson_two_n_", n, "_MR_", MR, "_rho_", rho.0, "_court_", r, ".txt", sep = "")
    res <- unlink(filename1)   # 或者 file.remove(filename1)
    if (res != 0) warning(sprintf("删除失败：%s", filename1))
  }
}

apply(Theta.hat, 2, mean)[1:p] - Theta.true