rm(list = ls())

library(DebiasedDoseResponse)
library(plyr)
library(here)
library(nprobust)
source("data/synthetic2.R")

run.experiment <- function(n, seed){
  # Generate synthetic data from Section 4 ------------------------------------
  setwd(".")
  b <- c(-1, -1, 1, 1);
  g1 <- c(-1, -1, -1, 1, 1); g2 <- c(3, -1, -1, 1, 1); g3 <- 1; g4 <- 3
  data <- generate.data(n, seed=seed, beta0=b,
                        gamma.1=g1, gamma.2=g2,
                        gamma.3=g3, gamma.4=g4)
  Y <- data$Y; A <- data$A; W <- data$W
  A0.eval <- seq(0.5, 0.9, by=0.01)
  bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 50)
  bw.seq.short <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 20)

  # Fit parametric nuisance estimators ----------------------------------------
  ll <- function(b, u, w) {
    l <- lambda(w, b)
    sum(log(g0(u, l))) # MLE
  }
  # correctly specified density model g
  opt <- optim(par=c(0,0,0,0), fn=ll, u=A, w=W, control=list(fnscale=-1))
  ghat.correct <- function(a, w) g0(a, lambda(as.matrix(w), opt$par))
  # incorrectly specified density model g
  opt.ic <- optim(par=c(0,0), fn=ll, u=A, w=W[,1:2], control=list(fnscale=-1))
  ghat.incorrect <- function(a, w) g0(a, lambda(as.matrix(w[, 1:2]), opt.ic$par))

  # correctly specified outcome regression model mu
  lm.df <- data.frame(Y=Y, A=A, A2=A^2, TA=sin.T(A), W=W)
  mu.correct <- glm(Y ~ W.1 + W.2 + W.3 + W.4 +
                      A * W.1 + A*W.2 + A*W.3 + A*W.4 +
                      A + A2 + TA, data=lm.df, family='binomial')
  muhat.correct <- function(a, w){
    w <- as.matrix(w)
    expit(cbind(1, w, a, a^2, sin.T(a), a*w) %*% coefficients(mu.correct))
  }
  # incorrectly specified outcome regression model mu
  mu.incorrect <- glm(Y ~ W.1 + W.2 +
                        A * W.1 + A*W.2 +
                        A + A2 + TA, data=lm.df, family='binomial')
  muhat.incorrect <- function(a, w){
    w <- as.matrix(w[,1:2])
    expit(cbind(1, w, a, a^2, sin.T(a), a*w) %*% coefficients(mu.incorrect))
  }

  # Iterate over different combination of nuisance estimators -----------------
  save.df <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(save.df) <- c("a1", "a2", "theta.n", "if.se", "t0",
                         "method", "nuisance.mu", "nuisance.g")
  for(muright in c("Correct", "Incorrect")){
    mu.name <- paste0("muhat.", tolower(muright))
    mu.n <- get(mu.name)
    for(gright in c("Correct", "Incorrect")){
      if(muright == "Incorrect" & gright == "Incorrect") {next}
      g.name <- paste0("ghat.", tolower(gright))
      g.n <- get(g.name)

      # Debiased inference with DPI method ------------------------------------
      est.res <- debiased_ate_inference(Y, A, W, mu.n, g.n,
                                        eval.pts.1=A0.eval[-1],
                                        eval.pts.2=A0.eval[1],
                                        bandwidth.method="DPI")
      save.df <- rbind(save.df,
                       data.frame(a1=est.res$eval.pts.1,
                                  a2=est.res$eval.pts.2,
                                  theta.n=est.res$theta.eff,
                                  if.se = est.res$if.val,
                                  t0 = sapply(est.res$eval.pts.1, theta0)-
                                    sapply(est.res$eval.pts.2, theta0),
                                  method="DPI",
                                  nuisance.mu=muright,
                                  nuisance.g=gright))

      # Debiased inference with LOOCV method ----------------------------------
      est.res <- debiased_ate_inference(Y, A, W, mu.n, g.n,
                                        eval.pts.1=A0.eval[-1],
                                        eval.pts.2=A0.eval[1],
                                        bandwidth.method="LOOCV",
                                        bw.seq=bw.seq.short)
      save.df <- rbind(save.df,
                       data.frame(a1=est.res$eval.pts.1,
                                  a2=est.res$eval.pts.2,
                                  theta.n=est.res$theta.eff,
                                  if.se = est.res$if.val,
                                  t0 = sapply(est.res$eval.pts.1, theta0)-
                                    sapply(est.res$eval.pts.2, theta0),
                                  method="LOOCV",
                                  nuisance.mu=muright,
                                  nuisance.g=gright))

      # Debiased inference with LOOCV(h=b) method -----------------------------
      est.res <- debiased_ate_inference(Y, A, W, mu.n, g.n,
                                        eval.pts.1=A0.eval[-1],
                                        eval.pts.2=A0.eval[1],
                                        bandwidth.method="LOOCV(h=b)",
                                        bw.seq=bw.seq)
      save.df <- rbind(save.df,
                       data.frame(a1=est.res$eval.pts.1,
                                  a2=est.res$eval.pts.2,
                                  theta.n=est.res$theta.eff,
                                  if.se = est.res$if.val,
                                  t0 = sapply(est.res$eval.pts.1, theta0)-
                                    sapply(est.res$eval.pts.2, theta0),
                                  method="LOOCV(h=b)",
                                  nuisance.mu=muright,
                                  nuisance.g=gright))
    }
  }
  save.df$n <- n; save.df$seed <- seed
  save.df
}

# Rscript enters here
args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1]); base.seed <- as.numeric(args[2])

# Create directory for saving results
base.dir <- paste0(getwd(), "/result/")
if (!(dir.exists(base.dir))){
  dir.create(base.dir)
}
base.dir <- paste0(getwd(), "/result/effect/")
if (!(dir.exists(base.dir))){
  dir.create(base.dir)
}
if (!(dir.exists(paste0(base.dir, n)))){
  dir.create(paste0(base.dir, n))
}

# Run 100 experiments from seed ~ seed+100
for (i in 1:100){
  out.df <- run.experiment(n, base.seed+i-1)
  save(out.df, file=paste0(base.dir,  n, '/', base.seed, '_', i, '.Rdata'))
}
