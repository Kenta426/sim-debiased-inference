rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("fit_nuisances.R")
source("debiased_ctseff.R")
source("debiased_ctseff_finite.R")
source("my.lprobust.R")
 
library(ggplot2)
library(latex2exp)



data <- read.csv("data/CMR_1990_2010_y.csv")
Y <- data$CMR
A <- data$PM2.5
X <- data[, -c(1,2,3,4)]
X$population_2000 <- log(X$population_2000)

n <- length(A)
bw.seq <- seq(1 * n^{-1/5}, 5 * n^{-1/5}, length.out=50)
# point-wise grid
x0.vals <- seq(3.5, 10, by=0.1)
# uniform grid 
x0.vals.unif <- seq(3.5, 10, length.out=500)
# fit super learner
mu.hat <- fit.regression(Y, A, X, cdf=TRUE)
g.hat <- fit.density(A, X, cdf=TRUE)


est.bc.eff <- debiased.ctseff.finite(Y, A, X, 6, 8, x0.vals, mu.hat, g.hat, 
                     verbose=FALSE, kernel.type="epa")
est.bc.eff$mu2-est.bc.eff$mu1 + 1.96 * est.bc.eff$se.infl.robust
est.bc.eff$mu2-est.bc.eff$mu1 - 1.96 * est.bc.eff$se.infl.robust

est.bc.eff <- debiased.ctseff.finite(Y, A, X, 5, 9, x0.vals, mu.hat, g.hat, 
                                     verbose=FALSE, kernel.type="epa")
est.bc.eff$mu2-est.bc.eff$mu1 + 1.96 * est.bc.eff$se.infl.robust
est.bc.eff$mu2-est.bc.eff$mu1 - 1.96 * est.bc.eff$se.infl.robust

# GP approximation 
est.bc.unif <- unif.debiased.ctseff(Y, A, X, bw.seq, x0.vals, mu.hat, g.hat, 
                                    verbose=FALSE, kernel.type="epa", 
                                    boots.coverage=0.95)

# plotting results
plot.data.y <- data.frame(PM2.5=est.bc.unif$x, CMR=est.bc.unif$mu, 
                        type="Debiased Local Linear", sd=est.bc.unif$se.infl.robust, 
                        unif=approx(est.bc.unif$x, 
                                    est.bc.unif$ep.unif, xout = x0.vals)$y,
                        ylabel="Cardiovascular Mortality Rate  (per 100K)")

data <- read.csv("data/CMR_1990_2010_ratio.csv")
Y <- data$CMR
A <- data$PM2.5
X <- data[, -c(1,2,3,4)]
X$population_2000 <- log(X$population_2000)

# fit super learner
mu.hat <- fit.regression(Y, A, X, cdf=TRUE)
g.hat <- fit.density(A, X, cdf=TRUE)

# GP approximation
est.bc.unif <- unif.debiased.ctseff(Y, A, X, bw.seq, x0.vals, mu.hat, g.hat,
                                    verbose=FALSE, kernel.type="epa",
                                    boots.coverage=0.95)

# plotting results
plot.data.ratio <- data.frame(PM2.5=est.bc.unif$x, CMR=est.bc.unif$mu,
                              type="Debiased Local Linear", sd=est.bc.unif$se.infl.robust,
                              unif=approx(est.bc.unif$x,
                                          est.bc.unif$ep.unif, xout = x0.vals)$y,
                        ylabel="Mortality Rate Ratio")


data <- read.csv("data/CMR_1990_2010_diff.csv")
Y <- data$CMR
A <- data$PM2.5
X <- data[, -c(1,2,3,4)]
X$population_2000 <- log(X$population_2000)

# fit super learner
mu.hat <- fit.regression(Y, A, X, cdf=TRUE)
g.hat <- fit.density(A, X, cdf=TRUE)

# GP approximation
est.bc.unif <- unif.debiased.ctseff(Y, A, X, bw.seq, x0.vals, mu.hat, g.hat,
                                    verbose=FALSE, kernel.type="epa",
                                    boots.coverage=0.95)

# plotting results
plot.data.diff <- data.frame(PM2.5=est.bc.unif$x, CMR=est.bc.unif$mu,
                              type="Debiased Local Linear", sd=est.bc.unif$se.infl.robust,
                              unif=approx(est.bc.unif$x,
                                          est.bc.unif$ep.unif, xout = x0.vals)$y,
                              ylabel="Mortality Rate Difference")

plot.data <- rbind(plot.data.y, plot.data.ratio)
plot.data <- rbind(plot.data, plot.data.diff)
# plot.data <- plot.data.y
ggplot(plot.data, aes(x=PM2.5, y=CMR, ymin=CMR-1.96*sd, ymax=CMR+1.96*sd)) +
  theme_minimal() + 
  geom_line() +
  facet_wrap(~ylabel, 
             strip.position = "left", scales = "free")+
  geom_pointrange(aes(size="95% Pointwise CIs"), fatten=0.01) +
  scale_color_discrete(name="Estimator") +
  geom_line(aes(x=PM2.5, y=CMR+unif, linetype="95% Uniform band")) +
  geom_line(aes(x=PM2.5, y=CMR-unif, linetype="95% Uniform band")) +
  ylab(NULL) +
  xlab(label = TeX("$PM_{2.5} \\,  Concentration \\, (\\mu g/m^3)$")) + 
  scale_linetype_manual("",values=c("95% Uniform band"=2))+
  scale_size_manual("",values=c("95% Pointwise CIs"=0.4))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 11),
        strip.background = element_blank(),
        strip.placement = "outside") 


ggsave("PM25pointwise2.png", width=4, height=4)
