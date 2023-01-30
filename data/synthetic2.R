# bias^2/variance/coverage of bias-corrected method when the nuisances are parametric models
# numerical experiment
expit <- function(x) 1 / (1 + exp(-x))
lambda <- function(w, beta, kappa=0.1) c(kappa + 2 * (1 - kappa) * expit(w %*% beta))
g0 <- function(u, lambda) c(lambda + 2 * (1 - lambda) * u)
G0inv <- function(u, lambda) { 
  ret <- rep(NA, length(u))
  ones <- lambda == 1
  ret[ones] <- u[ones]
  ret[!ones] <- (-lambda[!ones] + sqrt(lambda[!ones]^2 + 4 * u[!ones] * (1 - lambda[!ones])))/
    (2 * (1 -lambda[!ones]))
  ret
}
sin.T <- function(u){
  sin(3*(2*u-1)*pi/2)/(1 + u^2)
}

b <- c(-1, -1, 1, 1)
g1 <- c(-1, -1, -1, 1, 1); g2 <- c(3, -1, -1, 1, 1); g3 <- 1; g4 <- 3

generate.data <- function(n, seed=1234, beta0=b, 
                          gamma.1=g1, gamma.2=g2,
                          gamma.3=g3, gamma.4=g4){
  cat(paste0("setting seed ", seed, "\n"))
  set.seed(seed)
  W <- matrix(rnorm(n * 4, 0, 1), ncol=4)
  U0 <- G0inv(runif(n), lambda(W, beta0))
  mu0 <- function(u, w, gamma1=gamma.1, gamma2=gamma.2, gamma3=gamma.3, gamma4=gamma.4){
    expit(gamma1[1] + w %*% gamma1[-1] + (gamma2[1] + w %*% gamma2[-1]) * u + 
            gamma3 * u^2 + gamma4 *sin.T(u))
  }
  Y <- rbinom(n, 1, mu0(U0, W))
  list(Y=Y, A=U0, W=W)
}

# 
gamma.1 <- c(-1, -1, -1, 1, 1)
gamma.2 <- c(3, -1, -1, 1, 1)
gamma.3 <- 1
gamma.4 <- 3
theta0 <- function(u, gamma1=gamma.1, gamma2=gamma.2, gamma3=gamma.3, gamma4=gamma.4){
  int.f <- function(z){
    c <- gamma1[-1] + u * gamma2[-1]
    expit(gamma1[1] + gamma2[1] * u + gamma3 * u^2 + gamma4 * sin.T(u) +
            z * sqrt(sum(c^2))) * dnorm(z)
  }
  integrate(int.f, -Inf, Inf)$value
}

# library(ggplot2)
# library(latex2exp)
# # 
# As <- seq(0, 1, by=0.005)
# ggplot(data.frame(A=As, theta0=sapply(As, theta0))) +
#   theme_minimal() +
#   theme(text = element_text(size = 15))+
#   geom_line(aes(x=A, y=theta0)) +
#     ylab(label = TeX("$\\theta_0(a)$")) +
#     xlab(label = TeX("$a$"))
# ggsave("~/Documents/workspace/Debiased/debiased_ctseff/fig/truecurve1.png", height=3, width=5)
# 
# theta0.seq <- function(u){
#   sapply(u, theta0)
# }
# library(Deriv)
# theta0.1 <- Deriv(theta0.seq)
# theta0.1("0")
# 
# sm <- smooth.spline(As, sapply(As, theta0))
# sm.2 <- predict(sm,deriv=2)
# ggplot(data.frame(A=As, theta0=sapply(As, theta0),
#                   deriv2=sm.2$y)) +
#   theme_minimal() +
#   theme(text = element_text(size = 15))+
#   geom_line(aes(x=A, y=deriv2))+
#   ylab(label = TeX("$\\theta_0^{(2)}$(a)")) +
#   xlab(label = TeX("$a$")) + geom_hline(yintercept = 0, linetype="dotted")
# ggsave("~/Documents/workspace/Debiased/debiased_ctseff/fig/true_derive.png", height=3, width=5)
# 
# As <- seq(0, 1, length.out=1000)
# ggplot(data.frame(A=As, theta0=sapply(As, theta0))) +
#  theme_minimal() +
#  geom_line(aes(x=A, y=theta0)) +
#    ylab(label = TeX("$\\theta_0(a)$")) +
#    xlab(label = TeX("$a$"))


