# Numerical experiment code 
# originally written by Ted Westling
# modified by Kenta Takatsu

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x) - log(1-x)
qW0 <- function(w) prod(dnorm(w))  # product of all densities by normal distribution 

# lambda is a function of w
lambda <- function(w, b, kappa) c(kappa + 2 * (1 - kappa) * expit(w %*% b))
lambda <- function(w, b, kappa) c(expit(kappa + w %*% b))


# g models conditional density of u given w
g0 <- function(u, lambda) c(lambda + 2 * (1 - lambda) * u)

# inverse of distribution function
G0inv <- function(u, lambda) {
  if(length(lambda) == 1) {
    if(lambda == 1) return(u)
    (-lambda + sqrt(lambda^2 + 4 * u * (1 - lambda)))/(2 * (1 -lambda))
  } else {
    ret <- rep(NA, length(u))
    ones <- lambda == 1
    ret[ones] <- u[ones]
    ret[!ones] <- (-lambda[!ones] + sqrt(lambda[!ones]^2 + 4 * u[!ones] * (1 - lambda[!ones])))/(2 * (1 -lambda[!ones]))
    return(ret)
  }
}

# true mean function mu0(a, w)
mu0 <- function(u, w, gamma1, gamma2, gamma3, gamma4) {
  out <- c(expit(gamma1[1] + w %*% gamma1[-1] + u * (gamma2[1] + w %*% gamma2[-1])+ gamma3 * u^2 + gamma4 * u^3))
  return(out)
}
  
# mixuture of two normal centered at 2 and -2
f0 <- function(x) .5 * dnorm(x, -1) + .5 * dnorm(x, 1) 
F0 <- function(x) .5 * pnorm(x, -1) + .5 * pnorm(x, 1)

# f0 <- function(x) dnorm(x) 
# F0 <- function(x) pnorm(x)

F0.vals <- data.frame(x=seq(-10, 10, by=.01))
F0.vals$F0 <- F0(F0.vals$x)
F0inv <- function(x) approx(F0.vals$F0, F0.vals$x, xout=x, rule=2)$y  #qnormMix(x, mean1 = -2, sd1 = 1, mean2=2, sd2 = 1, p=.5)

# true dose-response curve
theta0 <- function(x, gamma1, gamma2, gamma3, gamma4) {
  sapply(F0(x), function(u0) {
    c <- gamma1[-1] + u0 * gamma2[-1]
    f <- function(z) expit(gamma1[1] + gamma2[1] * u0 + gamma3 * u0^2 + gamma4 * u0^3 + z * sqrt(sum(c^2))) * dnorm(z)
    integrate(f, -Inf, Inf)$value
  })
}



# observational response curve
marg.theta0 <- function(x, gamma1, gamma2, gamma3, gamma4, beta, kappa) {
  sapply(x, function(x0) {
    u0 <- F0(x0)
    c <- gamma1[-1] + u0 * gamma2[-1]
    f <- function(z) expit(gamma1[1] + gamma2[1] * u0 + 
                             gamma3 * u0^2 + gamma4 * u0^3 + 
                             z * sqrt(sum(c^2))) *  
      g0(u0, kappa + 2 * (1 - kappa) * expit(z * sqrt(sum(beta^2))))  * dnorm(z)
    integrate(f, -Inf, Inf)$value
  })
}

# mu.true <- function(x, w, gamma1, gamma2, gamma3, gamma4){
#   mu0(F0(x), w, gamma1, gamma2, gamma3, gamma4)
# }

generate.isotonic.data <- function(n, avals, dw = 4, 
                                  beta0 = c(-1,-1,1,1), 
                                  kappa0 = .1,
                                  gamma1.0 = c(-1, -1,-1,1,1),
                                  gamma2.0 = c(3, -1,-1,1,1),
                                  gamma3.0 = 3,
                                  gamma4.0 = 0){
  W <- matrix(rnorm(n * dw, 0, 1), ncol=dw)  # generate confounder d=4
  Lambda <- lambda(W, b=beta0, kappa=kappa0) # lambda
  U0 <- G0inv(runif(n), Lambda) # generate U
  g0.vals <- g0(U0, Lambda) # convert to density
  g0.true <- function(x, w){
    u <- F0(x)
    w <- as.matrix(w)
    Lambda <- lambda(w, b=beta0, kappa=kappa0) # lambda
    g0(u, Lambda)
  }
  mu0.vals <- mu0(U0, W, gamma1.0, gamma2.0, gamma3.0, gamma4.0) # true mean
  mu0.true <- function(x, w){
    u <- F0(x)
    w <- as.matrix(w)
    mu0(u, w, gamma1.0, gamma2.0, gamma3.0, gamma4.0);
  }
  mu0.vals.adj <- theta0(avals, gamma1.0, gamma2.0, gamma3.0, gamma4.0) # true mean (confounder adjusted)
  X <- F0inv(U0)
  obs <- data.frame(X, Y=mu0.vals, W=W)
  adj <- data.frame(X=avals, Y=mu0.vals.adj)
  list(obs=obs, adj=adj, mu.true=mu0.true, g.true=g0.true, Lambda=lambda, g=g0)
}


syth.data <- generate.isotonic.data(1000, seq(-4,4,by=0.1))
plot(syth.data$adj$X, syth.data$adj$Y)

