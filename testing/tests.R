d <- 1000
mu <- rep(1, d)
  # rep(0, d)
# Sigma <- diag(d)
Sigma <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb <- rep(0, d)
ub <- rep(Inf, d)

epmgp::axisepmgp(mu, Sigma, lb, ub)

epmgp::ptmvn(lb, ub, mu, Sigma)

mu <- c(0, 0)
Sigma <- rep(1, 2)

lb <- c(0, -1)
ub <- c(1, 1)

lb <- c(-Inf, -Inf)
ub <- c(0, 0)

lb <- c(0, 0)
ub <- c(Inf, Inf)
epmgp::trunc_norm_moments(lb, ub, mu, Sigma)

epmgp::erfcx(100)

# testing trunc_norm_moments
mu <- c(0, 0, 0, 0)
Sigma <- rep(1, 4)
lb <- c(0, -Inf, 0, -Inf)
ub <- c(1, 0, Inf, Inf)
epmgp::trunc_norm_moments(lb, ub, mu, Sigma)

# polytope regions ----

# 2d trapezoid
A <- matrix(c(2, -1, 0, 1, -1, 0, 0, -1), byrow = TRUE, ncol = 2)
b <- c(0, 1, 1, 1)
lb <- rep(-Inf, 4)
mu <- rep(0, 2)
Sigma <- diag(2)

epmgp::pmvn2(mu, Sigma, lb, b, A)

