d <- 2
mu <- rep(0, d)
Sigma <- .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb <- rep(0, d)
ub <- rep(Inf, d)

# log normalizer of natural parameterization
np_norm <- -.5*(d * log(2*pi) + determinant(Sigma)$modulus + 
                  t(mu) %*% solve(Sigma) %*% mu)

# normalizer
mom_norm <- -.5*(d * log(2*pi) + determinant(Sigma)$modulus)
