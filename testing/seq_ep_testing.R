# trapezoid problem ----

d <- 2
muf <- function(d) rep(0, d)
Sigmaf <- function(d) diag(d)
lbf <- function(d) rep(-Inf, 2*d)
ubf <- function(d) c(0, rep(2, 2*d-1))
Af <- function(d) {
  lower_bounds <- -diag(d)
  upper_bounds <- diag(d)
  upper_bounds[1, ] <- c(2, 1, rep(0, d-2))
  A <- rbind(upper_bounds, lower_bounds)
  return(A)
}

mu <- muf(d)
Sigma <- Sigmaf(d)
lb <- lbf(d)
ub <- ubf(d)
A <- Af(d)

eig_vals <- colSums(A^2)
Ad <- A[, d, drop=FALSE]
Xi <- Ad %*% t(Ad)


Matrix::rankMatrix(Xi)
Xi_eigen <- eigen(Xi)

Xi_inv <- Xi_eigen$vectors[, 1] %*% diag(1 / Xi_eigen$values[1], nrow = 1) %*% 
  t(Xi_eigen$vectors[, 1])

# this is pseduo inverse, check:
Xi_inv %*% Xi %*% Xi_inv
Xi %*% Xi_inv %*% Xi

# now why can't we write this out?

# plug this into seq_epmgp
epmgp::seq_epmgp(mu, Xi_inv, lb, ub, 1)

# what if we just ignore the problem dimension
idx <- c(1, 2, 4)
epmgp::seq_epmgp(mu[idx], Xi_inv[idx, idx], lb[idx], ub[idx], 1)

# compare to SOV
mu_t <- as.numeric(A %*% mu)
lb_t <- lb - mu_t
ub_t <- ub - mu_t
prob <- mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = A %*% Sigma %*% t(A))

epmgp::trunc_norm_moments(lb, ub, rep(0, 4), c(1, 1, 1, Inf))

seq_pmvn(mu = mu, Sigma = Sigma, lb = lb, ub = ub, A = A, max_steps = 2)


trunc_norm_moments(lb, ub, )

# random orthonormal polytope ----
d <- 2
muf <- function(d) rep(0, d)
Sigmaf <- function(d) .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lbf <- function(d) rep(-3, d)
ubf <- function(d) rep(3, d)
Af <- function(d) pracma::randortho(d, "orthonormal")

mu <- muf(d)
Sigma <- Sigmaf(d)
lb <- lbf(d)
ub <- ubf(d)
A <- Af(d)

eig_vals <- colSums(A^2)
Ad <- A[, d, drop=FALSE]
Xi <- Ad %*% t(Ad)

Xi_eigen <- eigen(Xi)
Xi_inv <- Xi_eigen$vectors[, 1] %*% diag(1 / Xi_eigen$values[1], nrow = 1) %*% 
  t(Xi_eigen$vectors[, 1])

result <- epmgp::seq_epmgp(mu, Xi_inv, lb, ub, 200)
exp(result$logZ)

# compare to SOV
mu_t <- as.numeric(A %*% mu)
lb_t <- lb - mu_t
ub_t <- ub - mu_t
prob <- mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = A %*% Sigma %*% t(A))

# identity cross ----
d <- 2
muf <- function(d) rep(0, d)
Sigmaf <- function(d) diag(d)
lbf <- function(d) rep(0, d)
ubf <- function(d) rep(Inf, d)
Abf <- function(d) {
  cross <- volesti::gen_cross(d, "H")
  A <- cross$A[1:(nrow(cross$A)/2), , drop=FALSE]
  return(A)
}

mu <- muf(d)
Sigma <- Sigmaf(d)
lb <- lbf(d)
ub <- ubf(d)
A <- Abf(d)

# eigenvalues of covariances
eig_vals <- colSums(A^2)

# pseudo inverse of x_2 | x_1, needed for EP
Ad <- A[, d, drop=FALSE]
Xi <- Ad %*% t(Ad)
Xi_eigen <- eigen(Xi)
Xi_inv <- Xi_eigen$vectors[, 1] %*% diag(1 / Xi_eigen$values[1], nrow = 1) %*% 
  t(Xi_eigen$vectors[, 1])

# compute (singular) covariance of x_i | x_{i-1}
Xis <- lapply(1:d, function(k) {
  Ad <- A[, k, drop = FALSE]
  Ad %*% t(Ad)  
})

# compute pseudo inverses of x_i | x_{i-1}
Xi_pinvs <- lapply(Xis, function(xi) {
  Xi_eigen <- eigen(xi)
  xi_pinv <- Xi_eigen$vectors[, 1] %*% diag(1 / Xi_eigen$values[1], nrow = 1) %*% 
    t(Xi_eigen$vectors[, 1])
})

# The prior precision
Kinv <- as.matrix(Matrix::bdiag(Xi_pinvs))

# do EP updates
result <- epmgp::seq_epmgp(mu, Xi_inv, lb, ub, 200)

# only the bottom right corner of approximate precision updates
Lambda <- Kinv
Lambda[((d-1)*d + 1):d^2, ((d-1)*d + 1):d^2] <- result$Lambda
Lambda_eigen <- eigen(Lambda)
# Lambda is not full rank, need pseudo determinant
logdet_lambda <- prod(Lambda_eigen$values[1:3])

# the estimate is wrong right now
seq_logprob <- result$logZ - .5 * sum(log(abs(eig_vals))) - .5*logdet_lambda
exp(seq_logprob)

# check the approximating mean
solve(result$Lambda, result$eta)

# compare to SOV
mu_t <- as.numeric(A %*% mu)
lb_t <- lb - mu_t
ub_t <- ub - mu_t

sov_prob <- mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = A %*% Sigma %*% t(A))

axisepmgp_result <- epmgp::moments(mu = mu_t, Sigma = A %*% Sigma %*% t(A), 
                                   lb = lb_t, ub = ub_t)
epmgp_result <- epmgp::moments2(mu, Sigma, lb, ub, A)


