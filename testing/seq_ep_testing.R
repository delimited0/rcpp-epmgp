# 2d identity positive orthant ----
mu <- rep(0, 2)
Sigma <- diag(2)
lb <- rep(0, 2)
ub <- rep(Inf, 2)
A <- diag(2)

# compute the covariance, which we can
Xi_list <- lapply(split(A, 1:ncol(A)), function(a) a %*% t(a))
K <- matrix(data = 0, nrow = nrow(A)*d, ncol = nrow(A)*d)
for (i in 1:length(Xi_list)) {
  idx <- (1 + (i-1)*nrow(A)):nrow(K)
  K[idx, idx] <- K[idx, idx] + kronecker(matrix(1, d - i + 1, d - i + 1), Xi_list[[i]])
}

m <- rep(0, nrow(K))
aug_result1 <- epmgp::aug_epmgp(m, K, A, lb, ub, 200)

axis_result <- epmgp::axisepmgp(mu, Sigma, lb, ub)
seq_result <- epmgp::seq_epmgp(mu, solve(Sigma), lb, ub, 200)

axis_result$Sigma
solve(seq_result$Lambda)

axis_result$mu
solve(seq_result$Lambda, seq_result$eta)

poly_prob <- epmgp::epmgp(mu, diag(2), A, lb, ub, 200)

aug_prob <- epmgp::aug_pmvn(mu, Sigma, lb, ub, A, moments = TRUE)

# 2d dense positive orthant ----
d <- 2
mu <- rep(0, 2)
Sigma <- .5 * diag(2) + 0.5 * rep(1, 2) %*% t(rep(1, 2))
lb <- rep(0, 2)
ub <- rep(Inf, 2)
A <- diag(2)
n_constr <- nrow(A)

# whiten the problem
R <- chol(Sigma) 
A_t <- t(R) %*% A
Sigma_t <- diag(d)

# compute the covariance, which we can
Xi_list <- lapply(1:ncol(A_t), function(d) {
    a <- A_t[, d, drop=FALSE]
    a %*% t(a)
  })
K <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
for (i in 1:length(Xi_list)) {
  idx <- (1 + (i-1)*nrow(A)):nrow(K)
  K[idx, idx] <- K[idx, idx] + kronecker(matrix(1, d - i + 1, d - i + 1), Xi_list[[i]])
}
m <- rep(0, nrow(K))

# also need the precision
Xi_pinvs <- lapply(1:ncol(A_t), function(d) {
  a <- A_t[, d, drop=FALSE]
  eig_val <- sum(a^2)
  unit_a <- a / sqrt(eig_val)
  unit_a %*% t(unit_a) / eig_val
})
K_inv <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
for (i in 1:(d-1)) {
  idx_diag <- ((i-1)*n_constr + 1):(i*n_constr)
  idx_shift <- (i*n_constr + 1):((i+1)*n_constr)
  K_inv[idx_diag, idx_diag] <- Xi_pinvs[[i]] + Xi_pinvs[[i+1]]
  K_inv[idx_diag, idx_shift] <- Xi_pinvs[[i+1]]
  K_inv[idx_shift, idx_diag] <- Xi_pinvs[[i+1]]
}
K_inv[idx_shift, idx_shift] <- Xi_pinvs[[i+1]]

# but what we get isn't actually a g-inverse of K
# do eigen decomposition to get g inverse
K_eigen <- eigen(K)
K_inv <- K_eigen$vectors[, 1:2] %*% diag(1 / K_eigen$values[1:2]) %*% t(K_eigen$vectors[, 1:2])

# axis_result <- epmgp::axisepmgp(mu, Sigma, lb, ub)
poly_result <- epmgp::epmgp(mu, Sigma, A, lb, ub, 200)
# whiten_result <- epmgp::epmgp(mu, Sigma_t, R, lb, ub, 200)
aug_result <- epmgp::aug_epmgp(m, K, lb, ub, 200, d)

# compute eigenvalues of x_i | x_{i-1}
eig_vals <- colSums(A_t^2)
eig_vals <- K_eigen$values[1:2]

# eigenvalues of approximate Sigma
eig_Sigma <- eigen(aug_result$Sigma, symmetric = TRUE)

# approximate precision
Prec <- eig_Sigma$vectors[, 1:2] %*% diag(1 / eig_Sigma$values[1:2]) %*% t(eig_Sigma$vectors[, 1:2])

# compute aug probability
mu_cavity <- aug_result$nu_cavity / aug_result$tau_cavity
tau_ratio <- aug_result$tau_site / aug_result$tau_cavity + 1

aug_log_prob <- 
  -.5 * sum(log(eig_vals)) +
  sum(
    aug_result$logz_hat + 
      .5 * log(tau_ratio) +
      .5 * ( mu_cavity^2 * aug_result$tau_site - 
             2 * mu_cavity * aug_result$nu_site -
             (aug_result$nu_site^2 / aug_result$tau_cavity) ) / tau_ratio
  ) +
  .5 * ( t(aug_result$mu) %*% Prec %*% aug_result$mu + 
           sum(log(eig_Sigma$values[1:2])) )

# compare relative error 
aug_rele <- (exp(aug_log_prob) - (1/3)) / (1/3)
poly_rele <- (exp(poly_result$logZ) - (1/3)) / (1/3)

aug_prob <- epmgp::aug_pmvn(mu, Sigma, lb, ub, A, moments = TRUE)

# 2d negative covariance ----
d <- 2
mu <- rep(0, 2)
Sigma <- matrix(c(1, -.5, -.5, 1), nrow = 2)
lb <- rep(0, 2)
ub <- rep(Inf, 2)
A <- diag(2)
n_constr <- nrow(A)

# whiten the problem
R <- chol(Sigma) 
A_t <- t(R) %*% A
Sigma_t <- diag(d)

# compute the covariance
Xi_list <- lapply(1:ncol(A_t), function(d) {
  a <- A_t[, d, drop=FALSE]
  a %*% t(a)
})
K <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
for (i in 1:length(Xi_list)) {
  idx <- (1 + (i-1)*nrow(A)):nrow(K)
  K[idx, idx] <- K[idx, idx] + kronecker(matrix(1, d - i + 1, d - i + 1), Xi_list[[i]])
}
m <- rep(0, nrow(K))

# also need the precision
Xi_pinvs <- lapply(1:ncol(A_t), function(d) {
  a <- A_t[, d, drop=FALSE]
  eig_val <- sum(a^2)
  unit_a <- a / sqrt(eig_val)
  unit_a %*% t(unit_a) / eig_val
})
K_inv <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
for (i in 1:(d-1)) {
  idx_diag <- ((i-1)*n_constr + 1):(i*n_constr)
  idx_shift <- (i*n_constr + 1):((i+1)*n_constr)
  K_inv[idx_diag, idx_diag] <- Xi_pinvs[[i]] + Xi_pinvs[[i+1]]
  K_inv[idx_diag, idx_shift] <- Xi_pinvs[[i+1]]
  K_inv[idx_shift, idx_diag] <- Xi_pinvs[[i+1]]
}
K_inv[idx_shift, idx_shift] <- Xi_pinvs[[i+1]]

# do the moment matching
aug_result <- epmgp::aug_epmgp(m, K, lb, ub, 200, d)

# update precision after moment matching
K_inv_ep <- K_inv
K_inv_ep[idx_shift, idx_shift] <- 
  # diag(as.numeric(aug_result$tau_site))
  K_inv_ep[idx_shift, idx_shift] + diag(as.numeric(aug_result$tau_site))

# mean precision after updating
eta <- m
eta[idx_shift] <- aug_result$nu_site 

# compute eigenvalues of x_i | x_{i-1}
eig_vals <- colSums(A_t^2)

# eigenvalues of approximate Sigma
eig_Sigma <- eigen(aug_result$Sigma, symmetric = TRUE)
eig_Lambda <- eigen(K_inv_ep, symmetric = TRUE)

# compute aug probability
mu_cavity <- aug_result$nu_cavity / aug_result$tau_cavity
tau_ratio <- aug_result$tau_site / aug_result$tau_cavity + 1

aug_log_prob <- 
  -.5 * sum(log(eig_vals)) + 
  sum(
    aug_result$logz_hat + 
      .5 * log(tau_ratio) +
      .5 * ( mu_cavity^2 * aug_result$tau_site - 
               2 * mu_cavity * aug_result$nu_site -
               (aug_result$nu_site / aug_result$tau_cavity)^2 ) / tau_ratio
  ) + 
  .5 * ( t(aug_result$mu) %*% K_inv_ep %*% aug_result$mu +
           sum(log(eig_Sigma$values[1:2])) )

poly_result <- epmgp::epmgp(mu, Sigma, A, lb, ub, 200)
sov_result <- mvtnorm::pmvnorm(lb, ub, mu, Sigma)
augm_result <- epmgp::aug_pmvn(mu, Sigma, lb, ub, A, moments = TRUE)

# 3d dense positive orthant ----
d <- 3
mu <- rep(0, d)
Sigma <- .5 * diag(d) + 0.5 * rep(1, d) %*% t(rep(1, d))
lb <- rep(0, d)
ub <- rep(Inf, d)
A <- diag(d)
n_constr <- nrow(A)

# whiten the problem
R <- chol(Sigma) 
A_t <- t(R) %*% A
Sigma_t <- diag(d)

# compute the covariance, which we can
Xi_list <- lapply(1:ncol(A_t), function(d) {
  a <- A_t[, d, drop=FALSE]
  a %*% t(a)
})
K <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
for (i in 1:length(Xi_list)) {
  idx <- (1 + (i-1)*nrow(A)):nrow(K)
  K[idx, idx] <- K[idx, idx] + kronecker(matrix(1, d - i + 1, d - i + 1), Xi_list[[i]])
}
m <- rep(0, nrow(K))

# also need the precision
Xi_pinvs <- lapply(1:ncol(A_t), function(d) {
  a <- A_t[, d, drop=FALSE]
  eig_val <- sum(a^2)
  unit_a <- a / sqrt(eig_val)
  unit_a %*% t(unit_a) / eig_val
})
K_inv <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
for (i in 1:(d-1)) {
  idx_diag <- ((i-1)*n_constr + 1):(i*n_constr)
  idx_shift <- (i*n_constr + 1):((i+1)*n_constr)
  K_inv[idx_diag, idx_diag] <- Xi_pinvs[[i]] + Xi_pinvs[[i+1]]
  K_inv[idx_diag, idx_shift] <- Xi_pinvs[[i+1]]
  K_inv[idx_shift, idx_diag] <- Xi_pinvs[[i+1]]
}
K_inv[idx_shift, idx_shift] <- Xi_pinvs[[i+1]]

aug_result <- epmgp::aug_epmgp(m, K, lb, ub, 200, d)
poly_result <- epmgp::epmgp(mu, Sigma, A, lb, ub, 200)
augm_result <- epmgp::aug_pmvn(mu, Sigma, lb, ub, A, moments = TRUE)

# compute eigenvalues of x_i | x_{i-1}
eig_vals <- colSums(A_t^2)

# eigenvalues of approximate Sigma
eig_Sigma <- eigen(aug_result$Sigma, symmetric = TRUE)

# compute aug probability
mu_cavity <- aug_result$nu_cavity / aug_result$tau_cavity
tau_ratio <- aug_result$tau_site / aug_result$tau_cavity + 1

aug_log_prob <- 
  -.5 * sum(log(eig_vals)) + 
  sum(
    aug_result$logz_hat + 
      .5 * log(tau_ratio) +
      .5 * ( mu_cavity^2 * aug_result$tau_site - 
               2 * mu_cavity * aug_result$nu_site -
               (aug_result$nu_site / aug_result$tau_cavity)^2 ) / tau_ratio
  ) +
  .5 * ( t(aug_result$mu[1:9]) %*% 
           aug_result$Sigma[1:9, 1:9] %*%
           aug_result$mu[1:9] + 
           sum(log(eig_Sigma$values[1:3])) )

# normalizing constant of unconstrained augmented normal



# dense positive orthant comparisons for various dimensions ----
dims <- c(2, 5, 10, 50, 100)
sapply(dims, function(d) {
  
  mu <- rep(0, d)
  Sigma <- .5 * diag(d) + 0.5 * rep(1, d) %*% t(rep(1, d))
  lb <- rep(0, d)
  ub <- rep(Inf, d)
  A <- diag(d)
  
  # whiten the problem
  R <- chol(Sigma) 
  A_t <- A %*% t(R)
  Sigma_t <- diag(d)
  
  polyprob <- epmgp::moments2(mu, Sigma_t, lb, ub, A_t)
  augprob <- epmgp::aug_pmvn(mu, Sigma, lb, ub, A, moments = TRUE)
  
  actual <- 1 / (d+1)
  aug_rele <- (exp(augprob) - actual) / actual
  poly_rele <- (exp(polyprob$logZ) - actual) / actual
})

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

# whiten the problem (in this case already whitened)
R <- chol(Sigma) 
A_t <- A %*% t(R)
Sigma_t <- diag(d)

# compute the covariance, which we can
Xi_list <- lapply(1:ncol(A_t), function(d) {
  a <- A_t[, d, drop=FALSE]
  a %*% t(a)
})
K <- matrix(data = 0, nrow = nrow(A)*d, ncol = nrow(A)*d)
for (i in 1:length(Xi_list)) {
  idx <- (1 + (i-1)*nrow(A)):nrow(K)
  K[idx, idx] <- K[idx, idx] + kronecker(matrix(1, d - i + 1, d - i + 1), Xi_list[[i]])
}
m <- rep(0, nrow(K))

# compute eigenvalues of x_i | x_{i-1}
eig_vals <- colSums(A_t^2)

# call epmgp
aug_result <- epmgp::aug_epmgp(m, K, lb, ub, 200, p = d)

# eigenvalues of approximate Sigma
eig_Sigma <- eigen(aug_result$Sigma, symmetric = TRUE)

# compute aug probability
mu_cavity <- aug_result$nu_cavity / aug_result$tau_cavity
tau_ratio <- aug_result$tau_site / aug_result$tau_cavity + 1

aug_prob <- 
  -.5 * sum(log(eig_vals)) + 
  sum(
    aug_result$logz_hat + 
      .5 * log(tau_ratio) +
      .5 * ( mu_cavity^2 * aug_result$tau_site - 
               2 * mu_cavity * aug_result$nu_site -
               (aug_result$nu_site / aug_result$tau_cavity)^2 ) / tau_ratio
  ) +
  .5 * ( t(aug_result$mu) %*% aug_result$Sigma %*% aug_result$mu +
           sum(log(eig_Sigma$values[1:2])) )

# compare with others
mu_t <- as.numeric(A %*% mu)
lb_t <- lb - mu_t
ub_t <- ub - mu_t
sov_prob <- mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = A %*% Sigma %*% t(A))
aug_prob <- epmgp::aug_pmvn(mu, Sigma, lb, ub, A, moments = TRUE)
poly_prob <- epmgp::moments2(mu, Sigma, lb, ub, A)

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
ubf <- function(d) rep(Inf,  d)
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
# and their eigenvalues
eig_vals <- colSums(A^2)

# compute pseudo inverses of covariance of x_i | x_{i-1}
Xi_pinvs <- lapply(Xis, function(xi) {
  Xi_eigen <- eigen(xi)
  xi_pinv <- Xi_eigen$vectors[, 1] %*% diag(1 / Xi_eigen$values[1], nrow = 1) %*% 
    t(Xi_eigen$vectors[, 1])
})

# The prior precision
Kdiag <- list(Xi_pinvs[[2]], matrix(0, nrow = d, ncol = d))
Kinv <- as.matrix(Matrix::bdiag(Xi_pinvs)) + as.matrix(Matrix::bdiag(Kdiag))

# do EP updates
result <- epmgp::seq_epmgp(mu, Xi_inv, lb, ub, 200)
# result_poly <- epmgp::epmgp(mu, Xi_inv, )

# only the bottom right corner of approximate precision updates
Lambda <- Kinv
Lambda[((d-1)*d + 1):d^2, ((d-1)*d + 1):d^2] <- result$Lambda
Lambda_eigen <- eigen(Lambda)
# Lambda is not full rank, need pseudo determinant
logdet_lambda <- sum(log(Lambda_eigen$values[1:3]))

# the estimate is wrong right now
seq_logprob <- result$logZ - .5 * sum(log(abs(eig_vals))) + .5*logdet_lambda
exp(seq_logprob)

# check the approximating mean
solve(result$Lambda, result$eta)

# check approximating covariance
solve(result$Lambda)

# compare to SOV
mu_t <- as.numeric(A %*% mu)
lb_t <- lb - mu_t
ub_t <- ub - mu_t

sov_prob <- mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = A %*% Sigma %*% t(A))

# compare to the usual EPs
axisepmgp_result <- epmgp::moments(mu = mu_t, Sigma = A %*% Sigma %*% t(A), 
                                   lb = lb_t, ub = ub_t)
epmgp_result <- epmgp::moments2(mu, Sigma, lb, ub, A)

# understanding updates ----
d <- 2
mu <- rep(0, d)
Sigma <- .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb <- rep(0, d)
ub <- rep(Inf, d)
A <- diag(d)

axis_result <- epmgp::moments(lb, ub, mu, Sigma)
poly_result <- epmgp::moments2(mu, Sigma, lb, ub, A)


