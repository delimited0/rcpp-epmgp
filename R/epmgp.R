
#' @param lb lower bound
#' @param ub upper bound
#' @export
pmvn <- function(lb, ub, mu, Sigma, log = FALSE) {
  result <- axisepmgp(mu, Sigma, lb, ub)
  
  log_prob <- result$logZ
  if (log)
    return(log_prob)
  else
    return(exp(log_prob))
}

#' @param A r x d matrix of constraints
#' @export
pmvn2 <- function(mu, Sigma, lb, ub, A, log = FALSE, max_steps = 200) {
  result <- epmgp(mu, Sigma, t(A), lb, ub, max_steps)
  
  log_prob <- result$logZ
  if (log)
    return(log_prob)
  else
    return(exp(log_prob))
}

#' @export
aug_pmvn <- function(mu, Sigma, lb, ub, A, log = FALSE, max_steps = 200,
                     eig_thresh = 1e-8, moments = FALSE) {
  
  d <- length(mu)
  n_constr <- nrow(A)
  
  # whiten the problem
  R <- chol(Sigma)
  # A_t <- t(R) %*% A
  A_t <- A %*% t(R)
  Sigma_t <- diag(d)
  lb_t <- lb - mu
  ub_t <- ub - mu
  
  # augmented sequence mean, covariance
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
  
  # augmented precision
  # Xi_pinvs <- lapply(1:ncol(A_t), function(d) {
  #   a <- A_t[, d, drop=FALSE]
  #   eig_val <- sum(a^2)
  #   unit_a <- a / sqrt(eig_val)
  #   unit_a %*% t(unit_a) / eig_val
  # })
  # K_inv <- matrix(data = 0, nrow = nrow(A_t)*d, ncol = nrow(A_t)*d)
  # for (i in 1:(d-1)) {
  #   idx_diag <- ((i-1)*n_constr + 1):(i*n_constr)
  #   idx_shift <- (i*n_constr + 1):((i+1)*n_constr)
  #   K_inv[idx_diag, idx_diag] <- Xi_pinvs[[i]] + Xi_pinvs[[i+1]]
  #   K_inv[idx_diag, idx_shift] <- Xi_pinvs[[i+1]]
  #   K_inv[idx_shift, idx_diag] <- Xi_pinvs[[i+1]]
  # }
  # K_inv[idx_shift, idx_shift] <- Xi_pinvs[[i+1]]
  
  # compute eigenvalues of x_i | x_{i-1}
  # eig_vals <- colSums(A_t^2)
  K_eigen <- eigen(K)
  eig_idx <- K_eigen$values > eig_thresh
  eig_vals <- K_eigen$values[eig_idx]
  
  # EP
  result <- aug_epmgp(m, K, lb_t, ub_t, max_steps, d)
  
  # update precision after moment matching
  # Sigma_inv <- K_inv
  # Sigma_inv[idx_shift, idx_shift] <- 
  #   # diag(as.numeric(result$tau_site))
  #   K_inv[idx_shift, idx_shift] + diag(as.numeric(result$tau_site))
  
  # eigenvalues of approximate Sigma
  eig_Sigma <- eigen(result$Sigma, symmetric = TRUE)
  eig_idx <- eig_Sigma$values > eig_thresh
  
  # need a true pseudo inverse of covariance
  Prec <- 
    eig_Sigma$vectors[, eig_idx] %*% 
    diag(1 / eig_Sigma$values[eig_idx]) %*% 
    t(eig_Sigma$vectors[, eig_idx])

  # compute aug log probability
  mu_cavity <- result$nu_cavity / result$tau_cavity
  tau_ratio <- result$tau_site / result$tau_cavity + 1
  
  aug_log_prob <- 
    -.5 * sum(log(eig_vals)) + 
    sum(
      result$logz_hat + 
        .5 * log(tau_ratio) +
        .5 * ( mu_cavity^2 * result$tau_site - 
                 2 * mu_cavity * result$nu_site -
                 (result$nu_site^2 / result$tau_cavity) ) / tau_ratio
    ) +
    .5 * ( t(result$mu) %*% Prec %*% result$mu +
             # sum(result$tau_site * result$mu[idx_shift]^2) +
             sum(log(eig_Sigma$values[eig_idx])) )
    # .5*d * log(2*pi) + .5*sum(log(eig_Sigma$values[eig_idx]))
  
  if (log) {
    ret_prob <- aug_log_prob
  }
  else {
    ret_prob <- exp(aug_log_prob)
  }
  ret_prob <- as.numeric(ret_prob)
  
  if (!moments) {
    return(ret_prob)
  }
  else {
    ret_list <- list(
      "prob" = ret_prob,
      "mu" = result$mu,
      "Sigma" = result$Sigma,
      "Sigma_inv" = Prec
    )
  }
}

#' @export
moments <- function(lb, ub, mu, Sigma) {
  result <- axisepmgp(mu, Sigma, lb, ub)
  return(result)
}

#' @export
moments2 <- function(mu, Sigma, lb, ub, A, max_steps = 200) {
  result <- epmgp(mu, Sigma, t(A), lb, ub, max_steps)
  return(result)
}

#' @export
seq_pmvn <- function(mu, Sigma, lb, ub, A, log=FALSE, max_steps = 200) {
  
  d <- ncol(A)
  eig_vals <- colSums(A^2)
  Ad <- A[, d, drop=FALSE]
  Xi <- Ad %*% t(Ad)
  
  result <- seq_epmgp(mu, Xi, lb, ub, max_steps)
  
  logZ = result$logZ - .5*sum(log(abs(eig_vals)))
  
  if (log)
    return(logZ)
  else
    return(exp(logZ))
}

