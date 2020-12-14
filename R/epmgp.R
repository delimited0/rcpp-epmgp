#'
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
pmvn_nce <- function(mu, Sigma, lb, ub, A = diag(length(mu)), data_samples = NULL,
                     n_noise) {
  
  Prec <- chol2inv(chol(Sigma))
  
  if (is.null(data_samples)) {
    data_samples <- t(TruncatedNormal::rtmvnorm(n_noise, mu, Sigma, lb, ub))
  }
  else {
    data_samples <- t(data_samples)
  }
  n_data <- ncol(data_samples)
  log_nu <- log(n_noise / n_data)
  
  ep_approx <- epmgp::moments2(mu, Sigma, lb, ub, A)
  noise_samples <- t(mvtnorm::rmvnorm(n_noise, ep_approx$mu, ep_approx$Sigma))
  
  X_mu <- t(data_samples - mu)
  tmvn_lpdf_data <- -.5 * rowSums( (X_mu %*% Prec) * X_mu )
  mvn_lpdf_data <- mvtnorm::dmvnorm(t(data_samples), ep_approx$mu, 
                                    ep_approx$Sigma, log = TRUE)
  log_ratio_data <- tmvn_lpdf_data - mvn_lpdf_data
  
  X_A <- A %*% noise_samples
  indicators <- apply((X_A >= lb) & (X_A <= ub), 2, all)
  noise_samples <- noise_samples[, indicators]
  X_mu <- t(noise_samples - mu)
  tmvn_lpdf_noise <- -.5 * rowSums( (X_mu %*% Prec) * X_mu)
  mvn_lpdf_noise <- mvtnorm::dmvnorm(t(noise_samples), ep_approx$mu,
                                     ep_approx$Sigma, log = TRUE)
  log_ratio_noise <- tmvn_lpdf_noise - mvn_lpdf_noise
  
  f <- function(c) {
    llik_data <- sum(plogis(log_ratio_data + c - log_nu, log.p = TRUE))
    llik_noise <- sum(plogis(log_ratio_noise + c - log_nu, 
                             lower.tail = FALSE, log.p = TRUE))
    return (- (llik_data + llik_noise) / n_data)
  }
  
  result <- optimize(f, c(-10, 10))
  
  d <- length(mu)
  logprob <- - (result$minimum + .5 *  (d * log(2 * pi) + determinant(Sigma)$modulus))
  prob <- exp(logprob)
  
  
  # prob <- (1 / exp(result$minimum)) / 
  #         ( (2 * pi)^(d / 2) * exp(determinant(Sigma)$modulus)^.5 )
  
  return(prob)
}

