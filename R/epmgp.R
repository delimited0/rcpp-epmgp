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
  result <- epmgp(mu, Sigma, A, lb, ub, max_steps)
  return(result)
}
