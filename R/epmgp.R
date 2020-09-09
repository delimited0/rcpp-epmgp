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