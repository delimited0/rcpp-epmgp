#include "tmvn_sampler.h"

// X: d x n matrix of samples
// mu: mean
// R: d x d upper Cholesky factor of covariance
arma::vec lpdf_mvn(const arma::mat & X, const arma::vec & mu, 
                   const arma::mat & R) {
  
  int d = mu.n_elem;
  arma::mat Xmu = X.each_col() - mu;
  arma::mat Q = arma::solve(R.t(), Xmu);
  arma::vec q = arma::sum(Q * Q, 1);
  double c = d * std::log(2 * arma::datum::pi) + 2 * arma::sum(arma::log(arma::diagvec(R)));
  return -(c + q) / 2.;
}

double TmvnSampler::lpdf_tmvn(const arma::vec & x) {
  
  arma::vec eval = F * x + g;
  
  double log_density;
  if (arma::any(eval < 0))
    log_density = -arma::datum::inf;
  else {
    arma::vec xmu = x - this->ep_mean;
    arma::rowvec Q = xmu.t() * this->ep_chol_inv;
    double q = arma::as_scalar(Q * Q.t());
    double c = this->dim * std::log(2. * arma::datum::pi) + 
      2 * arma::sum(arma::log(arma::diagvec(this->ep_chol)));
    log_density = -(c + q) / 2.;
  }
  
  return log_density;
}

