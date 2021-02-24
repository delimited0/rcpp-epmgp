#include "util.h"

const double EPS_CONVERGE = 1e-5;

// [[Rcpp::export]]
Rcpp::List seq_epmgp(arma::vec m, arma::mat Kinv, arma::vec lb, arma::vec ub,
                     int max_steps) {
  
  arma::vec nu_site = arma::zeros(Kinv.n_rows);
  arma::vec tau_site = arma::zeros(Kinv.n_rows);
  
  // initialize q(x)
  arma::mat Lambda = Kinv;
  // arma::vec eta = Lambda * (lb + ub) / 2;
  arma::vec eta = arma::zeros(Lambda.n_rows);
  
  Rcpp::Rcout << "Lambda: " << Lambda << std::endl;
  Rcpp::Rcout << "eta: " << eta << std::endl;
  
  double logZ;
  arma::vec eta_last = -arma::datum::inf * arma::ones(arma::size(eta));
  bool converged = false;
  int k = 1;
  int n = lb.n_elem;
  
  // algorithm loop
  arma::vec tau_cavity;
  arma::vec nu_cavity;
  arma::mat L;
  arma::vec logz_hat;
  arma::vec sigma_hat;
  arma::vec mu_hat;
  
  while (!converged) {
    
    // make cavity distribution
    tau_cavity = arma::diagvec(Lambda) - tau_site;
    nu_cavity = eta - nu_site; 
    
    // compute moments using truncated normals
    Rcpp::List moments = trunc_norm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity);
    arma::vec sigma_hat_out = moments["sigma_hat"];
    arma::vec logz_hat_out = moments["logz_hat"];
    arma::vec mu_hat_out = moments["mu_hat"];
    logz_hat = logz_hat_out;
    sigma_hat = sigma_hat_out;
    mu_hat = mu_hat_out;
    
    Rcpp::Rcout << "sigma hat: " << sigma_hat << std::endl;
    Rcpp::Rcout << "tau cavity: " << tau_cavity << std::endl;
    Rcpp::Rcout << "nu cavity: " << nu_cavity << std::endl;
    
    // update the site parameters
    arma::vec delta_tau_site = (1 / sigma_hat) - tau_cavity - tau_site;
    arma::vec delta_nu_site = (mu_hat / sigma_hat) - nu_cavity - nu_site; 
    tau_site += delta_tau_site;
    // nu_site = (mu_hat / sigma_hat) - nu_cavity;
    nu_site += delta_nu_site;
    
    // enforce nonnegativity of tau_site
    if (arma::any(tau_site < 0)) {
      for (int i = 0; i < tau_site.n_elem; i++) {
        if (tau_site(i) > -1e-8) {
          tau_site(i) = 0.0;
        }
      }
    }
    
    Rcpp::Rcout << "delta tau site: " << delta_tau_site << std::endl;
    Rcpp::Rcout << "nu site: " << nu_site << std::endl;
    
    // update q(x) Lambda and eta
    Lambda.diag(0) += delta_tau_site;
    eta += delta_nu_site;
    
    // check convergence criteria
    if ((arma::norm(eta_last - eta)) < EPS_CONVERGE || k == max_steps)
      converged = true;
    else
      eta_last = eta;
    k++;
  }
  
  arma::vec mu_cavity = nu_cavity / tau_cavity;
  double lZ1 = .5*arma::sum(arma::log(tau_site / tau_cavity + 1));
  double lZ2 = .5*arma::sum(
    (arma::square(mu_cavity) % tau_site -
      2 * mu_cavity % nu_site -
      arma::square(nu_site) / tau_cavity) /
      (1 + tau_site / tau_cavity)
  );
  // double lZ2 = 0.5 * arma::as_scalar(
  //   nu_cavity.t() * 
  //     arma::solve(arma::diagmat(tau_site) + arma::diagmat(tau_cavity), 
  //                 tau_site % nu_cavity / tau_cavity - 2 * nu_site)
  // );
                   
  // double val;
  // double sign;
  // arma::log_det(val, sign, Lambda);
  double lZ3 = .5*arma::as_scalar(eta.t() * Lambda * eta);
  logZ = arma::sum(logz_hat) + lZ1 + lZ2 + lZ3;
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["logZ"] = logZ,
    Rcpp::_["eta"] = eta,
    Rcpp::_["Lambda"] = Lambda
  );
  
  return result;
}