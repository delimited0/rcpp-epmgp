#include "util.h"

const double EPS_CONVERGE = 1e-5;

// [[Rcpp::export]]
Rcpp::List seq_moment_epmgp(arma::vec m, arma::mat K, 
                            arma::vec lb, arma::vec ub,
                            int max_steps) {
  
  arma::vec nu_site = arma::zeros(K.n_rows);
  arma::vec tau_site = arma::zeros(K.n_rows);
  
  // initialize q(x)
  arma::mat Sigma = K;
  arma::vec mu = (lb + ub) / 2;
  for (int i = 0; i < mu.n_elem; i++) {
    if (std::isinf(mu(i))) {
      mu(i) = copysign(1.0, mu(i)) * 10;
    }
  }
  
  double logZ = arma::datum::inf;
  arma::vec mu_last = -arma::datum::inf * arma::ones(mu.n_elem);
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
    
    Rcpp::Rcout << "iteration " << k << " ==========" << std::endl;
    
    // make cavity distribution
    tau_cavity = arma::diagvec(Lambda) - tau_site;
    nu_cavity = eta - nu_site; 
    
    Rcpp::Rcout << "tau cavity: " << tau_cavity << std::endl;
    Rcpp::Rcout << "nu cavity: " << nu_cavity << std::endl;
    
    // compute moments using truncated normals
    Rcpp::List moments = trunc_norm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity);
    arma::vec sigma_hat_out = moments["sigma_hat"];
    arma::vec logz_hat_out = moments["logz_hat"];
    arma::vec mu_hat_out = moments["mu_hat"];
    logz_hat = logz_hat_out;
    sigma_hat = sigma_hat_out;
    mu_hat = mu_hat_out;
    
    Rcpp::Rcout << "mu hat: " << mu_hat << std::endl;
    Rcpp::Rcout << "sigma hat: " << sigma_hat << std::endl;
    
    // update the site parameters
    arma::vec delta_tau_site = (1 / sigma_hat - tau_cavity) - tau_site;
    arma::vec delta_nu_site = ((mu_hat / sigma_hat) - nu_cavity) - nu_site; 
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
    Rcpp::Rcout << "delta nu site: " << delta_nu_site << std::endl;
    Rcpp::Rcout << "tau site: " << tau_site << std::endl;
    Rcpp::Rcout << "nu site: " << nu_site << std::endl;
    
    // update q(x) Lambda and eta
    Lambda.diag(0) += delta_tau_site;
    eta += delta_nu_site;
    // Lambda = Kinv + arma::diagmat(tau_site);
    // eta = Kinv * m + nu_site;
    
    Rcpp::Rcout << "Lambda: " << Lambda << std::endl;
    Rcpp::Rcout << "eta:" << eta << std::endl;
    
    // check convergence criteria
    if ((arma::norm(eta_last - eta)) < EPS_CONVERGE || k == max_steps)
      converged = true;
    else
      eta_last = eta;
    k++;
  }
  
  // arma::vec mu_cavity = nu_cavity / tau_cavity;
  // double lZ1 = .5*arma::sum(arma::log(tau_site / tau_cavity + 1));
  // double lZ2 = .5*arma::sum(
    //   (arma::square(mu_cavity) % tau_site -
            //     2 * mu_cavity % nu_site -
            //     arma::square(nu_site) / tau_cavity) /
      //     (1 + tau_site / tau_cavity)
    // );
  // double lZ2 = 0.5 * arma::as_scalar(
    //   nu_cavity.t() * 
      //     arma::solve(arma::diagmat(tau_site) + arma::diagmat(tau_cavity), 
                         //                 tau_site % nu_cavity / tau_cavity - 2 * nu_site)
    // );
  
  // double val;
  // double sign;
  // arma::log_det(val, sign, Lambda);
  // double lZ3 = .5*arma::as_scalar(eta.t() * Lambda * eta);
  // logZ = arma::sum(logz_hat) + lZ1 + lZ2 + lZ3;
  
  Rcpp::List result = Rcpp::List::create(
    // Rcpp::_["logZ"] = logZ,
    Rcpp::_["tau_site"] = tau_site,
    Rcpp::_["tau_cavity"] = tau_cavity,
    Rcpp::_["nu_cavity"] = nu_cavity,
    Rcpp::_["eta"] = eta,
    Rcpp::_["Lambda"] = Lambda
  );
  
  return result;
}