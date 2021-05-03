#include "util.h"

const double EPS_CONVERGE = 1e-5;

// [[Rcpp::export]]
Rcpp::List aug_epmgp(arma::vec m, arma::mat K, arma::vec lb, arma::vec ub,
                     int max_steps, int p) {
  int cp = m.n_elem;
  // int p = lb.n_elem;  // probability dimension -- but not necessarily
  int c = cp / p;     // number of constraints
  
  arma::vec nu_site = arma::zeros(c);
  arma::vec tau_site = arma::zeros(c);
  
  arma::vec nu_cavity(c);
  arma::vec tau_cavity(c);
  arma::vec delta_tau_site(c);
  arma::vec delta_nu_site(c);
  arma::vec logz_hat(c);
  arma::vec mu_hat(c);
  arma::vec sigma_hat(c);
  
  arma::mat Sigma = K;
  arma::mat mu = m;
  
  double logz = arma::datum::inf;
  arma::vec mu_last = -arma::datum::inf * arma::ones(mu.n_elem);
  bool converged = false;
  int k = 0;
  
  while (!converged && k < max_steps) {
    for (int j = 0; j < c; j++) {
      
      // Rcpp::Rcout << "Factor " << j << std::endl;
      int jprime = c*(p - 1) + j;
      // Rcpp::Rcout << "jprime : " << jprime << std::endl;
      
      // make the cavity distribution
      tau_cavity(j) = 1 / Sigma(jprime, jprime) - tau_site(j);
      nu_cavity(j) = mu(jprime) / Sigma(jprime, jprime) - nu_site(j);
 
      // Rcpp::Rcout << "Made cavity distribution" << std::endl;
      // Rcpp::Rcout << "tau_cavity: " << tau_cavity << std::endl;
      // Rcpp::Rcout << "nu_cavity: " << nu_cavity << std ::endl;
      
      // compute moments using truncated normals
      Rcpp::List moments = trunc_norm_moments(
        lb.row(j), ub.row(j), 
        nu_cavity.row(j) / tau_cavity.row(j), 1. / tau_cavity.row(j));
      
      arma::vec sigma_hat_out = moments["sigma_hat"];
      arma::vec logz_hat_out = moments["logz_hat"];
      arma::vec mu_hat_out = moments["mu_hat"];
      
      sigma_hat(j) = sigma_hat_out(0);
      logz_hat(j) = logz_hat_out(0);
      mu_hat(j) = mu_hat_out(0);
      
      // Rcpp::Rcout << "Made moments" << std::endl;
      // Rcpp::Rcout << "sigma_hat: " << sigma_hat << std::endl;
      // Rcpp::Rcout << "logz_hat: " << logz_hat << std::endl;
      // Rcpp::Rcout << "mu_hat: " << mu_hat << std::endl;
      
      if (sigma_hat(j) == 0) {
        // the algorithm has found a 0 weight dimension, terminate
        converged = true;
        logz = -arma::datum::inf;
        mu = arma::datum::nan;
        Sigma = arma::datum::nan;
        break;
      }
      
      delta_tau_site(j) = (1 / sigma_hat(j) - tau_cavity(j)) - tau_site(j);
      delta_nu_site(j) = (mu_hat(j) / sigma_hat(j) - nu_cavity(j)) - nu_site(j);
      tau_site(j) += delta_tau_site(j);
      nu_site(j) += delta_nu_site(j);
      
      // Rcpp::Rcout << "Updated sites" << std::endl;
      // Rcpp::Rcout << "delta_tau_site: " << delta_tau_site << std::endl;
      // Rcpp::Rcout << "delta_nu_site: " << delta_nu_site << std::endl;
      // Rcpp::Rcout << "tau_site: " << tau_site << std::endl;
      // Rcpp::Rcout << "nu_site: " << nu_site << std::endl;
      
      if (tau_site(j) < 0) {
        // if result negative, either due to numerical precision or error
        if (tau_site(j) > -1e-6)
          tau_site(j) = 0;
      }
      
      // update q(x) (Sigma and mu)
      arma::vec sc = Sigma.col(jprime);
      Sigma -= 
        (delta_tau_site(j) / (1 + delta_tau_site(j) * sc(jprime))) *
        (sc * sc.t());
      mu += 
        ((delta_nu_site(j) - delta_tau_site(j) * mu(jprime))) / 
        (1 + delta_tau_site(j) * sc(jprime)) * sc;
      
      // Rcpp::Rcout << "sc: " << sc << std::endl;
      // Rcpp::Rcout << "Sigma: " << Sigma << std::endl;
      // Rcpp::Rcout << "mu: " << mu << std::endl;
    }
    
    // check convergence criteria
    if (arma::norm(mu_last - mu, 2) < EPS_CONVERGE)
      converged = true;
    else
      mu_last = mu;
    
    k++;
  }
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["mu"] = mu,
    Rcpp::_["Sigma"] = Sigma,
    Rcpp::_["logz_hat"] = logz_hat,
    Rcpp::_["nu_cavity"] = nu_cavity,
    Rcpp::_["tau_cavity"] = tau_cavity,
    Rcpp::_["nu_site"] = nu_site,
    Rcpp::_["tau_site"] = tau_site
  );
  
  return result;
}