// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// axisepmgp
Rcpp::List axisepmgp(arma::vec m, arma::mat K, arma::vec lb, arma::vec ub);
RcppExport SEXP _epmgp_axisepmgp(SEXP mSEXP, SEXP KSEXP, SEXP lbSEXP, SEXP ubSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ub(ubSEXP);
    rcpp_result_gen = Rcpp::wrap(axisepmgp(m, K, lb, ub));
    return rcpp_result_gen;
END_RCPP
}
// epmgp
Rcpp::List epmgp(arma::vec m, arma::mat K, arma::mat C, arma::vec lb, arma::vec ub, int max_steps);
RcppExport SEXP _epmgp_epmgp(SEXP mSEXP, SEXP KSEXP, SEXP CSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(epmgp(m, K, C, lb, ub, max_steps));
    return rcpp_result_gen;
END_RCPP
}
// sample_epmh
arma::mat sample_epmh(int n_samples, arma::vec ep_mean, arma::mat ep_chol, arma::mat F, arma::vec g, arma::vec initial);
RcppExport SEXP _epmgp_sample_epmh(SEXP n_samplesSEXP, SEXP ep_meanSEXP, SEXP ep_cholSEXP, SEXP FSEXP, SEXP gSEXP, SEXP initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ep_mean(ep_meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ep_chol(ep_cholSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type g(gSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initial(initialSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_epmh(n_samples, ep_mean, ep_chol, F, g, initial));
    return rcpp_result_gen;
END_RCPP
}
// range_intersection
arma::vec range_intersection(arma::vec first, arma::vec second);
RcppExport SEXP _epmgp_range_intersection(SEXP firstSEXP, SEXP secondSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type first(firstSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type second(secondSEXP);
    rcpp_result_gen = Rcpp::wrap(range_intersection(first, second));
    return rcpp_result_gen;
END_RCPP
}
// sample_epess
arma::mat sample_epess(int n_samples, arma::vec ep_mean, arma::mat ep_chol, arma::mat F, arma::vec g, int J, int N, arma::vec initial);
RcppExport SEXP _epmgp_sample_epess(SEXP n_samplesSEXP, SEXP ep_meanSEXP, SEXP ep_cholSEXP, SEXP FSEXP, SEXP gSEXP, SEXP JSEXP, SEXP NSEXP, SEXP initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ep_mean(ep_meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ep_chol(ep_cholSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initial(initialSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_epess(n_samples, ep_mean, ep_chol, F, g, J, N, initial));
    return rcpp_result_gen;
END_RCPP
}
// test_wall_hit
arma::vec test_wall_hit(arma::vec nu, arma::vec initial, arma::vec ep_mean, arma::mat ep_chol, arma::mat F, arma::vec g, int J, int N);
RcppExport SEXP _epmgp_test_wall_hit(SEXP nuSEXP, SEXP initialSEXP, SEXP ep_meanSEXP, SEXP ep_cholSEXP, SEXP FSEXP, SEXP gSEXP, SEXP JSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ep_mean(ep_meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ep_chol(ep_cholSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(test_wall_hit(nu, initial, ep_mean, ep_chol, F, g, J, N));
    return rcpp_result_gen;
END_RCPP
}
// erfcx
double erfcx(double x);
RcppExport SEXP _epmgp_erfcx(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(erfcx(x));
    return rcpp_result_gen;
END_RCPP
}
// trunc_norm_moments
Rcpp::List trunc_norm_moments(arma::vec lb_in, arma::vec ub_in, arma::vec mu_in, arma::vec sigma_in);
RcppExport SEXP _epmgp_trunc_norm_moments(SEXP lb_inSEXP, SEXP ub_inSEXP, SEXP mu_inSEXP, SEXP sigma_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lb_in(lb_inSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ub_in(ub_inSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_in(mu_inSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_in(sigma_inSEXP);
    rcpp_result_gen = Rcpp::wrap(trunc_norm_moments(lb_in, ub_in, mu_in, sigma_in));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epmgp_axisepmgp", (DL_FUNC) &_epmgp_axisepmgp, 4},
    {"_epmgp_epmgp", (DL_FUNC) &_epmgp_epmgp, 6},
    {"_epmgp_sample_epmh", (DL_FUNC) &_epmgp_sample_epmh, 6},
    {"_epmgp_range_intersection", (DL_FUNC) &_epmgp_range_intersection, 2},
    {"_epmgp_sample_epess", (DL_FUNC) &_epmgp_sample_epess, 8},
    {"_epmgp_test_wall_hit", (DL_FUNC) &_epmgp_test_wall_hit, 8},
    {"_epmgp_erfcx", (DL_FUNC) &_epmgp_erfcx, 1},
    {"_epmgp_trunc_norm_moments", (DL_FUNC) &_epmgp_trunc_norm_moments, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_epmgp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
