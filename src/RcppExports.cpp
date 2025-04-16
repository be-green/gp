// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// linterp
arma::mat linterp(arma::mat& x, const arma::vec& z, const arma::vec& c);
RcppExport SEXP _gp_linterp(SEXP xSEXP, SEXP zSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(linterp(x, z, c));
    return rcpp_result_gen;
END_RCPP
}
// gh_quadrature
void gh_quadrature(double mu_n, double sigma_n, double mu_u, double sigma_u, arma::mat& nodes2D, arma::vec& weights2D);
RcppExport SEXP _gp_gh_quadrature(SEXP mu_nSEXP, SEXP sigma_nSEXP, SEXP mu_uSEXP, SEXP sigma_uSEXP, SEXP nodes2DSEXP, SEXP weights2DSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu_n(mu_nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< double >::type mu_u(mu_uSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type nodes2D(nodes2DSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type weights2D(weights2DSEXP);
    gh_quadrature(mu_n, sigma_n, mu_u, sigma_u, nodes2D, weights2D);
    return R_NilValue;
END_RCPP
}
// net_euler_diff
arma::vec net_euler_diff(arma::vec c_now, arma::vec x, arma::vec c_next, arma::mat const_scale_coh, arma::mat const_add_coh, arma::mat const_scale_consump, arma::vec weights, double R, double p_noinc, double beta, double rho);
RcppExport SEXP _gp_net_euler_diff(SEXP c_nowSEXP, SEXP xSEXP, SEXP c_nextSEXP, SEXP const_scale_cohSEXP, SEXP const_add_cohSEXP, SEXP const_scale_consumpSEXP, SEXP weightsSEXP, SEXP RSEXP, SEXP p_noincSEXP, SEXP betaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c_now(c_nowSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c_next(c_nextSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type const_scale_coh(const_scale_cohSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type const_add_coh(const_add_cohSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type const_scale_consump(const_scale_consumpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type p_noinc(p_noincSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(net_euler_diff(c_now, x, c_next, const_scale_coh, const_add_coh, const_scale_consump, weights, R, p_noinc, beta, rho));
    return rcpp_result_gen;
END_RCPP
}
// nw
Rcpp::List nw(double sigma_n, double sigma_u);
RcppExport SEXP _gp_nw(SEXP sigma_nSEXP, SEXP sigma_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    rcpp_result_gen = Rcpp::wrap(nw(sigma_n, sigma_u));
    return rcpp_result_gen;
END_RCPP
}
// consumption_rule
arma::mat consumption_rule(arma::vec& x, arma::vec& G, double sigma_n, double sigma_u, double gamma_0, double gamma_1, double R, double p_noinc, double beta, double rho);
RcppExport SEXP _gp_consumption_rule(SEXP xSEXP, SEXP GSEXP, SEXP sigma_nSEXP, SEXP sigma_uSEXP, SEXP gamma_0SEXP, SEXP gamma_1SEXP, SEXP RSEXP, SEXP p_noincSEXP, SEXP betaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< double >::type gamma_1(gamma_1SEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type p_noinc(p_noincSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(consumption_rule(x, G, sigma_n, sigma_u, gamma_0, gamma_1, R, p_noinc, beta, rho));
    return rcpp_result_gen;
END_RCPP
}
// simulate_assets
arma::vec simulate_assets(int N, double mu, double sigma);
RcppExport SEXP _gp_simulate_assets(SEXP NSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_assets(N, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// draw_bernoulli
arma::mat draw_bernoulli(double p, int N, int T);
RcppExport SEXP _gp_draw_bernoulli(SEXP pSEXP, SEXP NSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_bernoulli(p, N, T));
    return rcpp_result_gen;
END_RCPP
}
// simulate_income
Rcpp::List simulate_income(int N, int T, arma::vec& P_init, arma::vec& G, double sigma_n, double sigma_u, double p_noinc);
RcppExport SEXP _gp_simulate_income(SEXP NSEXP, SEXP TSEXP, SEXP P_initSEXP, SEXP GSEXP, SEXP sigma_nSEXP, SEXP sigma_uSEXP, SEXP p_noincSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type P_init(P_initSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    Rcpp::traits::input_parameter< double >::type p_noinc(p_noincSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_income(N, T, P_init, G, sigma_n, sigma_u, p_noinc));
    return rcpp_result_gen;
END_RCPP
}
// consume
arma::vec consume(arma::vec x, arma::vec& x_grid, arma::mat& cr, int t);
RcppExport SEXP _gp_consume(SEXP xSEXP, SEXP x_gridSEXP, SEXP crSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type cr(crSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(consume(x, x_grid, cr, t));
    return rcpp_result_gen;
END_RCPP
}
// simulate_lifecycle
Rcpp::List simulate_lifecycle(int N, int T, arma::vec& x_grid, arma::mat& N_shock, arma::mat& U_shock, arma::mat& P, arma::vec& init_a, arma::vec& G, double sigma_n, double sigma_u, double gamma_0, double gamma_1, double R, double p_noinc, double beta, double rho);
RcppExport SEXP _gp_simulate_lifecycle(SEXP NSEXP, SEXP TSEXP, SEXP x_gridSEXP, SEXP N_shockSEXP, SEXP U_shockSEXP, SEXP PSEXP, SEXP init_aSEXP, SEXP GSEXP, SEXP sigma_nSEXP, SEXP sigma_uSEXP, SEXP gamma_0SEXP, SEXP gamma_1SEXP, SEXP RSEXP, SEXP p_noincSEXP, SEXP betaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type N_shock(N_shockSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type U_shock(U_shockSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type init_a(init_aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< double >::type gamma_1(gamma_1SEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type p_noinc(p_noincSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_lifecycle(N, T, x_grid, N_shock, U_shock, P, init_a, G, sigma_n, sigma_u, gamma_0, gamma_1, R, p_noinc, beta, rho));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gp_linterp", (DL_FUNC) &_gp_linterp, 3},
    {"_gp_gh_quadrature", (DL_FUNC) &_gp_gh_quadrature, 6},
    {"_gp_net_euler_diff", (DL_FUNC) &_gp_net_euler_diff, 11},
    {"_gp_nw", (DL_FUNC) &_gp_nw, 2},
    {"_gp_consumption_rule", (DL_FUNC) &_gp_consumption_rule, 10},
    {"_gp_simulate_assets", (DL_FUNC) &_gp_simulate_assets, 3},
    {"_gp_draw_bernoulli", (DL_FUNC) &_gp_draw_bernoulli, 3},
    {"_gp_simulate_income", (DL_FUNC) &_gp_simulate_income, 7},
    {"_gp_consume", (DL_FUNC) &_gp_consume, 4},
    {"_gp_simulate_lifecycle", (DL_FUNC) &_gp_simulate_lifecycle, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_gp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
