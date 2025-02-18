#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Finds the roots of a function that takes a double vector as input
//' @param f function to find the roots of
//' @param intervals matrix of intervals to search over
//' @param tol tolerance for what is equivalent to zero
//' @param max_iter max number of iterations
//' @details Binary search (bisection) to find a root of f in the interval [a, b]
//' for a vector
arma::vec find_roots(std::function<arma::vec(arma::vec)>& f,
                     arma::vec a,
                     arma::vec b,
                     double tol = 1e-2,
                     int maxIter = 1000) {

  arma::vec fa = f(a);
  arma::vec fb = f(b);

  arma::vec mid(a.n_rows, arma::fill::zeros);
  arma::uvec a_idx;
  arma::uvec b_idx;
  arma::vec fmid;
  for (int iter = 0; iter < maxIter; ++iter) {
    mid = (a + b) / 2.0;
    fmid = f(mid);

    // Check if the mid value is close enough to a root.
    if (arma::sum(arma::abs(fmid)) < tol || arma::sum((b - a) / 2.0) < tol) {
      return mid;
    }

    // Decide which half-interval to keep.
    b_idx = arma::find(fa % fmid < 0);
    b(b_idx) = mid(b_idx);
    fb(b_idx) = fmid(b_idx);

    a_idx = arma::find(fa % fmid >= 0);
    a(a_idx) = mid(a_idx);
    fa(a_idx) = fmid(a_idx);

  }

  // If we reach here, max iterations were exceeded.
  Rcpp::Rcerr << "Warning: Maximum iterations reached. Returning best approximation." << std::endl;
  return mid;
}

//' Linear interpolation
//' @param z input parameter to linearly interpolate the x coordinates
//' @param c output vector to weight based on the input
//' @details Returns a matrix of linearly interpolated values
// [[Rcpp::export]]
arma::mat linterp(arma::mat x, const arma::vec& z, const arma::vec& c) {
 // Ensure both vectors have the same size and contain at least two points.
 if (z.n_elem != c.n_elem) {
   throw std::invalid_argument("Vectors z and c must have the same size.");
 }
 if (z.n_elem < 2) {
   throw std::invalid_argument("At least two points are required for interpolation.");
 }

 arma::mat y(x.n_rows, x.n_cols);  // Output vector of interpolated values

 for (arma::uword i = 0; i < x.n_rows; ++i) {
   for (arma::uword j = 0; j < x.n_cols; ++j) {
     double xi = x(i, j);

     // Extrapolate to the left if xi is less than the first element of z.
     if (xi <= z(0)) {
       double t = (xi - z(0)) / (z(1) - z(0));
       y(i, j) = c(0) + t * (c(1) - c(0));
       continue;
     }

     // Extrapolate to the right if xi is greater than the last element of z.
     if (xi >= z(z.n_elem - 1)) {
       double t = (xi - z(z.n_elem - 2)) / (z(z.n_elem - 1) - z(z.n_elem - 2));
       y(i, j) = c(z.n_elem - 2) + t * (c(z.n_elem - 1) - c(z.n_elem - 2));
       continue;
     }

     // Binary search: Find indices 'low' and 'high' such that z(low) <= xi <= z(high)
     arma::uword low = 0;
     arma::uword high = z.n_elem - 1;
     while (high - low > 1) {
       arma::uword mid = (low + high) / 2;
       if (z(mid) <= xi)
         low = mid;
       else
         high = mid;
     }

     // Linear interpolation between the points (z(low), c(low)) and (z(high), c(high))
     double t = (xi - z(low)) / (z(high) - z(low));
     y(i, j) = c(low) + t * (c(high) - c(low));
   }
 }
 return y;
}

 //' Computes the 2D Gaussian quadrature rule of order 12
 //' @param mu_n mean of permanent shocks
 //' @param sigma sd of permanent shocks
 //' @param mu_u mean of temp shocks
 //' @param sigma_u sd of temp shocks
 void gh_quadrature(double mu_n, double sigma_n,
                    double mu_u, double sigma_u,
                    arma::mat& nodes2D,
                    arma::vec& weights2D) {

   // nodes and weights for quadrature of order 12
   // see, e.g. https://www.efunda.com/math/num_integration/findgausshermite.cfm
   // or https://search.r-project.org/CRAN/refmans/fastGHQuad/html/ghQuad.html
   arma::vec nodes;
   arma::vec weights;

   weights = {2.65855168435631e-07, 8.57368704358781e-05,
              0.00390539058462906, 0.051607985615884,
              0.260492310264161, 0.57013523626248,
              0.57013523626248, 0.260492310264161,
              0.0516079856158839, 0.00390539058462906,
              8.57368704358789e-05, 2.65855168435631e-07};

   nodes = {-3.88972489786978, -3.02063702512089,
              -2.27950708050106, -1.5976826351526,
              -0.947788391240164, -0.314240376254359,
              0.314240376254359, 0.947788391240164,
              1.5976826351526, 2.27950708050106,
              3.02063702512089, 3.88972489786978};
   // Transform the 1D nodes and weights to integrate against a normal density.
   // For dimension i: x = sqrt(2)*sigma*node + mu, and weight = weight / sqrt(pi)
   arma::vec x_n = std::sqrt(2.0) * sigma_n * nodes + mu_n;
   arma::vec x_u = std::sqrt(2.0) * sigma_u * nodes + mu_u;

   // puts the normalization on the weights instead of outside the sum
   arma::vec w_n = weights / std::sqrt(arma::datum::pi);
   arma::vec w_u = weights / std::sqrt(arma::datum::pi);

   // Form the 2D quadrature rule by taking the Cartesian product of the 1D rules.
   int n = 12;
   int idx = 0;
   for (int i = 0; i < n; ++i) {
     for (int j = 0; j < n; ++j) {
       nodes2D(idx, 0) = x_n(i);
       nodes2D(idx, 1) = x_u(j);
       weights2D(idx) = w_n(i) * w_u(j);
       idx++;
     }
   }
}

//' Marginal utility
//' @param c consumption
//' @param rho risk aversion
arma::mat u_prime(arma::mat c, double rho) {
  arma::mat mu = pow(c, -rho);
  return mu;
}

//' Computes difference between expected marginal utility and
//' current marginal utility
// [[Rcpp::export]]
arma::vec net_euler_diff(arma::vec c_now,
                         arma::vec x,
                         arma::vec c_next,
                         arma::mat const_scale_coh,
                         arma::mat const_add_coh,
                         arma::mat const_scale_consump,
                         arma::vec weights,
                         double R,
                         double p_noinc,
                         double beta,
                         double rho) {

  arma::mat x_minus_c_now = repmat(x - c_now, 1, const_scale_coh.n_cols);
  arma::mat noincome_consumption_matrix;
  arma::mat income_consumption_matrix;

  // elementwise multiplication in arma is %
  // input to consumption function
  arma::mat noincome_next_coh = x_minus_c_now % const_scale_coh;

  noincome_consumption_matrix =
   linterp(noincome_next_coh, c_next, x) %
   const_scale_consump;

  arma::mat income_next_coh = x_minus_c_now % (const_scale_coh + const_add_coh);
  income_consumption_matrix =
   linterp(income_next_coh, c_next, x) %
   const_scale_consump;

  // marginal utility for each state / consumption pair
  arma::mat noincome_mu = u_prime(noincome_consumption_matrix, rho);
  arma::mat income_mu = u_prime(income_consumption_matrix, rho);

  // let's integrate into a column vector
  arma::vec e_mu_tp1 = p_noinc * (noincome_mu * weights) +
   (1 - p_noinc) * (income_mu * weights);
  arma::vec mu_t = u_prime(c_now, rho);
  arma::vec euler_diff = mu_t - beta * R * e_mu_tp1;
  return euler_diff;
}

// [[Rcpp::export]]
//' Find roots of Euler Equation
//' @param c_now this period's consumption
//' @param c_tp1 next period's consumption
//' @param x cash on hand
//' @param R risk free rate
//' @param G growth rate of income
//' @param sigma_n sd of log permanent income
//' @param sigma_u sd of log temp. income
//' @param beta time discount factor
arma::vec solve_euler(arma::vec c_next_grid,
                      arma::vec x,
                      double R,
                      double G,
                      double sigma_n,
                      double sigma_u,
                      double p_noinc,
                      double beta,
                      double rho) {

  // both shocks are mean zero
  // so change of variables has mu 0

  // assigned by ref via gh_quadrature
  arma::mat nodes(144, 2, arma::fill::zeros);
  arma::vec weights(144, arma::fill::zeros);
  gh_quadrature(
    0, sigma_n,
    0, sigma_u,
    nodes, weights
  );

  // these are ugly expressions from Gourinchas-Parker appendix A.1
  arma::mat const_scale_coh =
    repmat((R / G) * exp(sqrt(2) * sigma_n * nodes.col(0).t()), x.n_elem, 1);
  arma::mat const_add_coh =
    repmat(exp(-sqrt(2) * sigma_u * nodes.col(1).t()), x.n_elem, 1);
  arma::mat const_scale_consump =
    repmat(G * exp(-sqrt(2) * sigma_n * nodes.col(0).t()), x.n_elem, 1);

  std::function<arma::vec(arma::vec)> f = [=](arma::vec c) -> arma::vec {
    arma::vec diff = net_euler_diff(c, x, c_next_grid,
                                    const_scale_coh,
                                    const_add_coh,
                                    const_scale_consump,
                                    weights,
                                    R, p_noinc,
                                    beta, rho);
    return diff;
  };

  // no negative consumption, can't consume more than available cash
  arma::vec low(x.n_rows, arma::fill::zeros);
  arma::vec c_now = find_roots(f, low, x);
  return c_now;
}

// [[Rcpp::export]]
//' @param x cash on hand grid
//' @param sigma_n standard deviation of log permanent income
//' @param sigma_u standard deviation of log transitory income
//' @param gamma_0 intercept of retirement MPC
//' @param gamma_1 slope of retirement MPC in cash on hand
//' @param R risk free rate
//' @param G permanent income growth rate
//' @param T number of years
//' @param beta time discount factor
//' @export
arma::mat consumption_path(arma::vec x,
                           arma::vec G,
                           double sigma_n,
                           double sigma_u,
                           double gamma_0,
                           double gamma_1,
                           double R,
                           double p_noinc,
                           double beta,
                           double rho) {
 int N = x.n_elem;
 int T = G.n_elem;

 // consumption vector
 arma::mat c(N, T, arma::fill::zeros);

 // initialize at linear retirement rule
 c.col(T - 1) = gamma_0 + gamma_1 * x;

 // solve backwards by iterating on the euler
 for(int t = T - 1; t > 0; --t) {
   c.col(t - 1) =
     solve_euler(c.col(t),
                 x,
                 R,
                 G(t),
                 sigma_n,
                 sigma_u,
                 p_noinc,
                 beta,
                 rho);
 }
 return c;
}
