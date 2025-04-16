#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Finds the roots of a function that takes a double vector as input
//' @param f function to find the roots of
//' @param intervals matrix of intervals to search over
//' @param tol tolerance for what is equivalent to zero
//' @param max_iter max number of iterations
//' @details Binary search (bisection) to find a
//' root of f in the interval [a, b]
//' for a vector
arma::vec find_roots(std::function<arma::vec(arma::vec)> f,
                     arma::vec a,
                     arma::vec b,
                     double tol = 1e-5,
                     int maxIter = 1000) {

   arma::vec fa = f(a);
   arma::vec fb = f(b);
   arma::vec mid(a.n_elem, arma::fill::zeros);
   arma::vec fmid;

   for (int iter = 0; iter < maxIter; ++iter) {
     mid = (a + b) / 2.0;
     fmid = f(mid);

     // Check for convergence in each interval using the maximum value.
     if (arma::max(arma::abs(fmid)) < tol || arma::max((b - a) / 2.0) < tol) {
       return mid;
     }

     // Update intervals elementwise based on the sign change.
     arma::uvec b_idx = arma::find(fa % fmid < 0);
     arma::uvec a_idx = arma::find(fa % fmid >= 0);

     b.elem(b_idx) = mid.elem(b_idx);
     fb.elem(b_idx) = fmid.elem(b_idx);

     a.elem(a_idx) = mid.elem(a_idx);
     fa.elem(a_idx) = fmid.elem(a_idx);
   }

   Rcpp::Rcerr << "Warning: Maximum iterations reached. Returning best approximation." << std::endl;
   return mid;
}

//' Linear interpolation
//' @param z input parameter to linearly interpolate the x coordinates
//' @param c output vector to weight based on the input
//' @details Returns a matrix of linearly interpolated values
//' @export
// [[Rcpp::export]]
arma::mat linterp(arma::mat& x, const arma::vec& z, const arma::vec& c) {
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
//' @export
// [[Rcpp::export]]
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
//' @param c_now current consumption
//' @param x current cash on hand
//' @param c_next proposed next consumption
//' @param const_scale_coh cash on hand scale param (see paper)
//' @param const_add_coh cash on hand addition param (see paper)
//' @param const_scale_consump consumption scale param (see paper)
//' @param weights Weights for weighted average from Gauss-Hermite
//' quadrature
//' @param R risk free rate
//' @param p_noinc probability of income equal to 0
//' @param beta discount factor
//' @param rho CRRA risk aversion parameter
//' @export
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
 arma::mat c_next_mat = repmat(c_next, 1, const_scale_coh.n_cols);
 arma::mat noincome_consumption_matrix;
 arma::mat income_consumption_matrix;

 // when the u shock is equal to 0, u = 0
 // which makes exp(-\sqrt(2)\sigma_u u) = 1
 noincome_consumption_matrix =
   c_next_mat % x_minus_c_now % (const_scale_coh + 1) %
   const_scale_consump;

 income_consumption_matrix =
   c_next_mat % x_minus_c_now % (const_scale_coh + const_add_coh) %
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

//' Function just for testing quadrature
// [[Rcpp::export]]
Rcpp::List nw(double sigma_n, double sigma_u) {
  // assigned by ref via gh_quadrature
  arma::mat nodes(144, 2, arma::fill::zeros);
  arma::vec weights(144, arma::fill::zeros);
  gh_quadrature(
    0, sigma_n,
    0, sigma_u,
    nodes, weights
  );
  return Rcpp::List::create(
    Rcpp::Named("nodes") = nodes,
    Rcpp::Named("weights") = weights
  );
}

//' Find roots of Euler Equation
//' @param c_now this period's consumption
//' @param c_tp1 next period's consumption
//' @param x cash on hand
//' @param R risk free rate
//' @param G growth rate of income
//' @param sigma_n sd of log permanent income
//' @param sigma_u sd of log temp. income
//' @param beta time discount factor
//' @export
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
   repmat((R / G) * exp(sqrt(2) * sigma_n *
          nodes.col(0).t()), x.n_elem, 1);
 arma::mat const_add_coh =
   repmat(exp(-sqrt(2) * sigma_u *
              nodes.col(1).t()), x.n_elem, 1);
 arma::mat const_scale_consump =
   repmat(G * exp(-sqrt(2) * sigma_n *
          nodes.col(0).t()), x.n_elem, 1);

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

//' Get full consumption rule
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
// [[Rcpp::export]]
arma::mat consumption_rule(arma::vec& x,
                           arma::vec& G,
                           double sigma_n,
                           double sigma_u,
                           double gamma_0,
                           double gamma_1,
                           double R,
                           double p_noinc,
                           double beta,
                           double rho) {
 int N = x.n_elem;
 int T = G.n_elem + 1;

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
                 G(t - 1),
                 sigma_n,
                 sigma_u,
                 p_noinc,
                 beta,
                 rho);
 }
 return c;
}

//' Simulate N asset draws from log normal distribution
//' @param N number of draws
//' @param mu mean of distribution
//' @param sigma sd of distribution
//' @export
// [[Rcpp::export]]
arma::vec simulate_assets(int N, double mu, double sigma) {

  // standard normal
  arma::vec assets(N, arma::fill::randn);

  // simple transform
  assets = assets * sigma + mu;

  // assets are lognormal
  assets = exp(assets);

  return assets;
}

//' Draw a random matrix of bernoullis
//' @param p probability of a 0
//' @param N number of rows
//' @param T number of cols
//' @export
//' @details the probability is of a _zero_ not a one
//' because we are generating the probability of zero
//' income and it makes it easier, even though this is
//' non-standard.
// [[Rcpp::export]]
arma::mat draw_bernoulli(double p, int N, int T) {
  arma::mat m(N, T, arma::fill::randu);
  m.for_each([p](arma::mat::elem_type& val) { val = (val >= p); });
  return m;
}

//' Simulate income process
//' @param N number of simulated agents
//' @param T number of time periods
//' @param G vector of growth rates (levels)
//' @param sigma_n standard deviation of permanent income shocks
//' @param sigma_u standard deviation of temporary income shocks
//' @param p_noinc probably income is equal to zero
//' @export
// [[Rcpp::export]]
Rcpp::List simulate_income(
  int N,
  int T,
  arma::vec& P_init, // size N
  arma::vec& G, // size T - 1
  double sigma_n,
  double sigma_u,
  double p_noinc
) {

  // N x T matrix of income
  arma::mat inc_process(N, T);

  // shocks are log normal
  arma::mat n_shocks(N, T, arma::fill::randn);
  arma::mat u_shocks(N, T, arma::fill::randn);

  // scale appropriately
  n_shocks = n_shocks * sigma_n;
  u_shocks = u_shocks * sigma_u;

  // shocks where income is equal to 0
  arma::mat Z_shocks = draw_bernoulli(p_noinc, N, T);

  // need this for case when U = 0
  // otherwise log shock not well defined
  arma::mat U_shocks = exp(u_shocks) % Z_shocks;
  arma::mat N_shocks = exp(n_shocks);

  arma::mat Y(N, T, arma::fill::zeros);
  arma::mat P(N, T, arma::fill::zeros);

  P.col(0) = P_init;
  Y.col(0) = P_init % U_shocks.col(0);

  for(int t = 1; t < T; t++) {
    P.col(t) = P.col(t - 1) * G(t - 1) % N_shocks.col(t);
    Y.col(t) = P.col(t) % U_shocks.col(t);
  }
  return Rcpp::List::create(
    Rcpp::Named("Y") = Y,
    Rcpp::Named("P") = P,
    Rcpp::Named("U_shocks") = U_shocks,
    Rcpp::Named("N_shocks") = N_shocks
  );
}

//' Consume out of cash on hand
//' @param x cash on hand per unit of permanent income
//' @param x_grid grid we are using for linear interpolation
//' @param cr optimal consumption rule
//' solved for with consumption_rule function
//' @details This computes consumption given income
//' assets and return on assets by first computing cash
//' on hand and then passing to the optimal consumption
//' rule.
//' @export
// [[Rcpp::export]]
arma::vec consume(arma::vec x,
                  arma::vec& x_grid,
                  arma::mat& cr,
                  int t) {
  // interpolate on grid and convert back to consumption units
  arma::vec c = linterp(x, x_grid, cr.col(t));
  // enforce borrowing constraint
  c = min(c, x);
  return c;
}

//' Simulate consumption / savings lifecycle problem
//' @param N number of simulations
//' @param T number of time periods
//' @param x_grid grid of permanent income we use for consumption rule
//' @param N_shock matrix of permanent income shocks
//' @param U_shock matrix of temporary income shocks
//' @param G growth of permanent income
//' @param sigma_n sd of permanent income
//' @param sigma_u sd temporary income shocks
//' @param mu_a average log starting assets
//' @param sigma_a sd of log starting assets
//' @param p_noinc probability that income is zero
//' @param R gross return on assets
//' @param beta time discount factor
//' @param rho CRRA risk aversion
//' @export
// [[Rcpp::export]]
Rcpp::List simulate_lifecycle(
  int N,
  int T,
  arma::vec& x_grid,
  arma::mat& N_shock,
  arma::mat& U_shock,
  arma::mat& P,
  arma::vec& init_a,
  arma::vec& G,
  double sigma_n,
  double sigma_u,
  double gamma_0,
  double gamma_1,
  double R,
  double p_noinc,
  double beta,
  double rho
) {

  // consumption rule
  arma::mat cr = consumption_rule(
    x_grid, G,
    sigma_n, sigma_u,
    gamma_0, gamma_1,
    R, p_noinc,
    beta, rho
  );

  // assets
  // always lagged a period relative to income / consumption
  arma::mat x(N, T, arma::fill::zeros);

  // cash on hand constraint in the first period
  x.col(0) = init_a * R / (G(0) * N_shock.col(0)) + U_shock.col(0);

  // consumption over time
  arma::mat c(N, T, arma::fill::zeros);

  for(uint t = 0; t < T; t++) {

    // consumption for each period
    c.col(t) = consume(x.col(t), x_grid, cr, t);

    // assets we have access to next period
    // in last period we don't need to track this
    // this is eq 4 of gourinchas parker
    if(t < T - 1) {
      x.col(t + 1) =
        (x.col(t) - c.col(t)) * R / (G(t) * N_shock.col(t + 1)) +
        U_shock.col(t + 1);
    }
  }

  // column means of simulated consumption
  // this is what enters our loss function
  arma::mat C = c % P;

  arma::vec C_t = mean(C, 0).t();
  arma::vec x_t = mean(x, 0).t();
  return Rcpp::List::create(
    Rcpp::Named("Cbar") = C_t,
    Rcpp::Named("xbar") = x_t,
    Rcpp::Named("C") = C,
    Rcpp::Named("x") = x,
    Rcpp::Named("CR") = cr
  );
}
