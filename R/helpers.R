#' Create cash on hand grid
#' @param x_max max of the cash on hand grid, 40 in the paper
#' @param x_int internal stopping point for finer grid
#' @param n_points total number of points in the cash on hand grid
#' @details this just follows the paper, a log spaced grid might make
#' more sense
create_coh_grid = function(x_max, x_int, n_points) {
  n_int_points = round(n_points / 2)
  n_ext_points = n_points - n_int_points
  int_grid = seq(0.001, x_int, by = x_int / n_int_points)
  ext_grid =  seq(x_int + x_max / n_ext_points, x_max, by = x_max / n_ext_points)
  c(int_grid, ext_grid)
}

#' Solve model given parameters
#' gamma_0 intercept for retirement consumption
#' gamma_1 consumption in retirement as a linear function of
#' assets
#' beta discount factor
#' rho risk aversion
#' @export
loss = function(gamma_0,
                gamma_1,
                beta,
                rho) {

  out = sim_given_params(
    gamma_0,
    gamma_1,
    beta,
    rho
  )

  lc = log(gp::c / 10000)

  sum((lc - log(out$Cbar))^2)

}

sim_given_params = function(gamma_0,
                            gamma_1,
                            beta,
                            rho) {

  # set based on the paper
  N = 2e4
  T_max = 40
  p_noinc = 0.00302
  sigma_n = 0.0212
  sigma_u = 0.044
  r = 0.0344
  mu_a = -2.794
  sigma_a = 1.784

  x = create_coh_grid(40, x_int = 2, n_points = 100)

  init_p = matrix(rep(18690.96, N), ncol = 1)

  G = exp(g)

  # simulate income realizations conditional on params
  income =
    simulate_income(N, T_max,
                    init_p, G,
                    sigma_n, sigma_u, p_noinc)

  # actual income
  Y = income$Y

  # permanent income
  P = income$P

  # draw initial assets
  init_a = simulate_assets(N, mu_a, sigma_a) / P[,1]

  c_out = simulate_lifecycle(
    N,
    T_max,
    x,
    income$N_shocks,
    income$U_shocks,
    P,
    init_a,
    G,
    sigma_n,
    sigma_u,
    gamma_0,
    gamma_1,
    exp(r),
    p_noinc,
    beta,
    rho
  )
  c_out
}

solve_model = function(init) {
  optim(init, \(theta)  {
    loss(
      theta[1],
      theta[2],
      theta[3],
      theta[4]
    )
  }, lower = c(0, 0, 0.7, 1.3),
  upper = c(0.8, 0.8, 1, 5), method = "L-BFGS-B")
}
