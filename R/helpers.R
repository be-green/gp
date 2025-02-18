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
