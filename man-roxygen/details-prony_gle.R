#' @details The Prony-GLE model satisfies the integro-differential equation
#' \deqn{
#' F_t - \int_{-\infty}^t \kappa(t-s) \dot X_s ds = 0, \qquad ACF_F(t) = k_B T \kappa(t),
#' }
#' where \eqn{T} is temperature, \eqn{k_B} is Boltzmann's constant, and the memory kernel is a sum of exponentials:
#' \deqn{
#' \kappa(t) = \nu sum_{k=1}^K exp(-\lambda_k t).
#' }
#' The solution process is of the form
#' \deqn{
#' X_t = C_0 B_t + \sum_{i=1}^{K-1} C_i W_{it},
#' }
#' where \eqn{B_t} is Brownian motion and \eqn{d W_{it} = -\rho_i W_{it} dt + d B_{it}} are Ornstein-Uhlenbeck processes all independent of each other and of \eqn{B_t}.
