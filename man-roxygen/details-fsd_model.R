#' @details Let `X_t` denote the true position of an fBM process at time `t`.  The Savin-Doyle localization error model describes the measured position `Y_t` as
#' ```
#' Y_t = sigma * eps_t + 1/(tau*dt) * int_0^(tau*dt) X_{t - s} ds,
#' ```
#' where `eps_t ~iid N(0,1)` is a Gaussian white noise process.
