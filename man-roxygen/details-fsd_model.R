#' @details Let `X_t` denote the true position of an fBM process at time `t`.  The Savin-Doyle localization error model describes the measured position `Y_n` at time `t = n dt` as
#' ```
#' Y_n = sigma * eps_n + 1/(tau*dt) * int_0^(tau*dt) X_{n dt - s} ds,
#' ```
#' where `eps_n ~iid N(0,1)` is a Gaussian white noise process.  The resulting MSD as a function of `tau` and `sigma2 = sigma^2` is what gets passed to the [csi_model] base class to construct the `fsd_model` derived class.
