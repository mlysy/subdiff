#' Predefined drift functions for `csi_model` objects.
#'
#' @name drift_functions
#' @aliases drift_linear drift_none drift_quadratic
#'
#' @template args-phi
#' @template args-dt
#' @template args-N
#'
#' @details Models inheriting from the `csi_model` base class must specify a drift function having argument signature `function(phi, dt, N)` and returning an `N x n_drift` matrix of drift basis functions for the trajectory *increments* (as opposed to the original positions).  The predefined drift functions provided here can account for linear drift, quadratic drift, or no drift at all.
NULL

#' @rdname drift_functions
#' @export
drift_linear <- function(phi, dt, N) dt

#' @rdname drift_functions
#' @export
drift_none <- function(phi, dt, N) 0

#' @rdname drift_functions
#' @export
drift_quadratic <- function(phi, dt, N) cbind(dt, diff((0:N*dt)^2))
