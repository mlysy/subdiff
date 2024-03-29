#' MSD of fBM + Dynamic Error
#'
#' @param alpha Subdiffusion exponent
#' @param tau Width of averaging time-window.
#' @param t Vector of time points \eqn{ \{\Delta t, 2\Delta t, ..., N\Delta t\} }
#' @details 
#' this function returns the MSD of \eqn{Y_t}, the integral of fBM process \eqn{X_t} with subdiffusion 
#' exponent \eqn{\alpha} \deqn{Y_t = \int_{0}^{\tau} X(t-s)ds}. The expression of the MSD is
#' \deqn{\frac{(t+\tau)^\alpha + (t-\tau)^\alpha - 2t^\alpha - 2\tau^\alpha}{(\alpha+1)(\alpha+2)}}
#' @examples 
#' fdyn.msd(alpha = 0.8, tau = 1/600, t = (1:200) * 1/60)
#' @export
fdyn.msd <- function(alpha, tau, t){
  tau <- t/tau
  alpha2 <- alpha+2
  eta <- ((tau+1)^alpha2 + (tau-1)^alpha2 - 2*tau^alpha2 - 2)/alpha2
  eta * tau^alpha/(alpha+1)
}
