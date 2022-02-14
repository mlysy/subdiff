#' @return A vector of estimated parameters on the transformed scale. If `vcov == TRUE`, a list with components:
#' \describe{
#' \item{coef}{A vector of estimated parameters on transformed scale.}
#' \item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
#' }
#' If `ad_only == TRUE`, instead of the transformed scale parameters, returns an estimate (and possibly the estimated convariance) of `(alpha, D)`.
