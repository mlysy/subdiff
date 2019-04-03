#' @return An S3 object of type `csi_class`, consisting of a list wuth elements:
#' \itemize{
#' \item theta_names: A string vector of parameter names.
#' \item acf: ACF function with three arguments: vector of ACF parameters \eqn{\theta}; interobservation time \eqn{\Delta t = 1/fps}; length \eqn{N}.
#' \item theta_trans: Transformation function that converts constrained ACF parameter \eqn{\theta} into unrestricted parameters \eqn{\gamma}.
#' \item theta_itrans: Inverse Transformation function that converts unrestricted parameters \eqn{\gamma} into constrained ACF parameter \eqn{\theta}.
#' \item penalty: Penalty function for unrestricted parameters \eqn{\gamma}.
#' }
