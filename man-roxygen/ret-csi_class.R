#' @return An list of class `csi_class` with components:
#' \describe{
#' \item{`theta_names`}{A vector of parameter names.}
#' \item{`acf`}{An ACF function with three arguments: vector of parameters; interobservation time \eqn{\Delta t}; length `N`.}
#' \item{`theta_trans`}{Transformation function that converts original ACF parameter into unrestricted parameters.}
#' \item{`theta_itrans`}{Inverse Transformation function that converts unrestricted parameters into original ACF parameter.}
#' \item{`penalty`}{Penalty function for unrestricted parameters.}
#' }
