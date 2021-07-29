# **subdiff** 0.0.1.9002 (development version)

## Breaking changes

- The argument `dT` in all functions is replaced by `dt`.

- `*_fit()` functions have a completely different argument signature: 
    ```
	function(dX, dT, Tz, var_calc = TRUE) -> function(dX, dt, drift, vcov = TRUE)
	```
	In particular, the Toeplitz object `Tz` can no longer be supplied, one can now choose the type of drift, and `var_calc` is renamed to `vcov`.

- Model objects (`fbm_model`, `farma_model`, etc) now use R6 instead of S3 classes and inherit from base class `csi_model`.  Thus they are instantiated via `$new()` instead of just e.g., `fbm_model()`.

- Fitting functions (`fbm_fit()`, etc.) output parameter names in a (standardized) computational basis. 

- The `floc` prefix (`floc_model()`, `floc_acf()`, `floc_fit()`, etc.) has been changed to `fsd`, which stands for fractional Savin-Doyle.

- The functions `trans_Sigma()` and `itrans_Sigma()` have different meanings.

- `csi_resid()` no longer has arguments `mu` and `dt`, in favor of an arbitrary drift term supplied by `drift`.

## New features

- Model objects are formalized via R6 classes.  The main purpose is to provide a simple and flexible framework for adding user-defined models, parameter estimation routines, etc.

- `farma_model` (and `farma_acf()`, etc) now supports arbitrary order for the autoregressive component.

- `ls_fit()` now provides standard errors for `alpha` and `logD`.  Also, drift subtraction is performed much more accurately via linear regression than as previously via mean increment value.

- `*_fit()` for subdiffusion models gives the option of returning estimates (and standard errors) for `alpha` and `logD` only.
