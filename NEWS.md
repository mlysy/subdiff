# **subdiff** 0.0.1.9002 (development version)

## Breaking changes

- The argument `dT` in all functions is replaced by `dt`.

- `*_fit()` functions have a completely different argument signature: 
    ```
	function(dX, dT, Tz, var_calc = TRUE) -> function(dX, dt, drift, vcov = TRUE)
	```
	In particular, the Toeplitz object `Tz` can no longer be supplied, one can now choose the type of drift, and `var_calc` is renamed to `vcov`.

