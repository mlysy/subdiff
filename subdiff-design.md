## Design of subdiff models


### generic methods

- `summary(model)`: point est (regular scale) standard error (regular scale) CI (asymmetric)
- `coef`: overload estimates
- `sigma`: standard errors
- `confint`: confidence intervals (flags: regular/comp, sym/asym)
- `simulate`: simulate a bunch of paths
- `residuals`: get white noise from process
- `vcov`: computational or human basis? ans: both, so that we can convert back and forth

### `model` object

- acf
- parameter boundaries
- param names?
- computational/regular conversion functions (optional, otherwise determined by bounds)
- gradient of acf (optional)
- penalty (optional)
- fit (optional: function, or perhaps word "optim" or "nlm"?)

### `fitted_model` object

- mle
- variance matrix
- conversion functions for all parameters (incl. error)
- data


### fit method (external)

- interface: something like `csi_fit(model = fbm, data = x, ...)`
- optimization routine:

    - `optimize` (1-d)
    - `optim` (multi d)
    - `nlm` (multi d)

    which of optim/nlm to use by default? maybe depends on gradient?
	
- detrend t/f, 
- hf error removal t/f (use MA(1) model from aoas)
- `control` argument: passes other control arguments to desired fitting method

### other things for package

- documentation: examples, vignette(s) i.e., tutorials, github page
- q > 2 dimensions.  computational basis for q > 3?

###  TODO

- create subrepo with other packages `SuperGauss` and `LMN`.
- create an fbm S3 "object", i.e., a list with elements above.
- create each of the generics.
- create an informal tutorial/vignette for math + code of what we're doing with fbm.
