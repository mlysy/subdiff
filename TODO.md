# TODO List

- [x] Change `dX` argument to `Xt`.

- [x] Update examples to new API.

- [x] `farma_fit()` and `farma_model$new()` specify order via `order` and `p,q` respectively.  Use `order` to be consistent with `stats::arima()`.

- [ ] Adopt the usual conventions `MSD(t) = sum E[(X_i(t) - X_i(0))^2]` and `MSD(t) = 2d D t^alpha`.

    Files/functions affected: vignette, `csi_model$msd()` (documentation update), `{fbm/farma/fsd}_model$get_subdiff()`, `ls_fit()` and `msd_fit()`.  Also makes sense not to standardize in `src/msd_empirical.cpp`.  

- [x] Clean up `msd_fit()`, `ls_fit()`, and `msd_ls()`.  Proposal is to remove the latter (put in `inst/proj`), and instead have `ls_fit()` optionally accept `msd` argument if these have already been calculated by `msd_fit()`.  `ls_fit()` can still be given weights though.

	Update: No weights, but now providing `msd_fit()`, `ls_fit()`, and `ls_msd_fit()`, where the latter uses the MSD directly.  Also, `msd_ls()` is an internal function because it is still used by `msd_subdiff()`.  See below.
	
- [x] Standardize unconstrained parametrization for variances.

- [x] Fix `.get_nq` for values of `q > 2`.

- [ ] Create methods for S3 model objects.

- [x] Check over documentation, including function and argument names.  In particular, `order` argument to `farma` models is inconsistent...

- [ ] Remove `old` from repo (version controlled anyways).

- [x] Fix `*_resid()` functions.  Currently using depreciated `trans_alpha`, etc.

- [ ] Remove dependency of `msd_subdiff()` on `msd_ls()`.

- [ ] `msd_fit()` and `msd_subdiff()` should compare true MSD to estimate.  Also, should use GLE model.

- [ ] Check `rouse_sub()`.

- [x] Use regression method for `msd_fit(demean = TRUE)`.

# Formatting

- [ ] Use roxygen markdown for all documentation.

- [ ] Check all documentation for consistent use of sentence case.

- [ ] Explicitly name all arguments in function calls instead of simply relying on order.

- [ ] Format references to `@details` section as e.g., `(see 'Details')`.  The point is that Details should be between single quotes as opposed to double-stars (i.e., bold).  Same for referring to other sections.

- [ ] Fix/create examples for `csi_model`, `fbm_model`, etc.

