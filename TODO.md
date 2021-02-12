# TODO List

- [ ] Change `dX` argument to `Xt`.

- [ ] Adopt the usual conventions `MSD(t) = sum E[(X_i(t) - X_i(0))^2]` and `MSD(t) = 2d D t^alpha`.

- [ ] Clean up `msd_fit()`, `ls_fit()`, and `msd_ls()`.  Proposal is to remove the latter (put in `inst/proj`), and instead have `ls_fit()` optionally accept `msd` argument if these have already been calculated by `msd_fit()`.  `ls_fit()` can still be given weights though.

- [x] Standardize unconstrained parametrization for variances.

- [x] Fix `.get_nq` for values of `q > 2`.

- [ ] Create methods for S3 model objects.

- [x] Check over documentation, including function and argument names.  In particular, `order` argument to `farma` models is inconsistent...

- [ ] Remove `old` from repo (version controlled anyways).

- [x] Fix `*_resid()` functions.  Currently using depreciated `trans_alpha`, etc.

- [ ] `msd_fit()` and `msd_subdiff()` should compare true MSD to estimate.  Also, should use GLE model.

- [ ] Check `rouse_sub()`.

- [ ] Use regression method for `msd_fit(demean = TRUE)`.

# Formatting

- [ ] Use roxygen markdown for all documentation.

- [ ] Check all documentation for consistent use of sentence case.

- [ ] Explicitly name all arguments in function calls instead of simply relying on order.

- [ ] Format references to `@details` section as e.g., `(see 'Details')`.  The point is that Details should be between single quotes as opposed to double-stars (i.e., bold).  Same for referring to other sections.

- [ ] Fix/create examples for `csi_model`, `fbm_model`, etc.

