# Tools for Subdiffusion Modeling in Passive Particle-Tracking Rheology

## Overview of functions

### ACF Fitting

For the increments of parametric Gaussian CSI location-scale models.

- fBM model: `fbm_fit`, `fbm_acf`
- fARMA(p,q) model: `farma_fit`, `farma_acf`. Currenlty only `p = 0/1` implemented (but arbitary q).  Useful special case is `fma_fit` which is fARMA(0,1).
- Savin-Doyle Localization Error model: `floc_fit`, `floc2_fit` has different penalty.  `fdyn_acf` is fBM + dynamic error, since adding static error is trivial.
- Downsampling: `fds_fit`.  Can probably be depreciated...
- Helper functions:
    - `Tz_fit`: this does the optimization, but can perhaps be refactored once models are considered to be "objects".
    - `trans_{alpha,rho,tau,Sigma}`: parameter transformations.  These should probably be named by what they do, instead of the parameter they refer to.

### Semi-Parametric Fitting

Based on first estimating the empirical MSD of a CSI model.xs

- `msd_fit`: compute the empirical MSD.
- `msd_ls`: least-squares MSD estimator from empirical MSD.
- `msd_subdiff`: estimates an "effective subdiffusion" time window by defining an objective function between the true MSD and a straight line on the log-log scale.  Also returns the effective subdiffusion parameters.

### GLE Model

No fitting functions for this yet, but could easily be done for derived ACF/MSD of a couple of variants.

- Discrete Rouse-GLE: `rouse_acf`, `rouse_msd`.
    - `rouse_sub` is an "effective subdiffusion" estimator sho shojld probably use `msd_subdiff`.
    - Helper functions: `prony_coef`, `prony_msd`, `prony_acf` for arbitrary Prony series, with Rouse specificying a power lway coefficient spacing.  `mode_poly` for polynomial mode-finding involved.
- Continuous Rouse-GLE: `crouse_msd`.

### Goodness-of-Fit

Tools for diagnosing model fit.

- Model residuals: `fbm_resid`, `fma_resid`, `far_resid` (depreciated), `fsd_resid` (to be soon depreciated), `fd_resid` (rename to `floc_resid`).  The workhorse function is `csi_resid` which applies to any Gaussian CSI location-scale model.
- P-values for normal residuals: `ad_pval`, `sw_pval`, `bc_pval`.  Note that `ad_pval` and `bc_pval` assume that we are testing against N(0,1), whereas `sw_pval` is against N(mu,sig^2), which is not quite correct.
