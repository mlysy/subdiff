# subdiff: Subdiffusive Modeling in Passive Particle-Tracking Microrheology

*Martin Lysy, Yun Ling*

<!-- badges: start -->
[![R-CMD-check](https://github.com/mlysy/subdiff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mlysy/subdiff/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

---

### Description

Tools for implementing various models for particle subdiffusion in biological fluids.  In addition to the well-known semiparametric least-squares estimator based on the mean square displacement (MSD), the package provides functions for simulation, inference, and goodness-of-fit for two fully-parametric subdiffusion models: fractional Brownian motion and the generalized Langevin equation (GLE) model with Rouse memory kernel.  A generic model class allows users to easily implement custom subdiffusion models, with a ready-made framework for drift, high-frequency error correction, and efficient maximum likelihood estimation.

### Installation

1.  Install the [**Rcpp**](https://CRAN.R-project.org/package=Rcpp) package.  To do this properly you will also need to install a C++ compiler.  Instructions for this are available [here](https://teuder.github.io/rcpp4everyone_en/020_install.html).  

    To test that this step is done correctly, quit + restart R, then run the following command:
	
    ```r
    Rcpp::evalCpp("2 + 2")
    ```

    If the output is `4` then **Rcpp** is installed correctly.

2.  Install the [**devtools**](https://CRAN.R-project.org/package=devtools) package.

3.  In an R session run the following commands:

    ```r
    devtools::install_github("mlysy/subdiff",
                             force = TRUE, INSTALL_opts = "--install-tests")
    ```

4.  Once the packages are installed, you can test that everything works correctly by first quitting + restarting R, then running the command:

    ```r
    testthat::test_package("subdiff", reporter = "progress")
    ```
	
	Occasionally due to random chance a few of the tests may fail.  However, if everything is installed correctly, then rerunning this a few times will eventually produce zero test failures.  

### Usage

Please see the package vignette: `vignette("subdiff")`.

