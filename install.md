---
title: "Pre-Release Installation Instructions"
---

1.  Install the R package [**Rcpp**](http://www.rcpp.org/).  This may require you to set up your R to compile C++ source code.

2.  Install the following packages from CRAN.  You can do this from within [RStudio](https://www.rstudio.com/), or from the command line within an R session as follows:

    ```r
	install.packages(c("devtools", "testthat", 
	                   "numDeriv", "ggplot2", "SuperGauss"), 
					   INSTALL_opts = "--install-tests")
	```
	
	If on linux, the last of these requires you to first install the the C library [FFTW](http://www.fftw.org/).

3.  Install the following package from GitHub as follows:

    ```r
	devtools::install_github("mlysy/optimCheck", 
	                         INSTALL_opts = "--install-tests")
	```

4.  Install the following packages from source by opening R in the folder where the tarballs are saved, then:

    ```r
    install.packages(c("LMN_0.0.0.9001.tar.gz", "subdiff_0.0.1.9000.tar.gz"),
                     repos = NULL, type = "source", 
					 INSTALL_opts = "--install-tests")
    ```

5.  Test that the installation was successful by running the unit tests:

    ```r
	testthat::test_package("SuperGauss", reporter = "progress")
	testthat::test_package("LMN", reporter = "progress")
	testthat::test_package("subdiff", reporter = "progress")
	```
	
	It's probably safest to quit and restart R after installation before attempting this.  Also, some of the **subdiff** tests may occasionally fail because the optimization routine is not sufficiently close to the mode.  If you repeat the test and it passes (and perhaps a different one fails) then it's likely that everything is OK.
