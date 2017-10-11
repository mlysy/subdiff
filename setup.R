#--- setup for Toeplitz package --------------------------------------------

# for recompiling package
# first quit R, then setwd() to where setup.R is found. then:
pkg.path <- getwd()
pkg.path <- "D:/GitHub/SubDiff"

# regenerates Rcpp interface (i.e., RcppExports)
Rcpp::compileAttributes(pkgdir = pkg.path)
devtools::document(pkg = pkg.path)
devtools::install(pkg = pkg.path, args = "--clean") # installs the package
devtools::build(pkg = pkg.path) # builds a tar.gz file

# check
devtools::check("SubDiff")

# restart R before testing changes
testthat::test_package("SubDiff")

# cran check

# generating the pdf manual using rd file
pack <- "SubDiff"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))

build(pkg = pkg.path)

# build windows binary
# NOTE: this will also install the package
build(binary = TRUE)

pkg.path <- "c:/Users/Jerome/Documents/R/test/fftw"
cmd <- file.path(R.home(component = "bin"),
                 paste0("R CMD INSTALL ", pkg.path))
compiled <- system(cmd)

# First Time Installation -------------------------------------------------

require(Rcpp)

pkg.name <- "SubDiff"
pkg.path <- "D:/GitHub/SubDiff"

# for safest results, delete ALL traces of previous package versions
remove.packages(pkg.name)
unlink(file.path(pkg.path, pkg.name), recursive = TRUE)

# minimum number of files to get the package started
Rcpp.package.skeleton(name = pkg.name, path = pkg.path,
                      example_code = FALSE, force = TRUE)

# install package
cmd <- file.path(R.home(component = "bin"), paste0("R CMD INSTALL ", pkg.name))
compiled <- system(cmd)
