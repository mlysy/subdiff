#--- r setup --------------------------------------------------------------

# keep track of labels
if(!file.exists("internal_labels.rds")) {
  saveRDS(list(labels = NULL), "internal_labels.rds")
}
.internal_labels <- readRDS("internal_labels.rds")
.internal_labels$appendix <- FALSE
.internal_labels$section <- NULL
.internal_labels$subsection <- NULL
.internal_labels$subsubsection <- NULL
.internal_labels$figure <- NULL

# link to packages
pkg_link <- function(pkg, link) {
  if(link == "github") {
    link <- paste0("https://github.com/mlysy/", pkg)
  } else if(link == "cran") {
    link <- paste0("https://CRAN.R-project.org/package=", pkg)
  }
  paste0("[**", pkg, "**](", link, ")")
}
cran_link <- function(pkg) pkg_link(pkg, "cran")
github_link <- function(pkg) pkg_link(pkg, "github")

#--- section/figure/table numbering: --------------------------------------

# - counter increments via function
# - label + reference saved in html_label

# section numbering
section <- function(x, label) {
  if(!.internal_labels$appendix) {
    counter <- set_counter("section")
    ref <- counter
  } else {
    counter <- set_counter("section", 1)
    ref <- LETTERS[counter]
  }
  set_counter("subsection", 0)
  set_counter("subsubsection", 0)
  if(missing(label)) label <- default_label(x)
  sec_title <- paste0(ref, " ", x, " {#", label, "}")
  set_label(ref = ref, label = label, anchor = label)
  sec_title
}
subsection <- function(x, label) {
  counter <- set_counter("subsection")
  set_counter("subsubsection", 0)
  if(missing(label)) label <- default_label(x)
  sec_count <- .internal_labels$section
  sec_count <- if(.internal_labels$appendix) LETTERS[sec_count] else sec_count
  ref <- paste0(sec_count, ".", counter)
  sec_title <- paste0(ref, " ", x, " {#", label, "}")
  set_label(ref = ref, label = label, anchor = label)
  sec_title
}
subsubsection <- function(x, label) {
  counter <- set_counter("subsubsection")
  if(missing(label)) label <- default_label(x)
  sec_count <- .internal_labels$section
  sec_count <- if(.internal_labels$appendix) LETTERS[sec_count] else sec_count
  sub_count <- .internal_labels$subsection
  ref <- paste0(sec_count, ".", sub_count, ".", counter)
  sec_title <- paste0(ref, " ", x, " {#", label, "}")
  set_label(ref = ref, label = label, anchor = label)
  sec_title
}
appendix <- function() .internal_labels$appendix <<- TRUE

# figure and table captions
fig_label <- function(x, label) {
  counter <- set_counter("figure")
  if(missing(label)) label <- default_label(x)
  ref <- counter
  set_label(ref = ref, label = label, anchor = label)
  paste0("<center> ", html_label(label), "Figure ", counter, ": ", x, " </center>")
}
tab_label <- function(x, label) {
  counter <- set_counter("table")
  if(missing(label)) label <- default_label(x)
  ref <- counter
  set_label(ref = ref, label = label, anchor = label)
  paste0("<center> ", html_label(label), "Table ", counter, ": ", x, " </center>")
}

# low level label functions
html_label <- function(label) {
  paste0('<a name="', label, '"></a>')
}
set_label <- function(ref, label, anchor) {
  labels <- .internal_labels$labels
  if(label %in% colnames(labels)) {
    labels[,label] <- c(ref = ref, anchor = anchor)
  } else {
    labels <- cbind(c(ref = ref, anchor = anchor), labels)
    colnames(labels)[1] <- label
  }
  .internal_labels$labels <<- labels
  if(!identical(labels, readRDS("internal_labels.rds")$labels)) {
    saveRDS(.internal_labels, file = "internal_labels.rds")
  }
  ## attr(set_label, "labels") <<- c(attr(set_label, "labels"), label = ref)
  invisible(NULL)
}
default_label <- function(x) {
  label <- tolower(x)
  label <- gsub("[[:blank:]]+", "", label)
  gsub("[[:punct:]]+", "", label)
}
ref_label <- function(label) {
  ## ref <- as.character(attr(set_label, "labels")[label])
  labels <- .internal_labels$labels
  if(label %in% colnames(labels)) {
    ref <- labels["ref",label]
    ref <- gsub("([[:digit:]]+[.[:digit:]]*)", "$\\1$", ref)
    anchor <- labels["anchor",label]
  } else {
    ref <- "???"
    anchor <- ""
  }
  paste0("[", ref, "](#", anchor, ")")
}
set_counter <- function(element, n) {
  ## counter <- attr(get(fun_name, envir = globalenv()), "counter")
  if(missing(n)) {
    counter <- .internal_labels[[element]]
    if(is.null(counter)) counter <- 0
    counter <- counter + 1
  } else {
    counter <- n
  }
  ## attr(get(fun_name, envir = globalenv()), "counter") <<- counter
  .internal_labels[[element]] <<- counter
  counter
}
