
# helper function for histogram
hist_plot <- function(Z, xlab, main = "", lgd) {
  hist(Z, breaks = 25, freq = FALSE, main = main,
       xlab = "", ylab = "")
  curve(dnorm, add = TRUE, col = "red")
  title(xlab = xlab, ylab = "Density", line = 2.5)
  if(!missing(lgd)) legend("topright", legend = parse(text = lgd))
}

# helper function for qq plot
qq_plot <- function(Z, main = "", lgd) {
  qq <- qqnorm(Z, plot.it = FALSE)
  plot(0, type = "n", xlim = range(qq$x), ylim = range(qq$y),
       xlab = "", ylab = "", main = main)
  abline(a = 0, b = 1, lty = 2, col = "blue")
  points(qq$x, qq$y, pch = 16, cex = .5)
  title(xlab = "Theoretical Quantiles",
        ylab = "Sample Quantiles", line = 2.5)
  if(!missing(lgd)) legend("bottomright", legend = parse(text = lgd))
}

# multiple lines on same plot.
# essentially takes range first to draw an empty plots with correct limits.
.multi_plot <- function (x, y, xlim, ylim, add = FALSE, axes = TRUE,
                        log = "", ...) {
  # plot limits
  if(missing(x)) x <- 1:nrow(y)
  if(missing(xlim)) xlim <- range(x, na.rm = TRUE, finite = TRUE)
  if(missing(ylim)) ylim <- range(y, na.rm = TRUE, finite = TRUE)
  # new plot
  if(!add) {
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = log)
  }
  # axes
  if(axes) {
    box()
    axis(side = 1)
    axis(side = 2)
  }
  # content
  invisible(do.call(mapply, args = c(list(FUN = lines,
                                          x = as.data.frame(x),
                                          y = as.data.frame(y)),
                                     list(...))))
}

#' Add textbox to plot.
#'
#' Behaves exactly like \code{\link[graphics]{legend}}, except without symbols.
#' @param x, y, legend, bty, bg, box.lwd, box.lty, box.col, cex, xjust, yjust, x.intersp, y.intersp, adj, text.width, text.col, text.font, trace, plot, ncol, horiz, title, inset, xpd, title.col, title.adj All arguments behave exactly as \code{\link[graphics]{legend}} counterpart.
#' @return See \code{\link[graphics]{legend}}.
#' @details The implementation is essentially \code{legend(fill = NULL, ...}}
textbox <- function(x, y = NULL, legend,
                    bty = "o", bg = par("bg"),
                    box.lwd = par("lwd"), box.lty = par("lty"),
                    box.col = par("fg"), cex = 1,
                    xjust = 0, yjust = 1,
                    x.intersp = 1, y.intersp = 1, adj = c(0, 0.5),
                    text.width = NULL, text.col = par("col"), text.font = NULL,
                    trace = FALSE, plot = TRUE,
                    ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd,
                    title.col = text.col, title.adj = 0.5) {
  if (missing(legend) && !missing(y) && (is.character(y) ||
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  title <- as.graphicsAnnot(title)
  if (length(title) > 1) stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  n.leg <- if (is.call(legend)) 1 else length(legend)
  if (n.leg == 0) stop("'legend' is of length 0")
  auto <- if (is.character(x)) {
            match.arg(x, c("bottomright", "bottom", "bottomleft",
                           "left", "topleft", "top", "topright",
                           "right", "center"))
          } else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y, setLab = FALSE)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) stop("invalid coordinate lengths")
  } else nx <- 0
  xlog <- par("xlog")
  ylog <- par("ylog")
  rect2 <- function(left, top, dx, dy, ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    text(x, y, ...)
  }
  if (trace)
    catn <- function(...) do.call("cat", c(lapply(list(...),
                                                  formatC), list("\n")))
  cin <- par("cin")
  Cex <- cex * par("cex")
  if (is.null(text.width)) {
    text.width <- max(abs(strwidth(legend, units = "user",
                                   cex = cex, font = text.font)))
  } else if (!is.numeric(text.width) || text.width < 0) {
    stop("'text.width' must be numeric, >= 0")
  }
  xc <- Cex * xinch(cin[1L], warn.log = FALSE)
  yc <- Cex * yinch(cin[2L], warn.log = FALSE)
  if (xc < 0) text.width <- -text.width
  xchar <- xc
  xextra <- 0
  yextra <- yc * (y.intersp - 1)
  ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
  ychar <- yextra + ymax
  if (trace) catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, ychar))
  n.legpercol <- if (horiz) {
                   if (ncol != 1)
                     warning(gettextf("horizontal specification overrides: Number of columns := %d",
                                      n.leg), domain = NA)
                   ncol <- n.leg
                   1
                 } else ceiling(n.leg/ncol)
  if (is.na(auto)) {
    if (xlog)
      x <- log10(x)
    if (ylog)
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    wd <- diff(x)
    ht <- diff(y)
    w0 <- wd/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust))
      xjust <- 0.5
    if (missing(yjust))
      yjust <- 0.5
  } else {
    ## ht <- (n.legpercol + (!is.null(title))) * ychar + yc
    ## w0 <- text.width + (x.intersp + 1) * xchar
    ht <- (n.legpercol + (!is.null(title))) * ychar + .4*yc
    w0 <- text.width + .5*x.intersp * xchar
    wd <- ncol * w0 + 0.5 * xchar
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user",
                                               cex = cex) + 0.5 * xchar)) > abs(wd)) {
      xextra <- (tw - wd)/2
      wd <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * wd
      top <- y + (1 - yjust) * ht
    } else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = ,
                     right = usr[2L] - wd - insetx, bottomleft = ,
                     left = , topleft = usr[1L] + insetx, bottom = ,
                     top = , center = (usr[1L] + usr[2L] - wd)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = ,
                    bottomleft = usr[3L] + ht + insety, topleft = , top = ,
                    topright = usr[4L] - insety, left = , right = ,
                    center = (usr[3L] + usr[4L] + ht)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace)
      cat("  rect2(", left, ",", top, ", w=", wd, ", h=",
           ht, ", ...)\n", sep = "")
    rect2(left, top, dx = wd, dy = ht, col = bg, density = NULL,
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  ## xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1),
  ##                                             rep.int(n.legpercol, ncol)))[1L:n.leg]
  xt <- left + (w0 * rep.int(0:(ncol - 1),
                             rep.int(n.legpercol, ncol)))[1L:n.leg]
  yt <- top - 0.5 * yextra - .75*ymax - (rep.int(1L:n.legpercol, ncol)[1L:n.leg] -
                                     1 + (!is.null(title))) * ychar
  xt <- xt + .5 * xchar
  ## xt <- xt + x.intersp * xchar
  if (plot) {
    if (!is.null(title))
      text2(left + wd * title.adj, top - ymax, labels = title,
            adj = c(title.adj, 0), cex = cex, col = title.col)
    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = wd, h = ht, left = left, top = top),
                 text = list(x = xt, y = yt)))
}


