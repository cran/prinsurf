#' @export
print.prinsurf <- function(x, ...) {
  cat(sprintf("Principal surface fit: %d samples, %d variables\n",
              nrow(x$lambda), length(x$varnames)))
  cat(sprintf("  loess span %.2f, converged in %d iterations\n", x$span, x$iterations))
  cat(sprintf("  variables: %s\n", paste(x$varnames, collapse = ", ")))
  invisible(x)
}

#' Plot a principal-surface biplot
#'
#' Draws the sample coordinates. If \code{vars} is given, draws one panel per
#' named variable, each showing the sample coordinates together with that
#' variable's contour lines - the surface's analogue of a biplot axis, read by
#' interpolating between contour lines rather than by projection onto a single
#' calibrated axis. With no \code{vars}, a single panel of sample coordinates is
#' drawn with no contours.
#'
#' @param x A \code{"prinsurf"} object.
#' @param vars Variables to contour, one panel each (default: none -- a bare
#'   scatter of sample coordinates).
#' @param group Optional factor to colour the sample points.
#' @param col_contour Colour for the contour lines.
#' @param nlevels Number of contour levels per panel.
#' @param pch,cex Point symbol and size for samples.
#' @param asp Aspect ratio of each panel; the default \code{1} keeps the two
#'   surface coordinates on a common scale, so distances in the biplot are
#'   comparable in every direction.
#' @param main Panel title(s), replacing the default (the variable's name, or no
#'   title when \code{vars} is not given). A single string titles every panel;
#'   a vector titles the panels in the order of \code{vars} and is recycled to
#'   that length. Use \code{""} for no panel titles.
#' @param outer_main A single title for the figure as a whole, drawn in the
#'   outer margin above the panels. Can be combined with \code{main}, which
#'   titles the individual panels.
#' @param ... Passed to each panel's initial \code{plot}.
#' @return Invisibly, \code{vars}.
#' @examples
#' set.seed(1)
#' s <- runif(120, -1, 1); t <- runif(120, -1, 1)
#' X <- cbind(x = s, y = t, z = 0.7 * s + 0.9 * t^2) +
#'      matrix(rnorm(360, 0, 0.03), 120, 3)
#' fit <- prinsurf(X, max.iter = 6)
#'
#' ## default: each panel is titled with its variable's name
#' plot(fit, vars = c("x", "z"))
#'
#' ## your own panel titles, and a title for the figure as a whole
#' plot(fit, vars = c("x", "z"),
#'      main = c("First coordinate", "Third coordinate"),
#'      outer_main = "Principal-surface biplot")
#'
#' ## a title on a single, contour-free panel
#' plot(fit, main = "Sample coordinates")
#' @export
plot.prinsurf <- function(x, vars = NULL, group = NULL,
                          col_contour = "grey40", nlevels = 6,
                          pch = 16, cex = 0.7, asp = 1,
                          main = NULL, outer_main = NULL, ...) {
  lam <- x$lambda
  rng <- apply(lam, 2, range); pad <- 0.28 * (rng[2, ] - rng[1, ])
  cols <- if (is.null(group)) "grey60" else
    grDevices::hcl.colors(nlevels(as.factor(group)), "Dark 3")[as.integer(as.factor(group))]

  npanel <- max(length(vars), 1)
  ## default panel titles are the variable names; user titles are recycled
  main <- if (is.null(main)) {
    if (length(vars)) vars else ""
  } else rep_len(as.character(main), npanel)

  pars <- list()
  if (npanel > 1) {
    nc <- ceiling(sqrt(npanel)); nr <- ceiling(npanel / nc)
    pars$mfrow <- c(nr, nc)
  }
  ## room above the panels for the figure title
  if (!is.null(outer_main)) pars$oma <- c(0, 0, 3, 0)
  if (length(pars)) { op <- do.call(graphics::par, pars); on.exit(graphics::par(op)) }

  for (i in seq_len(npanel)) {
    v <- if (length(vars)) vars[i] else NULL
    graphics::plot(lam, pch = pch, cex = cex, col = cols, axes = FALSE, asp = asp,
                   xlim = rng[, 1] + c(-pad[1], pad[1]), ylim = rng[, 2] + c(-pad[2], pad[2]),
                   xlab = expression(lambda[1]), ylab = expression(lambda[2]),
                   main = main[i], ...)
    graphics::box()
    if (!is.null(v)) {
      g <- .ps_grid(x, v)
      graphics::contour(g$g1, g$g2, g$M, add = TRUE, col = col_contour, nlevels = nlevels,
                        lwd = 0.7, labcex = 0.6)
    }
    if (!is.null(group) && i == 1)
      graphics::legend("topright", levels(as.factor(group)),
                       col = grDevices::hcl.colors(nlevels(as.factor(group)), "Dark 3"),
                       pch = pch, bty = "n", cex = 0.7)
  }
  if (!is.null(outer_main))
    graphics::title(main = outer_main, outer = TRUE, cex.main = 1.3)
  invisible(vars)
}

#' Fitted (reconstructed) values from a principal surface
#'
#' Returns the fitted surface values \eqn{\hat f_j(\lambda_i)} for every sample
#' and variable: the point on the surface at which sample \eqn{i} sits, written
#' back in the coordinates of the original variables. Row \eqn{i} is the
#' surface's reconstruction of sample \eqn{i} -- what the fit says the sample
#' would be if it lay exactly on the surface -- and the residual
#' \eqn{x_i - \hat f(\lambda_i)} is the part of the sample the surface does not
#' capture, which is what \code{\link{predictivity}} summarises.
#'
#' Values are in the working units used for fitting: centred, and divided by
#' each variable's standard deviation if the surface was fitted with
#' \code{scale = TRUE}. This is the scale on which residuals and
#' \code{\link{predictivity}} are computed. \code{\link{predict.prinsurf}}
#' differs in two ways: it reads values off the plotted contour grid rather than
#' evaluating the coordinate functions exactly, and it returns them on the
#' variables' original scales.
#' @param object A \code{"prinsurf"} object.
#' @param ... Ignored.
#' @return An \eqn{n \times p} matrix of fitted values, in working
#'   (centred, optionally scaled) units.
#' @seealso \code{\link{predictivity}} for the per-sample quality of this
#'   reconstruction, and \code{\link{predict.prinsurf}} for the values a reader
#'   obtains from the contours.
#' @export
fitted.prinsurf <- function(object, ...) {
  p <- length(object$varnames)
  out <- vapply(seq_len(p), function(j) .ps_eval(object, object$lambda, j),
                numeric(nrow(object$lambda)))
  colnames(out) <- object$varnames
  out
}

#' Predict all variables from the biplot contours
#'
#' Reads every variable's value for every sample off that variable's contour
#' lines at the sample's biplot position \eqn{\lambda_i} -- the same
#' interpolation used to draw the contours in \code{\link{plot.prinsurf}}.
#' Values come from the contour grid alone: a sample whose position is not
#' covered by the supported part of the grid has no contours to read, and is
#' returned as \code{NA} for every variable.
#' @param object A \code{"prinsurf"} object.
#' @param ... Ignored.
#' @return An \eqn{n \times p} matrix of values read from the contours, on the
#'   variables' original scales, with \code{NA} rows where the biplot cannot be
#'   read.
#' @seealso \code{\link{contour_predictive_error}} to measure this reading
#'   against the variables' actual values, and \code{\link{fitted.prinsurf}}
#'   for the surface's own reconstruction of the data.
#' @export
predict.prinsurf <- function(object, ...) {
  pred <- .contour_read(object, object$lambda)
  rownames(pred) <- rownames(object$X)
  ## return on the variables' original scales (undo the fit-time centring/scaling)
  sweep(sweep(pred, 2, object$scale, "*"), 2, object$center, "+")
}
