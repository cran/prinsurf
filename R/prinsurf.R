#' Fit a principal surface
#'
#' Fits a two-dimensional principal surface to a numeric data matrix using the
#' iterative expectation / projection algorithm of Hastie and Stuetzle (1989),
#' with the coordinate functions estimated by local regression (\code{loess}).
#'
#' @param X A numeric matrix or data frame, \eqn{n \times p}. Column names are
#'   used as variable names; if absent, \code{V1, V2, ...} are assigned.
#' @param max.iter Maximum number of expectation/projection iterations.
#' @param span The \code{loess} span (\eqn{\alpha}) used for the coordinate
#'   functions.
#' @param scale Logical; if \code{TRUE} each variable is standardised to unit
#'   standard deviation before fitting (recommended when variables are on very
#'   different scales). Variables are always centred.
#' @param verbose Logical; print the relative change and residual sum of
#'   squares at each iteration.
#'
#' @return An object of class \code{"prinsurf"}: a list with elements
#'   \code{lambda} (the \eqn{n \times 2} surface coordinates of the samples),
#'   \code{fj.mat} (the fitted surface coordinates of the samples in the working
#'   units), \code{models} (the per-variable \code{loess} coordinate functions),
#'   \code{varnames}, \code{center}, \code{scale}, \code{span} and
#'   \code{iterations}.
#'
#' @references Hastie, T. and Stuetzle, W. (1989) Principal curves.
#'   \emph{Journal of the American Statistical Association} \strong{84}, 502--516.
#'
#' @examples
#' set.seed(1)
#' s <- runif(120, -1, 1); t <- runif(120, -1, 1)
#' X <- cbind(x = s, y = t, z = 0.7 * s + 0.9 * t^2) +
#'      matrix(rnorm(360, 0, 0.03), 120, 3)
#' fit <- prinsurf(X, max.iter = 6)
#' fit
#' @export
prinsurf <- function(X, max.iter = 10, span = 0.6, scale = FALSE,
                     verbose = FALSE) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (ncol(X) < 3) stop("A principal surface needs at least 3 variables.")
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  vn <- colnames(X)
  center <- colMeans(X)
  scl <- if (scale) apply(X, 2, stats::sd) else rep(1, ncol(X))
  if (any(scl == 0)) stop("A variable has zero variance.")
  Xs <- sweep(sweep(X, 2, center), 2, scl, "/")
  fit <- .ps_fit(Xs, max.iter = max.iter, alpha = span, verbose = verbose)
  structure(list(
    call = match.call(),
    X = Xs, varnames = vn, center = center, scale = scl,
    lambda = fit$lambda, fj.mat = fit$fj.mat, models = fit$models,
    span = span, iterations = fit$iter
  ), class = "prinsurf")
}

## clamp a scalar to [0, 1], treating NA as 0
.clamp01 <- function(x) { if (is.na(x)) return(0); min(max(x, 0), 1) }

## core fit: returns lambda, fj.mat, per-variable loess models, iteration count
.ps_fit <- function(X, max.iter = 10, alpha = 0.6, verbose = FALSE) {
  n <- nrow(X); p <- ncol(X)
  X0  <- scale(X, scale = FALSE)
  SVD <- svd(X0)
  fj  <- X0 %*% SVD$v[, 1:2] %*% t(SVD$v[, 1:2]) +
         matrix(1, n, 1) %*% colMeans(X)                 # initial plane f(0)
  SVD2 <- svd(fj); lambda <- fj %*% SVD2$v[, 1:2]
  sumD <- sum((fj - X)^2); count <- 0L; finish <- FALSE
  while (!finish && count < max.iter) {
    ## expectation step: coordinate functions are loess smooths of X on lambda
    for (j in seq_len(p))
      fj[, j] <- stats::fitted(
        stats::loess(X[, j] ~ lambda[, 1] + lambda[, 2], span = alpha))
    ## projection step: nearest tangent-plane projection of each sample.
    ## The two nearest surface neighbours of each point do not depend on the
    ## sample being projected, so they are computed once per iteration.
    D <- as.matrix(stats::dist(fj)); diag(D) <- max(D)
    nbr <- t(apply(D, 1, function(r) order(r)[1:2]))
    B1 <- fj[nbr[, 1], ] - fj; B2 <- fj[nbr[, 2], ] - fj
    PP2 <- matrix(0, n, p)
    for (i in seq_len(n)) {
      PP <- matrix(0, n, p)
      for (k in seq_len(n)) {
        b1 <- B1[k, ]; b2 <- B2[k, ]; a <- X[i, ] - fj[k, ]
        b1b1 <- sum(b1 * b1); b2b2 <- sum(b2 * b2); b1b2 <- sum(b1 * b2)
        ab1 <- sum(a * b1);   ab2 <- sum(a * b2)
        dt <- b1b1 * b2b2 - b1b2^2
        c1 <- .clamp01((ab1 * b2b2 - ab2 * b1b2) / dt)
        d1 <- .clamp01((ab2 * b1b1 - ab1 * b1b2) / dt)
        PP[k, ] <- c1 * b1 + d1 * b2 + fj[k, ]
      }
      PP2[i, ] <- PP[which.min(rowSums(sweep(PP, 2, X[i, ])^2)), ]
    }
    fj <- PP2
    SVD2 <- svd(fj); lambda <- fj %*% SVD2$v[, 1:2]
    sn <- sum((fj - X)^2); e1 <- abs(sn - sumD) / sumD
    if (e1 < 1e-3) finish <- TRUE
    sumD <- sn; count <- count + 1L
    if (verbose)
      message(sprintf("  iter %d  rel.change %.5f  sumD %.3f", count, e1, sumD))
  }
  dff <- data.frame(l1 = lambda[, 1], l2 = lambda[, 2])
  models <- lapply(seq_len(p), function(j) {
    d <- dff; d$y <- X[, j]
    stats::loess(y ~ l1 + l2, data = d, span = alpha)
  })
  list(lambda = lambda, fj.mat = fj, models = models, iter = count)
}

## evaluate the fitted coordinate function of variable j at 2-D coordinates L
.ps_eval <- function(object, L, j) {
  L <- matrix(L, ncol = 2)
  stats::predict(object$models[[j]],
                 data.frame(l1 = L[, 1], l2 = L[, 2]))
}

## resolve a variable given as name or index to an integer column
.ps_var <- function(object, var) {
  if (is.character(var)) {
    v <- match(var, object$varnames)
    if (is.na(v)) stop("Unknown variable: ", var)
    v
  } else {
    v <- as.integer(var)
    if (v < 1 || v > length(object$varnames)) stop("Variable index out of range.")
    v
  }
}
