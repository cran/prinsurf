## Evaluation-grid axes over the range of the sample coordinates.
.ps_axes <- function(lam, N) {
  list(g1 = seq(min(lam[, 1]), max(lam[, 1]), length.out = N),
       g2 = seq(min(lam[, 2]), max(lam[, 2]), length.out = N))
}

## TRUE for grid nodes further than a neighbour-density-based radius from every
## sample, i.e. outside the support of the data. The mask depends only on
## lambda, not on the variable, so it is computed once and shared by all
## variables' grids.
.ps_unsupported <- function(lam, grid) {
  rad <- 2.5 * sqrt(diff(range(lam[, 1]))^2 + diff(range(lam[, 2]))^2) / sqrt(nrow(lam))
  nn <- apply(grid, 1, function(q) sqrt(min(colSums((t(lam) - q)^2))))
  nn >= rad
}

## Evaluation grid of one variable's fitted coordinate function, used to draw
## that variable's contours in plot.prinsurf(). Unsupported nodes are left NA so
## the contours stop at the support of the data.
.ps_grid <- function(object, var, N = 55) {
  VAR <- .ps_var(object, var)
  lam <- object$lambda
  ax <- .ps_axes(lam, N)
  grid <- as.matrix(expand.grid(l1 = ax$g1, l2 = ax$g2))
  fg <- .ps_eval(object, grid, VAR)
  fg[.ps_unsupported(lam, grid)] <- NA
  list(g1 = ax$g1, g2 = ax$g2, M = matrix(fg, N, N))
}

## Bilinear interpolation of a grid Mat (N x N, over g1 x g2) at point (x, y).
## NA outside the grid range or where a corner node is unsupported (NA).
.ps_bilin <- function(Mat, x, y, g1, g2) {
  N <- length(g1)
  if (is.na(x) || x < g1[1] || x > g1[N] || y < g2[1] || y > g2[N]) return(NA)
  i <- min(max(findInterval(x, g1), 1), N - 1); k <- min(max(findInterval(y, g2), 1), N - 1)
  ax <- (x - g1[i]) / (g1[i + 1] - g1[i]); ay <- (y - g2[k]) / (g2[k + 1] - g2[k])
  (1 - ax) * (1 - ay) * Mat[i, k] + ax * (1 - ay) * Mat[i + 1, k] +
    (1 - ax) * ay * Mat[i, k + 1] + ax * ay * Mat[i + 1, k + 1]
}

## Read every variable's value at 2-D coordinates L by bilinear interpolation of
## the same evaluation grids drawn as contours in plot.prinsurf() -- i.e. the
## values a reader would obtain by interpolating between the printed contour
## lines. Reading is strictly from the contour grid: positions whose surrounding
## grid nodes are unsupported are NA, since there are no contours to read there.
## The support mask is the same for every variable, so NAs occupy whole rows.
.contour_read <- function(object, L, N = 55) {
  lam <- object$lambda
  ax <- .ps_axes(lam, N)
  grid <- as.matrix(expand.grid(l1 = ax$g1, l2 = ax$g2))
  bad <- .ps_unsupported(lam, grid)
  L <- matrix(L, ncol = 2)
  out <- matrix(NA_real_, nrow(L), length(object$varnames),
                dimnames = list(rownames(L), object$varnames))
  for (j in seq_along(object$varnames)) {
    fg <- .ps_eval(object, grid, j); fg[bad] <- NA
    M <- matrix(fg, N, N)
    out[, j] <- vapply(seq_len(nrow(L)),
                       function(i) .ps_bilin(M, L[i, 1], L[i, 2], ax$g1, ax$g2),
                       numeric(1))
  }
  out
}
