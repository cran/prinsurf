#' Sample predictivity of a principal-surface biplot
#'
#' For each sample, the proportion of its squared length reconstructed by the
#' fitted surface, \eqn{1 - \lVert x_i - \hat f(\lambda_i) \rVert^2 /
#' \lVert x_i \rVert^2}. This is the principal-surface analogue of biplot sample
#' predictivity and is read from the sample's position on the surface.
#' @param object A \code{"prinsurf"} object.
#' @return A numeric vector of per-sample predictivities (with the overall mean as
#'   the attribute \code{"overall"}).
#' @export
predictivity <- function(object) {
  X <- object$X                       # working (centred / scaled) data
  Xhat <- stats::fitted(object)
  Xc <- scale(X, scale = FALSE)
  sst <- rowSums(Xc^2)
  pred <- 1 - rowSums((X - Xhat)^2) / sst
  attr(pred, "overall") <- mean(pred)
  pred
}

#' Contour reading error
#'
#' For each variable, the root-mean-square difference between the value read
#' from the variable's contour lines at each sample's biplot position (as
#' \code{\link{predict.prinsurf}} does) and the sample's actual value. This
#' measures how accurately a variable can be recovered by reading its contours,
#' combining the surface's lack of fit for that variable with the small loss
#' from interpolating between contour lines; lower is better, in the units of
#' the working data. Samples that fall outside the supported part of the
#' contour grid cannot be read and are excluded; their number is returned as the
#' attribute \code{"n.unread"}.
#' @param object A \code{"prinsurf"} object.
#' @return A named numeric vector of per-variable RMS contour-reading errors,
#'   with the mean over all variables as the attribute \code{"overall"} and the
#'   number of unreadable samples as the attribute \code{"n.unread"}.
#' @seealso \code{\link{predictivity}} for whole-sample reconstruction.
#' @export
contour_predictive_error <- function(object) {
  vn <- object$varnames
  pred <- .contour_read(object, object$lambda)
  out <- vapply(seq_along(vn), function(j)
    sqrt(mean((object$X[, j] - pred[, j])^2, na.rm = TRUE)), numeric(1))
  names(out) <- vn
  attr(out, "overall") <- mean(out)
  ## the support mask is shared by all variables, so NAs occupy whole rows
  attr(out, "n.unread") <- sum(is.na(pred[, 1]))
  out
}
