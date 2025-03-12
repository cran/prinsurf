#' Principal Surface
#'
#' A function to compute principal surfaces based on input data containing continuous variables.
#'
#' @param X A data frame or matrix containing continuous variables.
#' @param max.iter Integer. Maximum number of iterations for the principal surface algorithm.
#' @param alpha Numeric. The span argument passed to the `loess()` function.
#' @param N Integer. The resolution for the interpolated grid surface, creating an \eqn{N^2 \times p} matrix.
#' @param print_iterations Logical. Should the iterations in the principal surface algorithm be printed? Defaults to `FALSE`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{fj.mat}}{A numeric \eqn{n \times p} matrix of the final principal surface fitted values.}
#'   \item{\code{lambda.j}}{A numeric representation of the samples in two dimensions.}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' surface <- principal.surface(iris[,1:3],max.iter = 3)}
#' surface <- principal.surface(iris[1:50,1:3],max.iter = 3)
principal.surface <- function(X, max.iter = 10, alpha = 0.6, N=50,print_iterations=FALSE)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  X0 <- scale(X, scale = F)
  SVD <- svd(X0)

  fj.mat <- X0 %*% SVD$v[, 1:2] %*% t(SVD$v[, 1:2]) + matrix(1,
                                                             ncol = 1, nrow = n) %*% apply(X, 2, mean)
  f0.mat <- fj.mat
  SVD2 <- svd(fj.mat)
  lambda.j <- fj.mat %*% SVD2$v[, 1:2]

  sumD <- sum(diag((fj.mat - X) %*% t(fj.mat - X)))
  eps <- 0.001
  count <- 0
  finish <- F

  while (!finish && count < max.iter) {
    for (j in 1:p) fj.mat[, j] <- stats::fitted(stats::loess(X[, j] ~ lambda.j[,
                                                                 1] + lambda.j[, 2], span = alpha))
    PP <- PP2 <- matrix(0, nrow = n, ncol = p)
    distances <- as.matrix(stats::dist(fj.mat))
    max_dist <- max(distances)
    distances[distances == 0] <- max_dist
    for (i in 1:n) {
      fA <- c() ; fB <- c() ; c.vec <- c() ; d.vec <- c()
      for (k in 1:n) {
        close_flam <- sort(distances[k, ])
        flamA <- as.numeric(names(close_flam)[1])
        flamB <- as.numeric(names(close_flam)[2])
        b1 <- fj.mat[flamA, ] - fj.mat[k, ]
        b2 <- fj.mat[flamB, ] - fj.mat[k, ]
        a <- X[i, ] - fj.mat[k, ]
        B = cbind(b1, b2)
        c1 = as.numeric((t(a) %*% b1 %*% (t(b2) %*% b2) -
                           t(a) %*% b2 %*% (t(b1) %*% b2))/det(t(B) %*%
                                                                 B))
        d1 = as.numeric((t(a) %*% b2 %*% (t(b1) %*% b1) -
                           t(a) %*% b1 %*% (t(b2) %*% b1))/det(t(B) %*%
                                                                 B))
        if (is.na(c1))
          c1 = 0
        if (c1 < 0)
          c1 = 0
        if (c1 > 1)
          c1 = 1
        if (is.na(d1))
          d1 = 0
        if (d1 < 0)
          d1 = 0
        if (d1 > 1)
          d1 = 1
        P <- c1 * b1 + d1 * b2
        PP[k, ] <- P + fj.mat[k, ]
        fA[k] = flamA ; fB[k] <- flamB ; c.vec[k] = c1 ; d.vec[k] = d1
      }
      dik <- as.matrix(stats::dist(rbind(X[i, ], PP)))[1, -1]
      small.d <- which.min(dik)
      PP2[i, ] <- PP[small.d, ]
    }
    fj.mat <- PP2

    SVD2 <- svd(fj.mat)
    lambda.j <- fj.mat %*% SVD2$v[,1:2]

    sumD.new <- sum(diag((fj.mat - X) %*% t(fj.mat - X)))
    eps1 <- abs(sumD.new - sumD)/sumD
    if (eps1 < eps)
      finish <- T
    sumD <- sumD.new
    count <- count + 1
    if(print_iterations==TRUE) print(c(round(count,0), eps1, sumD))
  }

  f.grid <- list(len = p)
  nn = N
  for (j in 1:p) {
    out <- akima::interp(x = lambda.j[, 1], y = lambda.j[, 2], z = fj.mat[,j],
                  xo = (seq(min(lambda.j[, 1]), max(lambda.j[,1]), length = nn)),
                  yo = (seq(min(lambda.j[, 2]), max(lambda.j[,2]), length = nn)),
                  linear = T, extrap = F,duplicate = "median")$z
    out[out > max(fj.mat[, j])] <- max(fj.mat[, j])
    out[out < min(fj.mat[, j])] <- min(fj.mat[, j])
    f.grid[[j]] <- out
  }


  f.grid2 = matrix(NA, nrow = nn^2, ncol = p)
  f.grid3 = matrix(NA, nrow = nn^2, ncol = p)
  index = 0
  for (i1 in 1:nn) for (i2 in 1:nn) {
    index = index + 1
    for (j in 1:p) {
      f.grid2[index, j] = f.grid[[j]][i1, i2]
      f.grid3[index, j] = f.grid[[j]][i2, i1]
    }
  }

  D.mat <- matrix(NA, nrow = (nn - 1)^2, ncol = nn)
  I1 = c()
  I2 = c()
  block = c()
  for (i1 in 2:nn) {
    for (i2 in 2:nn) {
      l11 <- i1 - 1
      l12 <- i1
      l21 <- i2 - 1
      l22 <- i2
      b1 <- b2 <- zero.point <- rep(NA, p)
      A <- X
      for (j in 1:p) {
        b1[j] <- f.grid[[j]][l12, l21] - f.grid[[j]][l11,
                                                     l21]
        b2[j] <- f.grid[[j]][l11, l22] - f.grid[[j]][l11,
                                                     l21]
        zero.point[j] <- f.grid[[j]][l11, l21]
        A[, j] <- X[, j] - f.grid[[j]][l11, l21]
      }
      mat <- rbind(0, b1, b1 + b2, b2)
      mat <- scale(mat, center = -zero.point, scale = F)
    }
  }
  lambda.grid <- f.grid2 %*% SVD2$v[,1:2]


  #rgl::points3d(f0.mat,col="red")
  rgl::points3d(fj.mat,col="seagreen")
  rgl::lines3d(f.grid2,col="lightgrey")
  rgl::lines3d(f.grid3,col="lightgrey")

  plot(lambda.j, col="seagreen",xlab = expression(lambda[1]),ylab=expression(lambda[2]),
       pch=16,yaxt="n",xaxt="n",asp=1)
  return(list(fj.mat = fj.mat,lambda.j = lambda.j))
}



