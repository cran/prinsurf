## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.height = 6, fig.width = 7,
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(prinsurf)
library(rgl)

## -----------------------------------------------------------------------------
data2 <- function(n,side)
{
  e1 <- runif(n,0,0.5)
  x <- 2*runif(n) - 1
  theta <- 2*pi*runif(n) - pi
  y <- sin(theta)*sqrt(1-x^3)
  z <- sin(theta)*sqrt(x^2)
  X = matrix(c(x+e1,y+e1,z+e1),nrow=n,ncol=3,byrow=F)
  return(X)
}
X <-  data2(100,1)
points3d(X,col="blue")

## ----echo=FALSE---------------------------------------------------------------
my.plot <- scene3d()
rglwidget(my.plot)

## ----example------------------------------------------------------------------
surface <- principal.surface(X,max.iter = 5,print_iterations = TRUE)

## ----echo=FALSE---------------------------------------------------------------
my.plot <- scene3d()
rglwidget(my.plot)

