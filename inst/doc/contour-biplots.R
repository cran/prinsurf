## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6.5, fig.height = 5.2, fig.align = "center")
set.seed(1)

## ----fit-iris-----------------------------------------------------------------
library(prinsurf)

fit <- prinsurf(iris[, 1:4], scale = TRUE)   # standardise heterogeneous units
fit

## ----plot-bare----------------------------------------------------------------
plot(fit, group = iris$Species)

## ----plot-iris----------------------------------------------------------------
plot(fit, vars = colnames(iris)[1:4], group = iris$Species)

## ----plot-subset--------------------------------------------------------------
plot(fit, vars = c("Petal.Length", "Petal.Width"), group = iris$Species)

## ----plot-titles--------------------------------------------------------------
plot(fit, vars = c("Petal.Length", "Petal.Width"), group = iris$Species,
     main = c("Petal length (cm)", "Petal width (cm)"),
     outer_main = "Iris contour biplot")

## ----predict-iris-------------------------------------------------------------
phat <- predict(fit)
head(round(phat, 2))
head(iris[, 1:4])

## ----contour-error------------------------------------------------------------
contour_predictive_error(fit)

## ----sample-pred--------------------------------------------------------------
pred <- predictivity(fit)
summary(pred)
attr(pred, "overall")

## PCA biplot (rank-2) sample predictivity on the same standardised data, for comparison
Z   <- scale(as.matrix(iris[, 1:4]))
V2  <- svd(Z)$v[, 1:2]; Zhat <- Z %*% V2 %*% t(V2)
pca <- mean(1 - rowSums((Z - Zhat)^2) / rowSums(Z^2))

c(principal_surface = round(attr(pred, "overall"), 3), pca_biplot = round(pca, 3))

## ----halfsphere---------------------------------------------------------------
## upper cap of a sphere: x = height, y, z = horizontals
n <- 200; u <- 2 * runif(n) - 1; th <- 2 * pi * runif(n) - pi
S <- cbind(x = u, y = sin(th) * sqrt(1 - u^2), z = cos(th) * sqrt(1 - u^2))
S <- S[S[, "x"] > -0.4, ]                 # keep the cap
S <- sweep(S, 2, colMeans(S))             # centre
sph <- prinsurf(S, max.iter = 8)
plot(sph, vars = colnames(S))

## ----pca-biplot-sphere, eval = requireNamespace("biplotEZ", quietly = TRUE)----
biplotEZ::biplot(data = S) |> biplotEZ::PCA() |> plot()

## ----sphere-predictivity------------------------------------------------------
ps  <- mean(predictivity(sph))

Zc  <- scale(S, scale = FALSE)
V2  <- svd(Zc)$v[, 1:2]; Zhat <- Zc %*% V2 %*% t(V2)
pca <- mean(1 - rowSums((Zc - Zhat)^2) / rowSums(Zc^2))

c(principal_surface = round(ps, 3), pca_biplot = round(pca, 3))

