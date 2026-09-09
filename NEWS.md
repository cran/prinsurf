# prinsurf 2.0

## New features

* Principal surface fits are now displayed as **contour biplots**. `plot()` draws
  each variable as the contour lines of its fitted surface coordinate function
  `f_j`, so a variable's value is read off its contours at a sample's
  position on the surface, in place of the linear axes of a conventional biplot.

* `plot()` gains arguments for tailoring the display: `vars` to select which
  variables are drawn, `group` to colour samples by a factor, `nlevels` and
  `col_contour` to control the contour axes, and `main` / `outer_main` to title
  individual panels and the overall figure.

* `predictivity()` reports **sample predictivity** — the proportion of each
  sample's squared length that is reconstructed by the fitted surface. The
  per-sample values are returned as a vector, with the mean over all samples in
  the `"overall"` attribute.

* `contour_predictive_error()` reports, for each variable, the root-mean-square
  difference between the value read from its contour lines and the sample's
  actual value. This measures how accurately a variable can be recovered by
  reading its contours; the mean over variables is in the `"overall"` attribute.

* `predict()` returns the values read from the contour axes at each sample's
  position, back-transformed to the variables' original scales.

* `fitted()` returns the fitted surface coordinates `f_j(lambda_i)`
  for every sample and variable, in the working (centred / scaled) units.

* `print()` summarises a fit: the number of samples and variables, the loess
  span, the number of iterations to convergence, and the variable names.

## Other changes

* The package is now released under the MIT licence (previously GPL-3).

* Added a package website at <https://raeesaganey91.github.io/prinsurf/> and a
  bug tracker at <https://github.com/RaeesaGaney91/prinsurf/issues>.

* The `Description` field now cites Hastie and Stuetzle (1989)
  <doi:10.1080/01621459.1989.10478797>.

* Added the vignette *Contour biplots*, working through a fit and its
  diagnostics end to end.


# prinsurf 1.0

* First CRAN release, fitting principal surfaces as the two-dimensional
  generalisation of the principal curves of Hastie and Stuetzle (1989).
