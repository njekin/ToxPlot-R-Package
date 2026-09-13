# toxplot 0.1.2

* Fix the CRAN archive-blocking example failure when `tcplFit()` returns
  vectors and NULL predictions: model tables now contain one row per fit,
  with list columns for concentration and prediction vectors.
* Correct the plate identifier extracted from assay data and use the exact
  primary assay name when checking the assay configuration.
* Use `linewidth` for ggplot2 lines and a proper margin object for legends.
* Close the PDF graphics device even when plot printing fails.
* Remove the unused tidyr import and development-only devtools suggestion.
* Build the vignette using the installed package and declare knitr as its builder.
* Add regression coverage for the complete archived failing example, inactive
  curves, insufficient concentrations, ranking, plotting and PDF cleanup.

# toxplot v0.1.1
------------------
## Bug fixes
- Fix crash when imported data contains NAs
- Removed vignette. Refer to github repo for usage details.

# toxplot v0.1.0
------------------
## Initial release
