## Resubmission addressing the September 22, 2026 review

Thank you for reviewing toxplot. All three requested changes are addressed:

1. DESCRIPTION now cites Filer et al. (2017)
   <doi:10.1093/bioinformatics/btw680> for the ToxCast fitting pipeline and
   Wang et al. (2018) <doi:10.1021/acs.est.7b06145> for toxicity-adjusted ranking.
2. save_plot_pdf.Rd now contains a value section documenting the invisible
   NULL return and the side effect of writing one PDF page per plot. The
   roxygen source was updated and the documentation regenerated. All eight
   exported functions have value documentation.
3. Progress output in R/functions_tcpl_based.R and R/utilities.R now uses
   message() and can be suppressed with suppressMessages(). The print()
   calls used to draw ggplot objects on the PDF device remain necessary.
   Regression tests verify silent operation under suppressMessages(), the
   invisible NULL return, and PDF device cleanup.

## Resubmission of archived package

This is toxplot 0.1.2, maintained by the original author Jun Wang
(njekin@gmail.com).

The package was archived on 2018-12-10 because examples failed while
converting tcplFit results to a data.frame:
"arguments imply differing number of rows: 1, 18, 0".

The conversion now retains vector and NULL outputs as list columns while
keeping scalar model statistics in a single row. Regression tests exercise
the full failing example, inactive and insufficient-concentration curves,
ranking, both plotting functions, and PDF device cleanup. The unused tidyr
import has been removed, ggplot2 line parameters updated, and the vignette
now builds against the installed package.

## Verified local environment (September 25, 2026)

* macOS, Apple silicon, R 4.2.3
* tcpl 3.3.1, ggplot2 4.0.3, ggthemes 6.0.0, dplyr 1.2.1
* R CMD build, including the HTML vignette
* R CMD check --as-cran (including PDF manual)

Local check result: 0 ERRORs, 0 WARNINGs, 2 NOTEs.

1. CRAN incoming feasibility identifies a new submission of an archived
   package. This is the intended reinstatement request. The same incoming
   NOTE reports HTTP 403 (Forbidden) responses while validating both DOIs.
   The citations have been independently verified against the published
   articles:
   https://pmc.ncbi.nlm.nih.gov/articles/PMC6697091/
   https://academic.oup.com/bioinformatics/article/33/4/618/2617576
2. The local system could not verify the current time.

All 82 regression assertions pass, including the new output-suppression
and return-value checks. Examples and rebuilding vignette outputs pass.
The PDF reference manual compiled successfully with TeX Live 2026.

## Cross-platform checks of the previous submission

GitHub Actions run 34739796890 completed successfully on 2026-09-13:
https://github.com/njekin/ToxPlot-R-Package/actions/runs/34739796890

These checks predate the present review fixes; they have not been rerun
for this resubmission. All five jobs reported Status: OK:

* Windows, R 4.6.1
* macOS, R 4.6.1
* Ubuntu, R 4.6.1
* Ubuntu, R 4.5.3
* Ubuntu, R-devel (2026-09-11 r90528)

The checked change was merged in pull request #1, merge commit
f3f6df6d8596731e5ff2b5f62af9b9db74e701fa.
These CI jobs used --no-manual. The full local check subsequently verified
the PDF manual successfully without --no-manual.
