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

## Verified local environment

* macOS, Apple silicon, R 4.2.3
* tcpl 3.3.1, ggplot2 4.0.3, ggthemes 6.0.0, dplyr 1.2.1
* R CMD build, including the HTML vignette
* R CMD check --as-cran --no-manual

Local check result: 0 ERRORs, 0 WARNINGs, 2 NOTEs.

1. CRAN incoming feasibility identifies a new submission of an archived
   package. This is the intended reinstatement request.
2. The local system could not verify the current time.

Examples, tests and rebuilding vignette outputs pass.
The PDF reference manual was not checked (--no-manual).

## Still required before submission

This file is a preparation draft, not a claim of completed cross-platform
validation. Run the supplied GitHub Actions workflow on current R release,
R-devel and R-oldrel across Linux, macOS and Windows, and update this file
with the actual results. Check the PDF reference manual on a TeX-equipped
machine. Do not submit until any remaining significant check issues are
resolved or explained.
