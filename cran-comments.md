## Test environments
* local OS X install, R 4.0.0
* local Windows 10 install, R 4.0.0


## R CMD check results
There were no ERRORs, WARNINGs or NOTES.

Below is the build record:

==> devtools::check(document = FALSE)

Setting env vars --------------------------------------------------------------
CFLAGS  : -Wall -pedantic
CXXFLAGS: -Wall -pedantic
Building toxplot --------------------------------------------------------------
'/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file --no-environ  \
  --no-save --no-restore --quiet CMD build '/Users/junwang/Dropbox/EPA  \
  works/toxplot' --no-resave-data --no-manual 

* checking for file ‘/Users/junwang/Dropbox/EPA works/toxplot/DESCRIPTION’ ... OK
* preparing ‘toxplot’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files
* checking for empty or unneeded directories
* looking to see if a ‘data/datalist’ file should be added
* building ‘toxplot_0.1.0.tar.gz’

Setting env vars --------------------------------------------------------------
_R_CHECK_CRAN_INCOMING_ : FALSE
_R_CHECK_FORCE_SUGGESTS_: FALSE
Checking toxplot --------------------------------------------------------------
'/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file --no-environ  \
  --no-save --no-restore --quiet CMD check  \
  '/var/folders/9x/8jptjf_d5kxgjnj9qxrkzp200000gn/T//RtmpMyJWKG/toxplot_0.1.0.tar.gz'  \
  --as-cran --timings --no-manual 

* using log directory ‘/Users/junwang/Dropbox/EPA works/toxplot.Rcheck’
* using R version 3.4.0 (2017-04-21)
* using platform: x86_64-apple-darwin15.6.0 (64-bit)
* using session charset: UTF-8
* using options ‘--no-manual --as-cran’
* checking for file ‘toxplot/DESCRIPTION’ ... OK
* checking extension type ... Package
* this is package ‘toxplot’ version ‘0.1.0’
* package encoding: UTF-8
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking whether package ‘toxplot’ can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking ‘build’ directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of ‘data’ directory ... OK
* checking data for non-ASCII characters ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking installed files from ‘inst/doc’ ... OK
* checking files in ‘vignettes’ ... OK
* checking examples ... NONE
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in ‘inst/doc’ ... OK
* checking re-building of vignette outputs ... OK
* DONE
Status: OK



R CMD check results
0 errors | 0 warnings | 0 notes

R CMD check succeeded

