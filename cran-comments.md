## Resubmission
This is a resubmission. In this version I have addressed the following notes:

* installed size is 8.4Mb\
It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is then above the 1MB threshold. 
However, it seems that this NOTE only appears under LINUX, but not under Windows or OSX. 
My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp.
Indeed, some functions of the penppml package have been written in C++ using Rcpp. They are needed to perform faster computation of summary statistics, OLS and ridge regressions.
Without the speed up gained from those C++ functions, this package would become impractical.

* checking dependencies in R code ... All declared Imports should be used.\
I excluded the packages not used.

* Note: found 2793 marked UTF-8 strings\
I changed strings to ASCII-characters.

* Additional issues: noLD\
One of the automated tests failed. I was not able to reproduce it with rhub, 
with the platform Debian Linux, R-devel, GCC, no long double (and under any other platform). 
The reason for the error seems dependent on a value drawn in the test. However,
the subset of the parameter space in which this occurs is not relevant for our relevant.
To address this problem, I set a seed in the relevant test.

* Corrected standardization in fastridge(), in R/utils.R

* Replaced lfe::demeanlist() with collapse::fhdwithin()\
lfe-package is not maintained anymore, therefore I switched to the fhdwithin()-function from the collapse package which performs the same operations. I then updated DESCRIPTION.

* checking dependencies in R code ... fixest: All declared Imports should be used.\
collapse-package requires loading fixest beforehand.\
When not including fixest in Imports, additional errors arise.

## Test environments
* Windows 10, 64 bit, R 4.1.2
* Debian Linux, R-devel, GCC, no long double

## R CMD check results

* Windows:
0 errors | 0 warnings | 0 notes

* Linux:
0 errors | 0 warnings | 2 notes
