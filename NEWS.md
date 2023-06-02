# SAPP 1.0.9-1

* Fixed an issue with the S3 method that was reported as NOTE in the pre-test.


# SAPP 1.0.9

* Removed C wrapper functions and registered entry points for the routines accessed by the `.Fortran` interface to call Fortran subroutines.

* Added a ‘`NEWS.md`’ file.


# SAPP 1.0.8

* Fixed the LTO (Link-Time Optimization) problems reported in CRAN package check results.

* Fixed the following problems reported in CRAN package check results.

	checking compiled code ... NOTE
 
	    File `SAPP/libs/SAPP.so':
 
		Found no calls to: `R_registerRoutines', `R_useDynamicSymbols'

* Fixed Fortran code according to the warning message when using gfortran with `-fpic  -g -O2 -Wall -mtune=native -c`.

* Revised the guide so that it can be referenced in `vignette("SAPP")`.


# SAPP 1.0.7

* Fixed Fortran code according to the warning message when using gfortran with `-Wall -pedantic`.


# SAPP 1.0.6

* Fixed a bug in `etasap()`.
  The order of the elements in the output value `param` is wrong.
  (Reported by Amato Kasahara.)


# SAPP 1.0.5

* Fixed a bug in `etasap()`.
  If approx = 0, the minimization function is not called.
  (Reported by Amato Kasahara.)


# SAPP 1.0.4

* Corrected ‘`inst/doc/index.html`‘ according to the error message from the W3C Markup Validator Service.


# SAPP 1.0.3

* Fixed the dimensions for array bounds (1) and (*) in eptrenf.f, etasimf.f, comsub.f, linsimf.f, aftpoi.f, pgraphf.f and simbvhf.f.

* Added new argument `ier` to the subroutines `linsimf` (linsimf.f) and `simbvhf` (simbvhf.f).


# SAPP 1.0.2

* Fixed a bug in `momori()`.
  Fixed the dimension of array ti in subroutine `momorif` (aftpoi.f). (Reported by Brian Ripley.)

* Fixed a bug in `pgraph()`.
  Fixed the dimension of array xx in subroutine `palmpr` (pgraphf.f).

* Added C wrapper functions for calling Fortran subroutines.

* Changed output to `tmpfile` from R functions (eptren(), linlin(), etasap() and momori()) instead of from Fortran subroutine.
 
  Added new argument `nlmax` to `eptren()`, `linlin()`, `etasap()` and `momori()` functions.


# SAPP 1.0.1

* Deleted an extra line with the unterminated string at the end of the file ‘`data/res2003JUL26.R`‘.
* Deleted the extra line of unterminated string at the end of the '`data/res2003JUL26.R` file.
  (Reported by William Dunlap and Duncan Murdoch.)

* Fixed a bug in `eptren()`. 
  Fixed the dimension of array amg in subroutine `eptrenf` (eptrenf.f). 

* Added new argument `kmax` to the subroutines `linlinf`, `dav`, `linear`, `davidn`, `funct` (linlinf.f) and `comfac` (comsub.f).
  In these subroutines, the dimension of the array lf(kmax,kmax) is declared as a dynamic array.