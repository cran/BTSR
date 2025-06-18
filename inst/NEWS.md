### 1.0.0 (June-2025)

* Minor changes
  - Renamed several modules and moved functions and subroutines for better file
    organization.
  - Removed `00_specfun` module; migrated `trigamma` and `psi` to `00_lib_utils`. 
    Other functions in `00_specfun` are no longer necessary. The COPYRIGHTS and 
    DESCRIPTION files have been updated to reflect the changes. 
  - Consolidated `03_barfima`, `03_karfima`, `03_garfima`, and `03_uwarfima` into 
    a single module named `05_generic`.
  - Simplified the internal structure of the `link.btsr` function.  

* Major Changes
  - Rewrote and renamed `01_RNG` to `01_distrib.f90`, now using R's native 
    `rand_unif`, `rbeta`, `dbeta`, `rgamma`, and `dgamma`. The COPYRIGHTS and 
    DESCRIPTION files have been updated to reflect these changes.
  - Modified subroutines in `02_base` to support time-varying `nu`.
  - Added input validation and variable conversion functions for backward 
    compatibility.
  - The `btsr.sim`, `btsr.extract`, and `btsr.fit` functions have been updated 
    to support the new models. Regression models can now be easily fitted using 
    these functions.
  - Fixed parameter initialization typo:  
    Old: `sigma2 <- sum(er^2)/((n1 - k) * (dlink)^2)`  
    New: `sigma2 <- sum(er^2/((n1 - k) * (dlink)^2))`
  - Changed the initialization of `xreg` from zero to the mean of the first 
    `p` observations.
  - Added a robust covariance calculation method that falls back to the outer 
    product of gradients if Hessian calculation fails.


### 0.1.5 (2023-09-23)

* Minor Revision

* Added the "_PACKAGE" special sentinel to fix the problem related to the fact 
  that @docType package no longer automatically adds a -package alias.

* Temporary disabled all subroutines that use integer(kind=8) (file 01_RNG.f90).

* Fixed the ``not yet implemented: derived type components with non default 
  lower bounds'' message due to mag01 component in 01_RNG.f90 file and changed 
  subroutines that use this component accordingly.


### 0.1.4 (2023-01-19)

* Minor fix

* We fixed some typos in src/BTSR_init.c and added the information 
  about this file to src/Makevars (this was missing in the last 
  version)
  
* Internal subroutines to simulate, extract and fit models were 
  modified, but this caused no change in the results from previous 
  versions. This change was done since using a construction such as
  'foo <- .check.model(model[1], "extract")' and then
  '.Fortran(foo,...)' gives an error during the registration process.
  Therefore, we created auxiliary functions for each model.
  
* We have fixed a typo in the FORTRAN subroutine used to predict BARC
  models. This should fix the memory allocation problem found in the 
  last version.

### 0.1.3 (2023-01-10)

* Minor fix

* We have added the file src/BTSR_init.c 


### 0.1.2 (2022-11-16)

* Minor fix

* 00_lb.f90 was merged with 00_lbfgsb.f90

* We have added more information to src/Makevars in order to 
  setup dependencies for parallel make.  

### 0.1.1 (2022-11-10)

* Minor fix

* We have added explicit dependencies to src/Makevars 
  to tell make the constraints on the order of compilation.  

### 0.1.0 

* First version

* Submission date:   2022-11-07
