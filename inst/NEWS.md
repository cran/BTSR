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

* Submission date: 	2022-11-07
