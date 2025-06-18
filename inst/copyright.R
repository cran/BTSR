#################################################################################
## R package BTSR is copyright Taiane Schaedler Prass and Guilherme Pumi,
## with the exceptions described in the sequel.
##
## This file is part of the R package BTSR.
##
## The R package BTSR is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package BTSR is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################

##------------------------------------------------------------------------------
src/00_lbfgsb.f90

Contains the subroutines dtrsl, dpofa, ddot which are part of LINPACK, with
authors J.J. Dongarra, Cleve Moler and G.W. Stewart. They were taken from the
Netlib archive now at www.netlib.org and do not clearly state their copyright
status.

Contains the subroutines of l-bfgs-b algorithm Written by
Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis Morales and Peihuang Lu-Chen.
The original code is now available at
http://users.iems.northwestern.edu/~nocedal/lbfgsb.html

Condition for Use: This software is freely available, but we expect that all publications
describing  work using this software, or all commercial products using it, quote
at least one of the references given below. This software is released under the
"New BSD License" (aka "Modified BSD License" or "3-clause license").

References

R. H. Byrd, P. Lu and J. Nocedal. A Limited Memory Algorithm for Bound Constrained
Optimization, (1995), SIAM Journal on Scientific and Statistical Computing , 16, 5,
pp. 1190-1208.

C. Zhu, R. H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN
routines for large scale bound constrained optimization (1997), ACM Transactions
on Mathematical Software, Vol 23, Num. 4, pp. 550 - 560.

J.L. Morales and J. Nocedal. L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN
routines for large scale bound constrained optimization (2011), to appear in ACM
Transactions on Mathematical Software.
##-----------------------------------------------------------------------------

##------------------------------------------------------------------------------
src/00_main.f90

The subroutines xtransform, xtransformstart were based on the Matlab subroutine
'fminsearchbnd' available at
https://www.mathworks.com/matlabcentral/fileexchange/8277-
      fminsearchbnd-fminsearchcon?focused=5216898&tab=function
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
src/00_lib_utils.f90

The function trigamma was taken from
https://people.sc.fsu.edu/~jburkardt/f_src/asa121/asa121.html
The FORTRAN90 version was written by John Burkardt.
This code is distributed under the GNU LGPL license.

The functions psi was taken from https://jblevins.org/mirror/amiller/specfunc.zip
These are special functions from the NSWC library. The file specfunc.zip was
compiled by Alan Miller. At https://github.com/jacobwilliams/nswc one reads:
"The NSWC Mathematics Subroutine Library is a collection of Fortran 77
routines specializing in numerical mathematics collected and developed
by the U.S. Naval Surface Warfare Center. This software is made available,
without cost, to the general scientific community. The 1993 edition is an
update of the 1990 edition. NSWC has made every effort to include only reliable,
transportable, reasonably efficient and easy to use code in this library. They
have thoroughly tested all the routines on a variety of machines ranging
from supercomputers to PC's."
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
src/01_Nelder.f90

Contains the subroutine minim	programmed by D.E. Shaw, with amendments by
R.W.M. Wedderburn. Further amended by Alan Miller. Further modified by
Taiane Schaedler Prass. The original code available at
https://jblevins.org/mirror/amiller/minim.f90
##------------------------------------------------------------------------------
