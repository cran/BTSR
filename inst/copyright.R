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
#################################################################################

##--------------------------------------------------------------------------------
src/00_lbfgsb.f90

Contains the subroutines
	dtrsl, dpofa, ddot
which are part of LINPACK, with authors J.J. Dongarra, Cleve Moler and G.W. Stewart.
They were taken from the Netlib archive now at www.netlib.org and do not 
clearly state their copyright status.

Contains the subroutines of l-bfgs-b algorithm Written by 
Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis Morales and Peihuang Lu-Chen.
The original code is now available at 
http://users.iems.northwestern.edu/~nocedal/lbfgsb.html

Condition for Use: This software is freely available, but we expect that 
all publications describing  work using this software, or all commercial 
products using it, quote at least one of the references given below. 
This software is released under the "New BSD License" (aka "Modified BSD License" 
or "3-clause license").
 
References
 
R. H. Byrd, P. Lu and J. Nocedal. A Limited Memory Algorithm for Bound Constrained
Optimization, (1995), SIAM Journal on Scientific and Statistical Computing , 16, 5,
pp. 1190-1208.
 
C. Zhu, R. H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN
routines for large scale bound constrained optimization (1997), ACM Transactions on
Mathematical Software, Vol 23, Num. 4, pp. 550 - 560.

J.L. Morales and J. Nocedal. L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN
routines for large scale bound constrained optimization (2011), to appear in ACM
Transactions on Mathematical Software.
##-----------------------------------------------------------------------------

##--------------------------------------------------------------------------------
src/00_main.f90

The subroutines 
	xtransform, xtransformstart
were based on the Matlab subroutine 'fminsearchbnd' available at
https://www.mathworks.com/matlabcentral/fileexchange/8277-
      fminsearchbnd-fminsearchcon?focused=5216898&tab=function
##--------------------------------------------------------------------------------

##--------------------------------------------------------------------------------
src/00_specfun.f90

The function
	trigama
was taken from
https://people.sc.fsu.edu/~jburkardt/f_src/asa121/asa121.html
The FORTRAN90 version was written by John Burkardt.
This code is distributed under the GNU LGPL license.

The function
	lngamma
was taken from https://jblevins.org/mirror/amiller/lanczos.f90
This function was written by Alan Miller. 
https://jblevins.org/mirror/amiller/ is an archived copy of the 
Fortran source code repository of Alan Miller previously located at 
http://users.bigpond.net.au/amiller/. It is hosted by Jason Blevins 
with permission. All code written by Alan Miller is released into 
the public domain. Code written by other authors or from other sources 
(e.g., academic journals) may be subject to other restrictions.

The functions
	alnrel, algdiv, gsumln, bcorr, betaln, rlog1, 
	gamln1, gamln, gam1, brcomp, ipmpar, dpmpar, psi
were taken from https://jblevins.org/mirror/amiller/specfunc.zip
These are special functions from the NSWC library. The file 
specfunc.zip was compiled by Alan Miller.
At https://github.com/jacobwilliams/nswc one reads: 
"The NSWC Mathematics Subroutine Library is a collection of Fortran 77 
routines specializing in numerical mathematics collected and developed 
by the U.S. Naval Surface Warfare Center. This software is made available, 
without cost, to the general scientific community. The 1993 edition is an 
update of the 1990 edition. NSWC has made every effort to include only reliable,
transportable, reasonably efficient and easy to use code in this library. They 
have thoroughly tested all the routines on a variety of machines ranging
from supercomputers to PC's."
##--------------------------------------------------------------------------------


##--------------------------------------------------------------------------------
src/01_Nelder.f90

Contains the subroutine 
	minim	
programmed by D.E. Shaw, with amendments by R.W.M. Wedderburn.
Further amended by Alan Miller. Further modified by Taiane Schaedler Prass.
	
The original code available at https://jblevins.org/mirror/amiller/minim.f90
##--------------------------------------------------------------------------------

##--------------------------------------------------------------------------------
src/01_RNG.f90

rng_seed_Blevins, rng_uniform_Blevins
were taken from https://jblevins.org/log/openmp
They were written by Jason Blevins.  	
	
rng_uniform_wh
was taken from https://people.math.sc.edu/Burkardt/f_src/asa183/asa183.f90
The FORTRAN90 version was written by John Burkardt. 
This code is distributed under the GNU LGPL license.

rng_seed_sgrnd, rng_uniform_Mersenne
were taken from https://jblevins.org/mirror/amiller/mt19937a.f90
Fortran translation by Hiroshi Takano. Code converted to Fortran90 by Alan Miller. 

rng_uniform_kiss32
is public domain code. It was taken from http://www.fortran.com/kiss.f90

rng_uniform_kiss64
was taken from http://lgge.osug.fr/meom/pages-perso/brankart/Outils/mod_kiss.f90
This module was written by Jean-Michel Brankart.

rng_array, rng_seed_rnstrt
were taken from https://jblevins.org/mirror/amiller/rand3.f90
FORTRAN 77 version by Steve Kifowit with modifications by Alan Miller 
based upon the code written by Knuth. The code was converted to FOTRAN90 
by Alan Miller. The FORTRAN77 code written by Donald E. Knuth is available at
https://www-cs-faculty.stanford.edu/~knuth/programs/frng.f

rng_seed_lfsr258, rng_uniform_Le
were taken from https://jblevins.org/mirror/amiller/lfsr258.f90
Fortran version by Alan Miller.

random_beta, standard_qnorm
were taken from https://jblevins.org/mirror/amiller/random.f90
https://jblevins.org/mirror/amiller/as241.f90
They were all written by Alan Miller. All code written by Alan Miller is 
released into the public domain.

dgamma_default, dpois_raw, bd0, stirlerr
are based on the code in "dgamma.c", "dpois.c", "bd0.c", "stirlerr.c", 
found in "R-4.1.0/src/nmath/". These codes were written by Catherine Loader.
