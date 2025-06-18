#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rmath.h>  // needed for some random number generators


/* -------------------------------------- */
/* C code from R to be passed to FORTRAN  */
/* -------------------------------------- */
/* Random seed related */
void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }

/* Random number generators */
double F77_SUB(unifrnd)(void) {
	return unif_rand(); } /* Uniform */
double F77_SUB(betarnd)(double *aa, double *bb){
	return rbeta(*aa,*bb); } /* Beta */
double F77_SUB(gammarnd)(double *a, double *scale){
	return rgamma(*a,*scale);} /* Gamma */

/* Density functions */
double F77_SUB(betadens)(double *x, double *a, double *b, int *give_log){
	return dbeta(*x, *a, *b, *give_log);} /* Beta */
double F77_SUB(gammadens)(double *x, double *shape, double *scale, int *give_log){
	return dgamma(*x, *shape, *scale, *give_log);} /* Gamma */


/* --------------------------- */
/* .Fortran calls */
/* --------------------------- */
extern void F77_NAME(dltest)(int * n, double * y, int * B, int * p, double * vals);

extern void F77_NAME(linkr)(int * link, double * par, int * n, int * ind, double * x, double * gx, double * dlink);

/* Simulating time series */
extern void F77_NAME(simbtsr)(int * code, int * length, int * order, double* ts, double * xreg1, double * xreg2, double * tstart, double * xstart, int * link, double * lconfig, int * np, double * pdist, int * xregar, double * alpha,  double * beta, double * phi, double * theta, double * d, int * rev);

extern void F77_NAME(simbarc)(int * length, int * order, double * ts, double * xreg, double * ystart, double * xstart, int * link, double * lconfig, int * map, int * xregar, double * alpha, double * beta, double * phi, double * theta, double * u0, int * rev);

/* Extracting components */
extern void F77_NAME(extractbtsr)(int * code, int * length, int * order, double * ts, double * xreg1, double * xreg2, double * tstart, double * xstart, double * xnew1, double * xnew2, double * forecast, int * link, double * lconfig, int * npar, double * par, int * xregar, int * nfix, double * alpha, int * flagsb, double * beta, int * flagsar, double * phi, int * flagsma, double * theta, double * d, int * np, double * pdist, int * extras, double * sll, double * U, double * K, int * nd, double * Dg, double * T, double * E, double * h, int * ierr);

extern void F77_NAME(extractbarc)(int * length, int * order, double * ts, double * xreg,  double * ystart, double * xstart, double * xnew, double * forecast, int * link, double * lconfig, int * map, int * npar, double * par, int * xregar, int * nfix, double * alpha, int * flagsb, double * beta, int * flagsar, double * phi, int * flagst, double * theta, double * u0, int * extras, double * sll, double * U, double * K, int * ierr);

extern void F77_NAME(gradient)(int * n, int * npar, int * nd, double * D, double * T, double * h, double * grad);

/* Fitting parameters */
extern void F77_NAME(optimbtsr)(int * method, int * code, int * length, int * order, double * ts, double * xreg1, double * xreg2, double * tstart, double * xstart, double * xnew1, double * xnew2, double * forecast, int * link, double * lconfig, int * npar, double * par, int * nbd, double * bounds, int * xregar, int * nfix, double * alpha, int * flagsb, double * beta, int * flagsar, double * phi, int * flagsma, double * theta, double * d, int * np, double * pdist, int * extras, double * sll, double * U, double * K, int * nd, double * Dg, double * T, double * E, double * h, int * cf1, int * nc2, double * cf2, int * neval, int * conv);

extern void F77_NAME(optimbarc)(int * method, int * length, int * order, double * ts, double * xreg, double * ystart, double * xstart, double * xnew, double * forecast, int * link, double * lconfig, int * map, int * npar, double * par, int * nbd, double * bounds, int * xregar, int * nfix, double * alpha, int * flagsb, double * beta, int * flagsar, double * phi, int * flagsma, double * theta, double * u0, int * extras, double * sll, double * U, double * K, int * cf1, int * nc2, double * cf2, int * neval, int * conv);

/* Prediction */
extern void F77_NAME(predictbtsr)(int * length, int * order, double * ts, double * xreg1, double * xreg2, double * xnew1, double * xnew2, double * forecast, int * link, double * lconfig, int * npar, double * par, int * xregar, int * nfix, double * alpha, int * flagsb, double * fvbeta, int * flagsar, double * fvar, int * flagsma, double * fvma, double * d);

extern void F77_NAME(predictbarc)(int * length, int * order, double * ts, double * xreg, double * xnew, double * forecast, int * link, double * lconfig, int * map, int * npar, double * par, int * xregar,  int * nfix, double * alpha, int * flagsb, double * fvbeta, int * flagsar, double * fvar, int * flagst, double * fvT, double * u0);

/* --------------------------- */
/* end of .Fortran calls */
/* --------------------------- */

static const R_FortranMethodDef FortranEntries[] = {
    {"dltest",           (DL_FUNC) &F77_NAME(dltest),                5},
    {"extractbtsr",      (DL_FUNC) &F77_NAME(extractbtsr),          37},
    {"extractbarc",      (DL_FUNC) &F77_NAME(extractbarc),          28},
		{"gradient",         (DL_FUNC) &F77_NAME(gradient),              7},
		{"linkr",            (DL_FUNC) &F77_NAME(linkr),                 7},
		{"optimbtsr",        (DL_FUNC) &F77_NAME(optimbtsr),            44},
		{"optimbarc",        (DL_FUNC) &F77_NAME(optimbarc),            35},
		{"predictbtsr",      (DL_FUNC) &F77_NAME(predictbtsr),          22},
		{"predictbarc",      (DL_FUNC) &F77_NAME(predictbarc),          21},
		{"simbtsr",          (DL_FUNC) &F77_NAME(simbtsr),              19},
		{"simbarc",          (DL_FUNC) &F77_NAME(simbarc),              16},
    {NULL, NULL, 0}
};

void R_init_BTSR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


