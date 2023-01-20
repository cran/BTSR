##----------------------------------------------------------
##    Generic functions: sim, extract and fit
##----------------------------------------------------------
##' @title
##' Generic functions to simulate, extract components and fit BTSR models
##'
##' @name btsr.functions
##' @order 1
##'
##' @description
##' These generic functions can be used to simulate, extract components
##' and fit any model of the class \code{btsr}. All functions are wrappers
##' for the corresponding function associated to the chosen model.
##' See \sQuote{The BTSR structure} and \sQuote{Common Arguments}.
##'
##' @details
##'
##' # The BTSR structure
##'
##' The general structure of the deterministic part of a BTSR model is
##'
##'  \deqn{g_1(\mu_t) = \alpha + X_t\beta +
##'  \sum_{j=1}^p \phi_j[g_2(y_{t-j}) - I_{xregar}X_{t-j}\beta] + h_t}
##'
##' where
##'  \itemize{
##'    \item \eqn{I_{xregar}} is 0, if \code{xreg} is not included in the AR part of the model and 1,
##' otherwise
##'
##'   \item the term \eqn{h_t} depends on the argument \code{model}:
##'    \itemize{
##'     \item for BARC models: \eqn{h_t =  h(T^{t-1}(u_0))}
##'     \item otherwise: \eqn{h_t =  \sum_{k = 1}^\infty c_k r_{t-k}}
##'    }
##'
##'   \item \eqn{g_1} and \eqn{g_2} are the links defined in \code{linkg}.
##'   Notice that \eqn{g_2} is only used in the AR part of the model and, typically,
##'   \eqn{g_1 = g_2}.
##'
##'   \item \eqn{r_t} depends on the \code{error.scale} adopted:
##'   \itemize{
##'     \item  if \code{error.scale = 0}: \eqn{r_t = y_t - \mu_t} (data scale)
##'     \item if \code{error.scale = 1}:  \eqn{r_t = g_1(y_t) - g_1(\mu_t)}
##'      (predictive scale)
##'   }
##'
##'    \item \eqn{c_k} are the coefficients of \eqn{(1-L)^d\theta(L)}.
##'    In particular, if \eqn{d = 0}, then \eqn{c_k = \theta_k}, for
##'    \eqn{k = 1, \dots, q}.
##' }
##'
##' # Common Arguments
##'
##' In what follows we describe some of the arguments that are
##' commom to all BTSR models. For more details on extra arguments,
##' see the corresponding function associated to the selected model.
##'
##' @md
NULL
#> NULL


##' @rdname btsr.functions
##' @order 2
##'
##' @details
##'
##' The function \code{btsr.sim} is used to generate random samples
##' from BTSR models. See \sQuote{The BTSR structure}.
##'
##' # Common Arguments
##'
##' ## Simulation Function
##'
##' Common arguments passed through \code{"..."} in \code{btsr.sim} are:
##'\itemize{
##'  \item \code{n} a strictly positive integer. The sample size of yt (after burn-in).
##'  Default for all models is 1.
##'
##'  \item \code{burn} a non-negative integer. length of "burn-in" period.
##'  Default for all models is 0.
##'
##'  \item \code{xreg} optionally, a vector or matrix of external regressors.
##'  For simulation purposes, the length of xreg must be \code{n+burn}.
##'  Default for all models is \code{NULL}
##'
##'  \item \code{coefs} a list with the coefficients of the model. Each model has
##'  its default. An empty list will result in an error. The arguments in this list
##'  are:
##'  \itemize{
##'   \item \code{alpha} optionally, A numeric value corresponding to the intercept.
##'    If the argument is missing, it will be treated as zero.
##'
##'   \item \code{beta} optionally, a vector of coefficients corresponding to the
##'   regressors in \code{xreg}. If \code{xreg} is provided but \code{beta} is
##'   missing in the \code{coefs} list, an error message is issued.
##'
##'   \item \code{phi} optionally, a vector of size \eqn{p}, corresponding to the
##'   autoregressive coefficients (including the ones that are zero), where \eqn{p}
##'   is the AR order.
##'
##'   \item \code{nu} the dispersion parameter. If missing, an error message is issued.
##'
##'   \item \code{rho, y.lower, y.upper, theta, d, u0} model specif arguments.
##'   See the documentation corresponding to each model.
##'  }
##'
##'  \item \code{y.start} optionally, a initial value for yt (to be used
##'  in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{g_2(y_t) = 0}, for \eqn{t < 1}.
##'
##'  \item \code{xreg.start} optionally, a vector of initial value for xreg
##'  (to be used in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{X_t = 0}, for \eqn{t < 1}. If \code{xregar = FALSE} this argument
##'  is ignored.
##'
##'  \item \code{xregar} logical; indicates if xreg is to be included in the
##'  AR part of the model.  See \sQuote{The BTSR structure}. Default is \code{TRUE}.
##'
##'  \item \code{error.scale} the scale for the error term. See also \sQuote{The BTSR structure}.
##'  Each model has its default.
##'
##'  \item \code{inf} the truncation point for infinite sums. Default is 1000.
##'  In practice, the Fortran subroutine uses \eqn{inf = q}, if \eqn{d = 0}.
##'  BARC models do not have this argument.
##'
##'  \item \code{linkg} character or a two character vector indicating which
##'  links must be used in the model.  See \sQuote{The BTSR structure}.
##'  If only one value is provided, the same link is used for \eqn{mu_t} and
##'  for \eqn{y_t} in the AR part of the model. Each model has its default.
##'
##'  \item \code{seed} optionally, an integer which gives the value of the fixed
##'  seed to be used by the random number generator. If missing, a random integer
##'  is chosen uniformly from 1,000 to 10,000.
##'
##'  \item \code{rngtype} optionally, an integer indicating which random number generator
##'  is to be used. Default is 2. The current options are:
##' \itemize{
##'     \item \code{0}: Jason Blevins algorithm. Available at <https://jblevins.org/log/openmp>
##'     \item \code{1}: Wichmann-Hill algorithm  (Wichmann and Hill, 1982).
##'     \item \code{2}: Mersenne Twister algorithm (Matsumoto and Nishimura, 1998).
##'     FORTRAN code adapted from <https://jblevins.org/mirror/amiller/mt19937.f90> and
##'     <https://jblevins.org/mirror/amiller/mt19937a.f90>
##'     \item \code{3}: Marsaglia-MultiCarry algorithm - kiss 32. Random number generator suggested
##'     by George Marsaglia in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
##'     in 1999.
##'     \item \code{4}: Marsaglia-MultiCarry algorithm - kiss 64. 	Based on the
##'      64-bit KISS (Keep It Simple Stupid) random number generator distributed by
##'      George Marsaglia in <https://groups.google.com/d/topic/comp.lang.fortran/qFv18ql_WlU>
##'     \item \code{5}: Knuth's 2002 algorithm (Knuth, 202). FORTRAN code adapted
##'     from <https://www-cs-faculty.stanford.edu/~uno/programs/frng.f>
##'     \item \code{6}: L'Ecuyer's 1999 algorithm - 64-bits (L'Ecuyer, 1999).
##'     FORTRAN code adapted from <https://jblevins.org/mirror/amiller/lfsr258.f90>
##'   }
##'   For more details on these algorithms see \code{\link[base]{Random}} and references
##'   therein.
##'
##' \item \code{debug} logical, if \code{TRUE} the output from FORTRAN is return (for
##' debuggin purposes).  Default is \code{FALSE} for all models.
##'}
##'
##'
##' @param model character; one of \code{"BARFIMA"}, \code{"GARFIMA"},
##' \code{"KARFIMA"}, \code{"BARC"}.
##'
##' @param complete logical; if \code{FALSE} the function returns only the simulated
##' time series yt, otherwise, additional time series are provided.
##' Default is \code{FALSE} for all models.
##'
##' @param ...  further arguments passed to the functions, according to
##' the model selected in the argument \code{model}. See \sQuote{Common Arguments}
##'
##' @return
##' The function \code{btsr.sim} returns the simulated time series yt  by default.
##' If \code{complete = TRUE}, a list with the following components
##' is returned instead:
##' \itemize{
##' \item \code{model}: character; one of \code{"BARFIMA"}, \code{"GARFIMA"},
##' \code{"KARFIMA"}, \code{"BARC"}. (same as the input argument)
##'
##' \item \code{yt}: the simulated time series
##'
##' \item \code{gyt}: the transformed time series \eqn{g2(y_t)}
##'
##' \item \code{mut}: the conditional mean
##'
##' \item \code{etat}: the linear predictor \eqn{g(\mu_t)}
##'
##' \item \code{error}: the error term \eqn{r_t}
##'
##' \item \code{xreg}: the regressors (if included in the model).
##'
##' \item \code{debug}: the output from FORTRAN (if requested).
##' }
##'
##' @seealso
##' \code{\link{BARFIMA.sim}},  \code{\link{GARFIMA.sim}},
##' \code{\link{KARFIMA.sim}},  \code{\link{BARC.sim}}
##'
##' @references
##' Knuth, D. E. (2002). The Art of Computer Programming. Volume 2,
##' third edition, ninth printing.
##'
##' L'Ecuyer, P. (1999). Good parameters and implementations for combined
##' multiple recursive random number generators. Operations Research, 47,
##' 159-164. <doi:10.1287/opre.47.1.159.>
##'
##' Matsumoto, M. and Nishimura, T. (1998). Mersenne Twister: A 623-dimensionally
##' equidistributed uniform pseudo-random number generator, ACM Transactions on
##' Modeling and Computer Simulation, 8, 3-30.
##'
##' Wichmann, B. A. and Hill, I. D. (1982). Algorithm AS 183: An Efficient
##' and Portable Pseudo-random Number Generator. Applied Statistics, 31, 188-190;
##' Remarks: 34, 198 and 35, 89. <doi:10.2307/2347988.>
##'
##'
##' @examples
##' # Generating a Beta model were mut does not vary with time
##' # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
##'
##' y <- btsr.sim(model= "BARFIMA", linkg = "linear",
##'                n = 1000, seed = 2021,
##'                coefs = list(alpha = 0.2, nu = 20))
##' hist(y)
##'
##' @export
##'
##' @md
btsr.sim <- function(model, complete = FALSE,...){

  temp <- list(...)
  cf = list(d = 0)
  if(!is.null(temp[['coefs']])){
    cf[names(temp$coefs)] <- temp$coefs
		cf$d = 0
	}	
  switch(EXPR = model[1],
         BARFIMA = BARFIMA.sim(complete = complete,...),
         GARFIMA = GARFIMA.sim(complete = complete,...),
         KARFIMA = KARFIMA.sim(complete = complete,...),
         UWARFIMA = UWARFIMA.sim(complete = complete,...),
         BARMA = BARFIMA.sim(coefs = cf,...),
         GARMA = GARFIMA.sim(coefs = cf,...),
         KARMA = KARFIMA.sim(coefs = cf,...),
         UWARMA = UWARFIMA.sim(coefs = cf,...),
         BARC = BARC.sim(complete = complete,...),
         "not available")
}


##' @rdname btsr.functions
##' @order 3
##'
##' @details
##'
##' The function \code{btsr.extract} allows the user to extract the
##' components \eqn{y_t}, \eqn{\mu_t},  \eqn{\eta_t = g(\mu_t)}, \eqn{r_t},
##' the log-likelihood, and the vectors and matrices used to calculate the
##' score vector and the information matrix associated to a given set of parameters.
##'
##' This function can be used by any user to create an objective function
##' that can be passed to optimization functions not available in BTSR Package.
##' At this point, there is no other use for which this function was intended.
##'
##' # Common Arguments
##'
##' ## Extracting Function
##'
##' Common arguments passed through \code{"..."} in \code{btsr.extract} are:
##'
##'\itemize{
##' \item \code{yt} a numeric vector with the observed time series. If missing, an error
##' message is issued.
##'
##' \item \code{xreg} optionally, a vector or matrix with the regressor's values.
##' Default is \code{NULL} for all models.
##'
##' \item \code{nnew} optionally, the number of out-of sample predicted values required.
##' Default is 0 for all models.
##'
##' \item \code{xnew} a vector or matrix, with \code{nnew} observations of the
##' regressors observed/predicted values corresponding to the period of
##' out-of-sample forecast. If \code{xreg = NULL}, \code{xnew} is ignored.
##'
##' \item \code{p} a non-negative integer. The order of AR polynomial.
##' If missing, the value of \code{p} is calculated from length(coefs$phi)
##' and length(fixed.values$phi).
##'
##' \item \code{q,r} a non-negative integer. The order of the MA polynomial and
##' the size of the vector of parameters for the map function (BARC only).
##' If missing, the argument is calcualted based on length(coefs$theta)
##' and length(fixed.values$theta).
##'
##'  \item \code{coefs} a list with the coefficients of the model. Each model has
##'  its default. Passing both, \code{coefs} and \code{fixed.values} empty
##'  will result in an error. The arguments in this list are
##'  \itemize{
##'   \item \code{alpha} a numeric value corresponding to the intercept.
##'   If missing, will be set as zero.
##'
##'   \item \code{beta} a vector of coefficients corresponding to the
##'   regressors in \code{xreg}. If \code{xreg} is provided but \code{beta} is
##'   missing in the \code{coefs} list, an error message is issued.
##'
##'   \item \code{phi} a vector with the non-fixed values in the vector of
##'   AR coefficients.
##'
##'   \item \code{nu} the dispersion parameter. If missing, an error message is issued.
##'
##'   \item \code{theta, d, u0} model specific arguments. See the documentation
##'   corresponding to each model.
##'  }
##'
##' \item \code{lags} optionally, a list with the lags that the values in \code{coefs} correspond to.
##' The names of the entries in this list must match the ones in \code{coefs}.
##' For one dimensional coefficients, the \code{lag} is obviously always 1 and can
##' be suppressed. An empty list indicates that either the argument \code{fixed.lags}
##' is provided or all lags must be used.
##'
##' \item \code{fixed.values} optionally, a list with the values of the coefficients
##' that are fixed. By default, if a given vector (such as the vector of AR coefficients)
##' has fixed values and the corresponding entry in this list is empty, the fixed values
##' are set as zero. The names of the entries in this list must match the ones
##' in \code{coefs}.
##'
##' \item \code{fixed.lags} optionally, a list with the lags that the fixed values
##' in \code{fixed.values} correspond to. The names of the entries in this list must
##' match the ones in \code{fixed.values}. ##' For one dimensional coefficients, the
##' \code{lag} is obviously always 1 and can be suppressed. If an empty list is provided
##' and the model has fixed lags, the argument \code{lags} is used as reference.
##'
##'  \item \code{y.start} optionally, a initial value for yt (to be used
##'  in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{g_2(y_t) = 0}, for \eqn{t < 1}.
##'
##'  \item \code{xreg.start} optionally, a vector of initial value for xreg
##'  (to be used in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{X_t = 0}, for \eqn{t < 1}. If \code{xregar = FALSE} this argument
##'  is ignored.
##'
##'  \item \code{xregar} logical; indicates if xreg is to be included in the
##'  AR part of the model.  See \sQuote{The BTSR structure}. Default is \code{TRUE}.
##'
##'  \item \code{error.scale} the scale for the error term. See also \sQuote{The BTSR structure}.
##'  Each model has its default.
##'
##'  \item \code{inf} the truncation point for infinite sums. Default is 1.
##'  BARC models do not have this argument.
##'
##'  \item \code{m} a non-negative integer indicating the starting time for the sum of the
##'  partial log-likelihoods, that is \eqn{\ell = \sum_{t = m+1}^n \ell_t}. Default is
##'  0.
##'
##'  \item \code{linkg} character or a two character vector indicating which
##'  links must be used in the model.  See \sQuote{The BTSR structure}.
##'  If only one value is provided, the same link is used for \eqn{mu_t} and
##'  for \eqn{y_t} in the AR part of the model. Each model has its default.
##'
##' \item \code{llk} logical, if \code{TRUE} the value of the log-likelihood function
##' is returned. Default is \code{TRUE} for all models.
##'
##' \item \code{sco} logical, if \code{TRUE} the score vector is returned.
##' Default is \code{FALSE} for all models.
##'
##' \item \code{info} logical, if \code{TRUE} the information matrix is returned.
##' Default is \code{FALSE} for all models.
##'
##' \item \code{extra} logical, if \code{TRUE} the matrices and vectors used to
##' calculate the score vector and the information matrix are returned.
##' Default is \code{FALSE} for all models.
##'
##' \item \code{debug} logical, if \code{TRUE} the output from FORTRAN is return (for
##' debuggin purposes).  Default is \code{FALSE} for all models.
##'}
##'
##'
##' @param model character; one of \code{"BARFIMA"}, \code{"GARFIMA"},
##' \code{"KARFIMA"}, \code{"BARC"}.
##' @param ...  further arguments passed to the functions, according to
##' the model selected in the argument \code{model}. See \sQuote{Common Arguments}
##'
##' @return
##' The function \code{btsr.extract} returns a list with the following components.
##' Each particular model can have additional components in this list.
##'
##' \itemize{
##' \item \code{model}: character; one of \code{"BARFIMA"}, \code{"GARFIMA"},
##' \code{"KARFIMA"}, \code{"BARC"}. (same as the input argument)
##'
##' \item \code{coefs}: the coefficients of the model passed through the
##' \code{coefs} argument
##'
##' \item \code{yt}: the observed time series
##'
##' \item \code{gyt}: the transformed time series \eqn{g_2(y_t)}
##'
##' \item \code{mut}: the conditional mean
##'
##' \item \code{etat}: the linear predictor \eqn{g_1(\mu_t)}
##'
##' \item \code{error}: the error term \eqn{r_t}
##'
##' \item \code{xreg}: the regressors (if included in the model).
##'
##' \item \code{forecast}: the out-of-sample forecast (if requested).
##'
##' \item \code{xnew}: the observations of the regressors observed/predicted
##' values corresponding to the period of out-of-sample forecast.
##' Only inlcudes if \code{xreg} is not \code{NULL} and \code{nnew > 0}.
##'
##' \item \code{sll}: the sum of the conditional log-likelihood (if requested)
##'
##' \item \code{sco}: the score vector  (if requested)
##'
##' \item \code{info}: the information matrix  (if requested)
##'
##' \item \code{Drho}, \code{T}, \code{E}, \code{h}: additional matrices and vectors
##' used to calculate the score vector and the information matrix.  (if requested)
##'
##' \item \code{yt.new}: the out-of-sample forecast  (if requested)
##'
##' \item \code{out.Fortran}: FORTRAN output  (if requested)
##'
##' }
##'
##' @seealso
##' \code{\link{BARFIMA.extract}},  \code{\link{GARFIMA.extract}},
##' \code{\link{KARFIMA.extract}},  \code{\link{BARC.extract}}
##'
##' @examples
##'  #------------------------------------------------------------
##'  # Generating a Beta model were mut does not vary with time
##'  # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
##'  #------------------------------------------------------------
##'
##'  m1 <- btsr.sim(model= "BARFIMA", linkg = "linear",
##'                 n = 100, seed = 2021, complete = TRUE,
##'                 coefs = list(alpha = 0.2, nu = 20))
##'
##'  #------------------------------------------------------------
##'  #  Extracting the conditional time series given yt and
##'  #  a set of parameters
##'  #------------------------------------------------------------
##'
##'  # Assuming that all coefficients are non-fixed
##'  e1 = btsr.extract(model = "BARFIMA", yt = m1$yt,
##'                   coefs = list(alpha = 0.2, nu = 20),
##'                   link = "linear", llk = TRUE,
##'                   sco = TRUE, info = TRUE)
##'
##'  # Assuming that all coefficients are fixed
##'  e2 = btsr.extract(model = "BARFIMA", yt = m1$yt,
##'                   fixed.values = list(alpha = 0.2, nu = 20),
##'                   link = "linear", llk = TRUE,
##'                   sco = TRUE, info = TRUE)
##'
##'  # Assuming at least one fixed coefficient and one non-fixed
##'  e3 = btsr.extract(model = "BARFIMA", yt = m1$yt,
##'                   fixed.values = list(alpha = 0.2, nu = 20),
##'                   link = "linear", llk = TRUE,
##'                   sco = TRUE, info = TRUE)
##'  e4 = btsr.extract(model = "BARFIMA", yt = m1$yt,
##'                   fixed.values = list(alpha = 0.2, nu = 20),
##'                   link = "linear", llk = TRUE,
##'                   sco = TRUE, info = TRUE)
##'
##'  #----------------------------------------------------
##'  # comparing the simulated and the extracted values
##'  #----------------------------------------------------
##'  cbind(head(m1$mut), head(e1$mut), head(e2$mut), head(e3$mut), head(e4$mut))
##'
##'  #----------------------------------------------------
##'  # comparing the log-likelihood values obtained (must be the all equal)
##'  #----------------------------------------------------
##'  c(e1$sll, e2$sll, e3$sll, e4$sll)
##'
##'  #----------------------------------------------------
##'  # comparing the score vectors:
##'  #----------------------------------------------------
##'  # - e1 must have 2 values: dl/dmu and dl/dnu
##'  # - e2 must be empty
##'  # - e3 and e4 must have one value corresponding
##'  #    to the non-fixed coefficient
##'  #----------------------------------------------------
##'   e1$score
##'   e2$score
##'   e3$score
##'   e4$score
##'
##'  #----------------------------------------------------
##'  # comparing the information matrices.
##'  #----------------------------------------------------
##'  # - e1 must be a 2x2 matrix
##'  # - e2 must be empty
##'  # - e3 and e4 must have one value corresponding
##'  #    to the non-fixed coefficient
##'  #----------------------------------------------------
##'   e1$info.Matrix
##'   e2$info.Matrix
##'   e3$info.Matrix
##'   e4$info.Matrix
##'
##' @export
##'
##' @md
btsr.extract <- function(model,...){

  temp <- list(...)
  fv = list(d = 0)
  if( !is.null(temp[['fixed.values']]))
    fv[names(temp$fixed.values)] <- temp$fixed.values

  switch(EXPR = model[1],
         BARFIMA = BARFIMA.extract(...),
         GARFIMA = GARFIMA.extract(...),
         KARFIMA = KARFIMA.extract(...),
         UWARFIMA = UWARFIMA.extract(...),
         BARMA = BARFIMA.extract(fixed.values = fv,...),
         GARMA = GARFIMA.extract(fixed.values = fv,...),
         KARMA = KARFIMA.extract(fixed.values = fv,...),
         UWARMA = UWARFIMA.extract(fixed.values = fv,...),
         BARC = BARC.extract(...),
         "not available")
}


##' @rdname btsr.functions
##' @order 4
##'
##' @details
##' The function \code{btsr.fit} fits a BTSR model to a given univariate time
##' series. For now, available optimization algorithms are \code{"L-BFGS-B"} and
##' \code{"Nelder-Mead"}. Both methods accept bounds for the parameters. For
##' \code{"Nelder-Mead"}, bounds are set via parameter transformation.
##'
##' # Common Arguments
##'
##' ## Fitting Function
##'
##' Common arguments passed through \code{"..."} in \code{btsr.fit} are the same as
##' in \code{\link{btsr.extract}} plus the following:
##'
##'\itemize{
##'
##'  \item \code{d} logical, if \code{TRUE}, the parameter \code{d} is included
##'  in the model either as fixed or non-fixed. If \code{d = FALSE} the value is
##'  fixed as 0. The default is \code{TRUE} for all models, except BARC that does
##'  not have this parameter.
##'
##'  \item \code{start} a list with the starting values for the non-fixed coefficients
##'  of the model. If an empty list is provided, the function \code{\link{coefs.start}}
##'  is used to obtain starting values for the parameters.
##'
##' \item \code{ignore.start} logical,  if starting values are not provided, the
##' function uses the default values and \code{ignore.start} is ignored.
##' In case starting values are provided and \code{ignore.start = TRUE}, those
##' starting values are ignored and recalculated. The default is \code{FALSE}.
##'
##' \item \code{lower, upper} optionally, list with the lower and upper bounds for the
##' parameters. The names of the entries in these lists must match the ones
##' in \code{start}. The default is to assume that the parameters are unbounded.
##' Only the bounds for bounded parameters need to be specified.
##'
##' \item \code{control} a list with configurations to be passed to the
##' optimization subroutines. Missing arguments will receive default values. See
##' \cite{\link{fit.control}}.
##'
##' \item \code{report} logical, if \code{TRUE} the summary from model estimation is
##' printed and \code{info} is automatically set to \code{TRUE}. Default is \code{TRUE}.
##'}
##'
##'
##' @param model character; one of \code{"BARFIMA"}, \code{"GARFIMA"},
##' \code{"KARFIMA"}, \code{"BARC"}.
##' @param ...  further arguments passed to the functions, according to
##' the model selected in the argument \code{model}. See \sQuote{Common Arguments}
##'
##' @return
##' The function \code{btsr.fit} returns a list with the following components.
##' Each particular model can have additional components in this list.
##'
##' \itemize{
##' \item \code{model}: character; one of \code{"BARFIMA"}, \code{"GARFIMA"},
##' \code{"KARFIMA"}, \code{"BARC"}. (same as the input argument)
##'
##' \item \code{convergence}: An integer code. 0 indicates successful completion.
##'  The error codes depend on the algorithm used.
##'
##' \item \code{message}: A character string giving any additional information
##' returned by the optimizer, or NULL.
##'
##' \item \code{counts}: an integer giving the number of function evaluations.
##'
##' \item \code{control}: a list of control parameters.
##'
##' \item \code{start}: the starting values used by the algorithm.
##'
##' \item \code{coefficients}: 	The best set of parameters found.
##'
##' \item \code{n}: the sample size used for estimation.
##'
##' \item \code{series}: the observed time series
##'
##' \item \code{gyt}: the transformed time series \eqn{g_2(y_t)}
##'
##' \item \code{fitted.values}: the conditional mean, which corresponds to
##' the in-sample forecast, also denoted fitted values
##'
##' \item \code{etat}: the linear predictor \eqn{g_1(\mu_t)}
##'
##' \item \code{error.scale}: the scale for the error term.
##'
##' \item \code{error}: the error term \eqn{r_t}
##'
##' \item \code{residuals}: the observed minus the fitted values. The same as
##' the \code{error} term if \code{error.scale = 0}.
##'
##' \item \code{sll}: the sum of the conditional log-likelihood (if requested)
##'
##' \item \code{info.Matrix}: the information matrix  (if requested)
##'
##' \item \code{configs}: a list with the configurations adopted to fit the model.
##' This information is used by the prediction function.
##'
##' \item \code{out.Fortran}: FORTRAN output  (if requested)
##'
##' \item \code{call}: a string with the description of the fitted model.
##'
##' }
##'
##' @seealso
##' \code{\link{BARFIMA.fit}},  \code{\link{GARFIMA.fit}},
##' \code{\link{KARFIMA.fit}},  \code{\link{BARC.fit}}
##'
##' @examples
##'
##' # Generating a Beta model were mut does not vary with time
##' # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
##'
##' y <- btsr.sim(model= "BARFIMA", linkg = "linear",
##'                n = 100, seed = 2021,
##'                coefs = list(alpha = 0.2, nu = 20))
##'
##' # fitting the model
##' f <- btsr.fit(model = "BARFIMA", yt = y, report = TRUE,
##'              start = list(alpha = 0.5, nu = 10),
##'              linkg = "linear", d = FALSE)
##'
##' @export
##'
##' @md

btsr.fit <- function(model,...){

  switch(EXPR = model[1],
         BARFIMA = BARFIMA.fit(...),
         GARFIMA = GARFIMA.fit(...),
         KARFIMA = KARFIMA.fit(...),
         UWARFIMA = UWARFIMA.fit(...),
         BARMA = BARFIMA.fit(d = FALSE,...),
         GARMA = GARFIMA.fit(d = FALSE,...),
         KARMA = KARFIMA.fit(d = FALSE,...),
         UWARMA = UWARFIMA.fit(d = FALSE,...),
         BARC = BARC.fit(...),
         "not available")
}
