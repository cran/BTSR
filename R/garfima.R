##----------------------------------------------------------
##   GARFIMA MODELS
##----------------------------------------------------------
##' @title
##' Functions to simulate, extract components and fit GARFIMA models
##'
##' @name GARFIMA.functions
##' @order 1
##'
##' @description
##' These functions can be used to simulate, extract components
##' and fit any model of the class \code{garfima}. A model with
##' class \code{garfima} is a special case of a model with class \code{btsr} .
##' See \sQuote{The BTSR structure} in \code{\link{btsr.functions}} for
##' more details on the general structure.
##'
##' The \eqn{\gamma}ARMA model, the gamma regression and a i.i.d. sample
##' from a gamma distribution can be obtained as special cases.
##' See \sQuote{Details}.
##'
##' @details
##' The \eqn{\gamma}ARMA model  and the gamma regression can be
##' obtained as special cases of the \eqn{\gamma}ARFIMA model.
##'
##' \itemize{
##'   \item \eqn{\gamma}ARFIMA: is obtained by default.
##'
##'   \item \eqn{\gamma}ARMA: is obtained by setting \code{d = 0}.
##'
##'   \item gamma regression: is obtained by setting \code{p = 0}, \code{q = 0}
##'   and \code{d = FALSE}. The \code{error.scale} is irrelevant.
##'   The second argument in \code{linkg} is irrelevant.
##'
##'   \item an i.i.d. sample from a Gamma distribution with parameters
##'   \code{shape} and \code{scale} (compatible with the one from \code{\link{rgamma}})
##'   is obtained by  setting \code{linkg = "linear"}, \code{p = 0}, \code{q = 0},
##'   \code{coefs$d = 0}, \code{d = FALSE} and, in the coefficient list,
##'   \code{alpha = shape*scale} and \code{nu = shape}. (\code{error.scale} and
##'   \code{xregar} are irrelevant)
##'}
##'
##' @md
NULL
#> NULL


##' @rdname GARFIMA.functions
##' @order 2
##'
##' @details
##' The function \code{GARFIMA.sim} generates a random sample from a \eqn{\gamma}ARFIMA(p,d,q)
##' model.
##'
##' @param n a strictly positive integer. The sample size of yt (after burn-in).
##'  Default is 1.
##'
##' @param burn a non-negative integer. The length of the "burn-in" period. Default is 0.
##'
##' @param xreg optionally, a vector or matrix of external regressors.
##'  For simulation purposes, the length of xreg must be \code{n+burn}.
##'  Default is \code{NULL}. For extraction or fitting purposes, the length
##'  of \code{xreg} must be the same as the length of the observed time series
##'  \eqn{y_t}.
##'
##' @param coefs a list with the coefficients of the model. An empty list will result
##' in an error. The arguments that can be passed through this list are:
##' \itemize{
##'   \item \code{alpha} optionally, a numeric value corresponding to the intercept.
##'    If the argument is missing, it will be treated as zero. See
##'     \sQuote{The BTSR structure} in \code{\link{btsr.functions}}.
##'
##'   \item \code{beta} optionally, a vector of coefficients corresponding to the
##'   regressors in \code{xreg}. If \code{xreg} is provided but \code{beta} is
##'   missing in the \code{coefs} list, an error message is issued.
##'
##'   \item \code{phi} optionally, for the simulation function this must be a vector
##'   of size \eqn{p}, corresponding to the autoregressive coefficients
##'   (including the ones that are zero), where \eqn{p} is the AR order. For
##'   the extraction and fitting functions, this is a vector with the non-fixed
##'   values in the vector of autoregressive coefficients.
##'
##'   \item \code{theta} optionally, for the simulation function this must be a vector
##'   of size \eqn{q}, corresponding to the moving average coefficients
##'   (including the ones that are zero), where \eqn{q} is the MA order. For
##'   the extraction and fitting functions, this is a vector with the non-fixed
##'   values in the vector of moving average coefficients.
##'
##'   \item \code{d} optionally, a numeric value corresponding to the long memory
##'   parameter. If the argument is missing, it will be treated as zero.
##'
##'   \item \code{nu} the dispersion parameter. If missing, an error message is issued.
##'
##' }
##'
##' @param y.start optionally, an initial value for yt (to be used
##'  in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{g_2(y_t) = 0}, for \eqn{t < 1}.
##'
##' @param xreg.start optionally, a vector of initial value for xreg
##'  (to be used in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{X_t = 0}, for \eqn{t < 1}. If \code{xregar = FALSE} this argument
##'  is ignored.
##'
##' @param xregar logical; indicates if xreg is to be included in the
##'  AR part of the model.  See \sQuote{The BTSR structure}. Default is \code{TRUE}.
##'
##' @param error.scale  the scale for the error term. See \sQuote{The BTSR structure}
##' in \code{\link{btsr.functions}}. Default is 0.
##'
##' @param complete logical; if \code{FALSE} the function returns only the simulated
##' time series yt, otherwise, additional time series are provided.
##' Default is \code{FALSE}
##'
##' @param inf the truncation point for infinite sums. Default is 1,000.
##'  In practice, the Fortran subroutine uses \eqn{inf = q}, if \eqn{d = 0}.
##'
##' @param linkg character or a two character vector indicating which
##'  links must be used in the model.  See \sQuote{The BTSR structure}
##'  in \code{\link{btsr.functions}} for details and \code{\link{link.btsr}}
##'  for valid links. If only one value is provided, the same link is used
##'  for \eqn{mu_t} and for \eqn{y_t} in the AR part of the model.
##'  Default is \code{c("log", "log")}. For the linear link, the constant
##'  will be always 1.
##'
##' @param seed optionally, an integer which gives the value of the fixed
##'  seed to be used by the random number generator. If missing, a random integer
##'  is chosen uniformly from 1,000 to 10,000.
##'
##' @param rngtype optionally, an integer indicating which random number generator
##'  is to be used. Default is 2: the Mersenne Twister algorithm. See \sQuote{Common Arguments}
##'  in \code{\link{btsr.functions}}.
##'
##' @param debug logical, if \code{TRUE} the output from FORTRAN is return (for
##' debugging purposes).  Default is \code{FALSE} for all models.
##'
##' @return
##' The function \code{GARFIMA.sim} returns the simulated time series yt  by default.
##' If \code{complete = TRUE}, a list with the following components
##' is returned instead:
##' \itemize{
##' \item \code{model}: string with the text \code{"GARFIMA"}
##'
##' \item \code{yt}: the simulated time series
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
##'
##' }
##'
##' @seealso
##' \code{\link{btsr.sim}}
##'
##' @examples
##' # Generating a Gamma model were mut does not vary with time
##' # yt ~ Gamma(a,b), a = nu (shape), b = mu/nu (scale)
##'
##' y <- GARFIMA.sim(linkg = "linear", n = 1000, seed = 2021,
##'                  coefs = list(alpha = 0.2, nu = 20))
##' hist(y)
##'
##' @export
##'
##' @md
GARFIMA.sim <- function(n = 1, burn = 0, xreg = NULL,
                        coefs = list(alpha = 0, beta = NULL, phi = NULL,
                                     theta = NULL, d = 0, nu = 20),
                        y.start = NULL, xreg.start = NULL,
                        xregar = TRUE, error.scale = 0, complete = FALSE,
                        inf = 1000, linkg = c("log", "log"), seed = NULL,
                        rngtype = 2, debug = FALSE){

  ####----------------------------------
  #### checking required parameters:
  ####----------------------------------
  if(is.null(coefs)) stop("coefs missing with no default")
  if(!"list" %in% class(coefs)) stop("coefs must be a list")

  ####----------------------------------
  #### checking configurations:
  ####----------------------------------
  cf <- .sim.configs(model = "GARFIMA",  xreg = xreg,
                     y.start = y.start, xreg.start = xreg.start,
                     linkg = linkg, n = n, burn = burn,
                     coefs = coefs, xregar = xregar,
                     error.scale = error.scale, seed = seed,
                     rngtype = rngtype, y.default = 0)

  out <- .btsr.sim(model = "GARFIMA", inf = inf, configs = cf,
                   complete = complete, debug = debug)
  class(out) <- c(class(out), "garfima")
  invisible(out)
}


##' @rdname GARFIMA.functions
##' @order 3
##'
##' @details
##'
##' The function \code{GARFIMA.extract} allows the user to extract the
##' components \eqn{y_t}, \eqn{\mu_t},  \eqn{\eta_t = g(\mu_t)}, \eqn{r_t},
##' the log-likelihood, and the vectors and matrices used to calculate the
##' score vector and the information matrix associated to a given set of parameters.
##'
##' This function can be used by any user to create an objective function
##' that can be passed to optimization algorithms not available in the BTSR Package.
##'
##' @param yt a numeric vector with the observed time series. If missing, an error
##' message is issued.
##'
##' @param nnew optionally, the number of out-of sample predicted values required.
##' Default is 0.
##'
##' @param xnew  a vector or matrix, with \code{nnew} observations of the
##' regressors observed/predicted values corresponding to the period of
##' out-of-sample forecast. If \code{xreg = NULL}, \code{xnew} is ignored.
##'
##' @param p a non-negative integer. The order of AR polynomial.
##' If missing, the value of \code{p} is calculated from length(coefs$phi)
##' and length(fixed.values$phi). For fitting, the default is 0.
##'
##' @param q a non-negative integer. The order of the MA polynomial.
##' If missing, the value of \code{q} is calculated from length(coefs$theta)
##' and length(fixed.values$theta). For fitting, the default is 0.
##'
##' @param lags optionally, a list with the lags that the values in \code{coefs} correspond to.
##' The names of the entries in this list must match the ones in \code{coefs}.
##' For one dimensional coefficients, the \code{lag} is obviously always 1 and can
##' be suppressed. An empty list indicates that either the argument \code{fixed.lags}
##' is provided or all lags must be used.
##'
##' @param fixed.values optionally, a list with the values of the coefficients
##' that are fixed. By default, if a given vector (such as the vector of AR coefficients)
##' has fixed values and the corresponding entry in this list is empty, the fixed values
##' are set as zero. The names of the entries in this list must match the ones
##' in \code{coefs}.
##'
##' @param fixed.lags optionally, a list with the lags that the fixed values
##' in \code{fixed.values} correspond to. The names of the entries in this list must
##' match the ones in \code{fixed.values}. ##' For one dimensional coefficients, the
##' \code{lag} is obviously always 1 and can be suppressed. If an empty list is provided
##' and the model has fixed lags, the argument \code{lags} is used as reference.
##'
##' @param m a non-negative integer indicating the starting time for the sum of the
##' partial log-likelihoods, that is \eqn{\ell = \sum_{t = m+1}^n \ell_t}. Default is
##' 0.
##'
##' @param llk logical, if \code{TRUE} the value of the log-likelihood function
##' is returned. Default is \code{TRUE}.
##'
##' @param sco logical, if \code{TRUE} the score vector is returned.
##' Default is \code{FALSE}.
##'
##' @param info logical, if \code{TRUE} the information matrix is returned.
##' Default is \code{FALSE}. For the fitting function, \code{info} is automatically
##' set to \code{TRUE} when \code{report = TRUE}.
##'
##' @param extra logical, if \code{TRUE} the matrices and vectors used to
##' calculate the score vector and the information matrix are returned.
##' Default is \code{FALSE}.
##'
##' @return
##' The function \code{GARFIMA.extract} returns a list with the following components.
##'
##' \itemize{
##'  \item \code{model}: string with the text \code{"GARFIMA"}
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
##' }
##'
##' @seealso
##' \code{\link{btsr.extract}}
##'
##' @examples
##'  #------------------------------------------------------------
##'  # Generating a Gamma model were mut does not vary with time
##'  # yt ~ Gamma(a,b), a = nu (shape), b = mu/nu (scale)
##'  #------------------------------------------------------------
##'
##'  m1 <- GARFIMA.sim(linkg = "linear",n = 100,
##'                    complete = TRUE, seed = 2021,
##'                    coefs = list(alpha = 0.2, nu = 20))
##'
##'  #------------------------------------------------------------
##'  #  Extracting the conditional time series given yt and
##'  #  a set of parameters
##'  #------------------------------------------------------------
##'
##'  # Assuming that all coefficients are non-fixed
##'  e1 = GARFIMA.extract(yt = m1$yt, coefs = list(alpha = 0.2, nu = 20),
##'                       link = "linear", llk = TRUE,
##'                       sco = TRUE, info = TRUE)
##'
##'  #----------------------------------------------------
##'  # comparing the simulated and the extracted values
##'  #----------------------------------------------------
##'  cbind(head(m1$mut), head(e1$mut))
##'
##'  #---------------------------------------------------------
##'  # the log-likelihood, score vector and information matrix
##'  #---------------------------------------------------------
##'  e1$sll
##'  e1$score
##'  e1$info.Matrix
##'
##' @export
##' @md
GARFIMA.extract <- function(yt, xreg = NULL, nnew = 0, xnew  = NULL,
                            p, q, coefs = list(),lags = list(),
                            fixed.values = list(), fixed.lags = list(),
                            y.start = NULL, xreg.start = NULL,
                            xregar = TRUE,  error.scale = 0, inf = 1000, m = 0,
                            linkg = c("log","log"), llk = TRUE, sco = FALSE,
                            info = FALSE, extra = FALSE,  debug = FALSE){

  if(is.null(coefs) & is.null(fixed.values))
    stop("Please, provide a list of coefficients")
  if(!is.null(coefs)){
    if(! "list" %in% class(coefs)) stop("coefs must be a list")}
  if(!is.null(fixed.values)){
    if(! "list" %in% class(fixed.values)) stop("fixed.values must be a list")}
  else{ fixed.values <- list()}

  if(missing(p)) p = length(coefs$phi) + length(fixed.values$phi)
  if(missing(q)) q = length(coefs$theta) + length(fixed.values$theta)

  cf <-  .extract.configs(model = "GARFIMA", yt = yt, y.start = y.start,
                          y.lower = 0, y.upper = Inf, openIC = c(TRUE, TRUE),
                          xreg = xreg, xnew = xnew, nnew = nnew,
                          xreg.start = xreg.start, linkg = linkg,
                          p = p, q = q, inf = inf, m = m, xregar = xregar,
                          error.scale = error.scale, coefs = coefs,
                          lags = lags, fixed.values = fixed.values,
                          fixed.lags = fixed.lags, llk = llk, sco = sco,
                          info = info, extra = extra)

  out <- .btsr.extract(model = "GARFIMA", yt = yt, configs = cf, debug = debug)
  class(out) <- c(class(out), "garfima")

  invisible(out)
}



##' @rdname GARFIMA.functions
##' @order 4
##'
##' @details
##' The function \code{GARFIMA.fit} fits a GARFIMA model to a given univariate time
##' series. For now, available optimization algorithms are \code{"L-BFGS-B"} and
##' \code{"Nelder-Mead"}. Both methods accept bounds for the parameters. For
##' \code{"Nelder-Mead"}, bounds are set via parameter transformation.
##'
##'
##' @param d logical, if \code{TRUE}, the parameter \code{d} is included
##'  in the model either as fixed or non-fixed. If \code{d = FALSE} the value is
##'  fixed as 0. The default is \code{TRUE}.
##'
##' @param start a list with the starting values for the non-fixed coefficients
##'  of the model. If an empty list is provided, the function \code{\link{coefs.start}}
##'  is used to obtain starting values for the parameters.
##'
##' @param ignore.start logical,  if starting values are not provided, the
##' function uses the default values and \code{ignore.start} is ignored.
##' In case starting values are provided and \code{ignore.start = TRUE}, those
##' starting values are ignored and recalculated. The default is \code{FALSE}.
##'
##' @param lower  optionally, list with the lower bounds for the
##' parameters. The names of the entries in these lists must match the ones
##' in \code{start}. The default is to assume that the parameters have no lower
##' bound except for \code{nu}, for which de default is 0. Only the bounds for
##' bounded parameters need to be specified.
##'
##' @param upper optionally, list with the upper bounds for the
##' parameters. The names of the entries in these lists must match the ones
##' in \code{start}. The default is to assume that the parameters have no upper
##' bound. Only the bounds for bounded parameters need to be specified.
##'
##' @param control a list with configurations to be passed to the
##' optimization subroutines. Missing arguments will receive default values. See
##' \cite{\link{fit.control}}.
##'
##' @param report logical, if \code{TRUE} the summary from model estimation is
##' printed and \code{info} is automatically set to \code{TRUE}. Default is \code{TRUE}.
##'
##' @param ... further arguments passed to the internal functions.
##'
##' @return
##' The function \code{btsr.fit} returns a list with the following components.
##' Each particular model can have additional components in this list.
##'
##' \itemize{
##'  \item \code{model}: string with the text \code{"GARFIMA"}
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
##' \item \code{residual}: the observed minus the fitted values. The same as
##' the \code{error} term if \code{error.scale = 0}.
##'
##' \item \code{forecast}: the out-of-sample forecast (if requested).
##'
##' \item \code{xnew}: the observations of the regressors observed/predicted
##' values corresponding to the period of out-of-sample forecast.
##' Only inlcudes if \code{xreg} is not \code{NULL} and \code{nnew > 0}.
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
##' \code{\link{btsr.fit}}
##'
##' @examples
##'
##' # Generating a Beta model were mut does not vary with time
##' # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
##'
##' y <- GARFIMA.sim(linkg = "linear", n = 100, seed = 2021,
##'                coefs = list(alpha = 0.2, nu = 20))
##'
##' # fitting the model
##' f <- GARFIMA.fit(yt = y, report = TRUE,
##'                  start = list(alpha = 0.5, nu = 10),
##'                  linkg = "linear", d = FALSE)
##'
##' @export
##'
##' @md
GARFIMA.fit <- function(yt, xreg = NULL, nnew = 0, xnew = NULL,
                        p = 0, d = TRUE,  q = 0, m = 0, inf = 1000,
                        start = list(), ignore.start = FALSE,
                        lags = list(), fixed.values = list(),
                        fixed.lags = list(), lower = list(nu = 0),
                        upper = list(nu = Inf),  linkg = c("log","log"),
                        sco = TRUE, info = FALSE, extra = FALSE, xregar = TRUE,
                        y.start = NULL, xreg.start = NULL,
                        error.scale = 0, control = list(), report = TRUE,
                        debug = FALSE,...){

  # default values for nu (merge with user provided values)
  lw <- list(nu = 0); up <- list(nu = Inf)
  lw[names(lower)] <- lower; up[names(upper)] <- upper
  lower <- lw; upper <- up

  if(report) info = TRUE
  cf <- .fit.configs(model = "GARFIMA", yt = yt, y.start = y.start,
                     y.lower = 0, y.upper = Inf, openIC = c(TRUE, TRUE),
                     xreg = xreg, xnew = xnew, nnew = nnew,
                     xreg.start = xreg.start, linkg = linkg,
                     p = p, d = d, q = q, inf = inf, m = m,
                     xregar = xregar, error.scale = error.scale,
                     start = start, ignore.start = ignore.start,
                     lags = lags, fixed.values = fixed.values,
                     fixed.lags = fixed.lags, lower = lower,
                     upper = upper, control = control,
                     sco = sco, info = info, extra = extra)
  if(!is.null(cf$conv)) return(invisible(out))

  out <- .btsr.fit(model = "GARFIMA", yt = yt, configs = cf, debug = debug)
  out$call <- .fit.print(model = "GARFIMA", p = cf$p, q = cf$q,
                         d = !(cf$d$nfix == 1 & cf$d$fvalues == 0),
                         nreg = cf$nreg)
  class(out) <- c(class(out), "garfima")
  if(report) print(summary(out))
  invisible(out)
}
