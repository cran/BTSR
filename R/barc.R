##---------------------------------------------------------------------
##
##  Simulating a BARC model
##
##---------------------------------------------------------------------
##' @title
##' Functions to simulate, extract components and fit BARC models
##'
##' @name BARC.functions
##' @order 1
##'
##' @description
##' These functions can be used to simulate, extract components
##' and fit any model of the class \code{barc}. A model with
##' class \code{barc} is a special case of a model with class \code{btsr} .
##' See \sQuote{The BTSR structure} in \code{\link{BARC.functions}} for
##' more details on the general structure. See \sQuote{Details}.
##'
##' @details
##'
##' Neither the beta regression or an i.i.d. sample
##' from a beta distribution can be obtained as special cases of the
##' \eqn{\beta}ARC model since the term \eqn{h(T(U_0))} is always present
##'
##' The model from Pumi et al. (2021) is obtained by setting
##' \code{xregar = TRUE} (so that the regressors are included in the AR
##' part of the model) and using the same link for \eqn{y_t} and \eqn{\mu_t}.
##'
##' # The map function
##'
##' The map function \eqn{T:[0,1] \to [0,1]} is a dynamical system, i.e.,
##' a function, potentially depending on a \eqn{r}-dimensional vector of
##' parameters \eqn{\theta}. Available choices are
##' \itemize{
##'   \item \code{map = 1}, \eqn{\theta = k}, for \eqn{k} integer greater
##'   or equal to 2.
##'
##'   \deqn{T(u) = (ku)(mod 1)}
##'
##'   \item \code{map = 2}, \eqn{0 \le \theta \le 1}
##'
##'   \deqn{T(u) = \frac{u}{\theta}I_( u < \theta) +
##'           \theta\frac{(u - \theta)}{(1 - \theta)}I(u \ge \theta)}
##'
##'   \item \code{map = 3} (logistic map), \eqn{ 0 \le \theta \le 4},
##'
##'   \deqn{T(u) = \theta(1-\theta)}
##'
##'   \item \code{map = 4} (Manneville-Pomeau map), \eqn{0 < \theta < 1}
##'
##'   \deqn{T(u) = (u + u^{1+\theta})(mod 1)}
##'
##'   \item \code{map = 5} (Lasota-Mackey's map),
##'
##'   \deqn{T(u) = \frac{u}{(1 - u)}I(u \le 0.5) + (2u - 1)I(u > 0.5)}
##' }
##' @references
##'
##' Pumi, G.; Prass, T.S. and Souza, R.R. (2021). A dynamic model for
##' double bounded time series with chaotic driven conditional averages.
##' Scandinavian Journal of Statistics. Vol 48 (1), 68-86.
##'
##' @seealso
##' \code{\link{btsr.sim}}, \code{\link{btsr.extract}}, \code{\link{btsr.fit}}
##'
##' @md
NULL
#> NULL


##' @rdname BARC.functions
##' @order 2
##'
##' @details
##' The function \code{BARC.sim} generates a random sample from a
##' \eqn{\beta}ARC(p) model.
##'
##' @param n a strictly positive integer. The sample size of yt (after burn-in).
##'  Default is 1.
##'
##' @param burn a non-negative integer. length of "burn-in" period. Default is 0.
##'
##' @param xreg optionally, a vector or matrix of external regressors.
##'  For simulation purposes, the length of xreg must be \code{n+burn}.
##'  Default is \code{NULL}. For extraction or fitting purposes, the length
##'  of \code{xreg} must be the same as the length of the observed time series
##'  \eqn{y_t}.
##'
##' @param map a non-negative integer from 1 to 5 corresponding to the map function.
##' Default is 4. See \sQuote{The map function}.
##'
##' @param coefs a list with the coefficients of the model. An empty list will result
##' in an error. The arguments that can be passed through this list are:
##' \itemize{
##'   \item \code{alpha} optionally, a numeric value corresponding to the intercept.
##'    If the argument is missing, it will be treated as zero. See
##'    \sQuote{The BTSR structure} in \code{\link{btsr.functions}}.
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
##'   \item \code{theta} the parameter (or vector of parameters) corresponding
##'   to the map function. If \code{map = 5} this value is ignored. For simulation,
##'   purposes, the default is \code{map = 4} and \code{theta = 0.5}.
##'
##'   \item \code{nu} the dispersion parameter. If missing, an error message is issued.
##'
##'   \item \code{u0} a numeric value in the interval \eqn{(0,1)}, corresponding
##'   to the value of the random variable \eqn{U_0}. For simulation purposes, the
##'   default is \code{u0 = pi/4}.
##'
##' }
##'
##' @param y.start optionally, a initial value for yt (to be used
##'  in the recursions). Default is \code{NULL}, in which case, the recursion assumes
##'  that \eqn{g_2(y_t) = 0}, for \eqn{t < 1}.
##'
##' @param xreg.start optionally, a vector of initial value for xreg
##'  (to be used in the recursions). Default is \code{NULL}, in which case, the
##'  recursion assumes that \eqn{X_t = 0}, for \eqn{t < 1}. If \code{xregar = FALSE}
##'  this argument is ignored.
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
##' @param linkg character or a two character vector indicating which
##'  links must be used in the model.  See \sQuote{The BTSR structure}
##'  in \code{\link{btsr.functions}} for details and \code{\link{link.btsr}}
##'  for valid links. If only one value is provided, the same link is used
##'  for \eqn{mu_t} and for \eqn{y_t} in the AR part of the model.
##'  Default is \code{c("linear", "linear")}
##'
##' @param linkh a character indicating which link must be associated to the
##'  the chaotic process.  See \sQuote{The BTSR structure}
##'  in \code{\link{btsr.functions}} for details and \code{\link{link.btsr}}
##'  for valid links. Default is \code{"linear"}.
##'
##' @param ctt.h numeric; the constant to be associated to the link \eqn{h},
##' when \code{linkh = "linear"}. Default is 1.
##'
##' @param seed optionally, an integer which gives the value of the fixed
##'  seed to be used by the random number generator. If missing, a random integer
##'  is chosen uniformly from 1,000 to 10,000.
##'
##' @param rngtype optionally, an integer indicating which random number generator
##'  is to be used. Default is 2. See \sQuote{Common Arguments}
##'  in \code{\link{btsr.functions}}.
##'
##' @param debug logical, if \code{TRUE} the output from FORTRAN is return (for
##' debuggin purposes).  Default is \code{FALSE} for all models.
##'
##' @return
##' The function \code{BARC.sim} returns the simulated time series yt  by default.
##' If \code{complete = TRUE}, a list with the following components
##' is returned instead:
##' \itemize{
##'   \item \code{model}: string with the text \code{"BARC"}
##'
##'   \item \code{yt}: the simulated time series
##'
##'   \item \code{mut}: the conditional mean
##'
##'   \item \code{etat}: the linear predictor \eqn{g(\mu_t)}
##'
##'   \item \code{error}: the error term \eqn{r_t}
##'
##'   \item \code{xreg}: the regressors (if included in the model).
##'
##'   \item \code{debug}: the output from FORTRAN (if requested).
##'
##' }
##'
##' @examples
##' m1 <- BARC.sim(linkg = "linear", linkh = "linear",
##'               n = 100, seed = 2021, complete = TRUE, ctt.h = 0.6,
##'              coefs = list(nu = 15, theta = 0.85, u0 = pi/4))
##'
##' plot.ts(m1$yt)
##' lines(m1$mut, col = "red")
##'
##' @export
##'
##' @md
BARC.sim <- function(n = 1, burn = 0, xreg = NULL, map = 4,
                     coefs = list(alpha = 0, beta = NULL, phi = NULL,
                                  theta = 0.5, nu = 20, u0 = pi/4),
                     y.start = NULL, xreg.start = NULL,
                     xregar = TRUE, error.scale = 0, complete = FALSE,
                     linkg = c("linear","linear"), linkh = "linear",
                     ctt.h = 1, seed = NULL, rngtype = 2, debug = FALSE){

  ##----------------------------------
  ## checking required parameters:
  ##----------------------------------
  if(is.null(coefs)) stop("coefs missing with no default")
  if(!"list" %in% class(coefs)) stop("coefs must be a list")
  if(is.null(coefs$u0)) stop("u0 is missing")

  ##--------------------------------------------
  ## checking the map, theta and the link "h"
  ##--------------------------------------------
  cb <- .barc.configs(map = map, theta = coefs$theta, linkh = linkh)

  ##--------------------------------------------
  ## checking remaining configurations:
  ##--------------------------------------------
  cf <- .sim.configs(model = "BARC",  xreg = xreg,
                     y.start = y.start, xreg.start = xreg.start,
                     linkg = linkg, n = n, burn = burn,
                     coefs = coefs, xregar = xregar,
                     error.scale = error.scale, seed = seed,
                     rngtype = rngtype, y.default = 0)

  out <- .barc.sim(u0 = coefs$u0, map = map, ctt.h = ctt.h,
                   configs = cf, bconfigs = cb, complete = complete,
                   debug = debug)
  class(out) <- c(class(out), "barc")
  invisible(out)
}

##------------------------------------------------------------------------------------
## internal function:  makes the calculations and reports only the relevant variables
##------------------------------------------------------------------------------------
.barc.sim <- function(u0, map, ctt.h, configs, bconfigs, complete, debug){

  out <- .Fortran("simbarcR",
                  n = configs$n,
                  burn = configs$burn,
                  nu = configs$nu,
                  alpha = configs$alpha,
                  nreg = configs$nreg,
                  beta = configs$beta,
                  p = configs$p,
                  phi = configs$phi,
                  r = bconfigs$r,
                  theta = bconfigs$theta,
                  u0 = u0,
                  map = as.integer(map),
                  link = c(configs$linkg, bconfigs$linkh),
                  ctt.h = ctt.h,
                  xreg  = configs$xreg,
                  xregar = configs$xregar,
                  yt = numeric(configs$n+configs$burn),
                  ystart = configs$y.start,
                  xstart = configs$xreg.start,
                  mut = numeric(configs$n+configs$burn),
                  etat = numeric(configs$n+configs$burn),
                  error = numeric(configs$n+configs$burn),
                  escale = configs$error.scale,
                  Ts = numeric(configs$n+configs$burn),
                  ns = length(configs$seed),
                  seed = configs$seed,
                  rngtype = configs$rngtype,
                  rev = 1L)

  if(out$rev == 1){
    warning("Revision Required. Try changing the link functions\n", immediate. = TRUE)
    return(invisible(out))
  }

  ##------------------------------------
  ## getting the final time series
  ##------------------------------------
  if(configs$burn == 0) u0.star <- u0
  else u0.star <- out$Ts[configs$burn+1]

  ##-----------------------------------------------
  ## if complete = TRUE returns the fulll model.
  ## otherwise only yt is returned
  ##-----------------------------------------------
  ini <- configs$burn + 1
  end <- configs$burn + configs$n
  if(complete){
    final <- list(model = "BARC",
                  yt = out$yt[ini:end],
                  mut = out$mut[ini:end],
                  u0 = u0.star,
                  Ts = out$Ts[ini:end],
                  etat = out$etat[ini:end],
                  error = out$error[ini:end],
                  xreg = out$xreg[ini:end,])
    if(out$nreg == 0) final$xreg <- NULL
    if(debug) final$out.Fortran <- out
  }
  else final <- out$yt[ini:end]

  invisible(final)
}


##' @rdname BARC.functions
##' @order 3
##'
##' @details
##'
##' The function \code{BARC.extract} allows the user to extract the
##' components \eqn{y_t}, \eqn{\mu_t},  \eqn{\eta_t = g(\mu_t)}, \eqn{r_t},
##' \eqn{T^t(u_0)}, the log-likelihood, and the vectors and matrices used to
##' calculate the score vector and the information matrix associated to a given
##' set of parameters.
##'
##' This function can be used by any user to create an objective function
##' that can be passed to optimization functions not available in BTSR Package.
##' At this point, there is no other use for which this function was intended.
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
##' @param r a non-negative integer. The size of the vector theta.
##' If missing, the value of \code{t} is calculated from length(coefs$theta)
##' and length(fixed.values$theta). For fitting, the default is 1.
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
##' @return
##' The function \code{BARC.extract} returns a list with the following components.
##'
##' \itemize{
##'   \item \code{model}: string with the text \code{"BARC"}.
##'
##'   \item \code{coefs}: the coefficients of the model passed through the
##'   \code{coefs} argument.
##'
##'   \item \code{yt}: the observed time series.
##'
##'   \item \code{gyt}: the transformed time series \eqn{g_2(y_t)}.
##'
##'   \item \code{mut}: the conditional mean.
##'
##'   \item \code{etat}: the linear predictor \eqn{g_1(\mu_t)}.
##'
##'   \item \code{error}: the error term \eqn{r_t}.
##'
##'   \item \code{xreg}: the regressors (if included in the model).
##'
##'   \item \code{TS}: the chaotic process \eqn{T^t(u0)}.
##'
##'   \item \code{sll}: the sum of the conditional log-likelihood (if requested).
##'
##'   \item \code{sco}: the score vector  (if requested).
##'
##'   \item \code{info}: the information matrix  (if requested).
##'
##'   \item \code{Drho}, \code{T}, \code{E}, \code{h}: additional matrices and vectors
##'   used to calculate the score vector and the information matrix.  (if requested).
##'
##'   \item \code{yt.new}: the out-of-sample forecast  (if requested).
##'
##'   \item \code{Ts.new}: the out-of-sample forecast for the chaotic
##'   process (if requested).
##'
##'   \item \code{out.Fortran}: FORTRAN output  (if requested).
##' }
##'
##' @seealso
##' \code{\link{btsr.extract}}
##'
##' @examples
##'  #------------------------------------------------------------
##'  # Generating a sample from a BARC model
##'  #------------------------------------------------------------
##'
##'  m1 <- BARC.sim(linkg = "linear", linkh = "linear",
##'                n = 100, seed = 2021, complete = TRUE, ctt.h = 0.6,
##'                coefs = list(nu = 15, theta = 0.85, u0 = pi/4))
##'
##'  #------------------------------------------------------------
##'  #  Extracting the conditional time series given yt and
##'  #  a set of parameters
##'  #------------------------------------------------------------
##'
##'   e1 = BARC.extract(yt = m1$yt, map = 4, ctt.h = 0.6,
##'                     coefs = list(nu = 15, theta = 0.85),
##'                     fixed.values = list(u0 = pi/4),
##'                     linkg = "linear", linkh = "linear", llk = TRUE,
##'                     sco = TRUE, info = TRUE)
##'
##'  #----------------------------------------------------
##'  # comparing the simulated and the extracted values
##'  #----------------------------------------------------
##'  cbind(head(m1$mut), head(e1$mut))
##'
##'  #---------------------------------------------------------
##'  # the log-likelihood, score vector and information matrix
##'  # score vector and information matrix are obtained
##'  # numerically.
##'  #---------------------------------------------------------
##'  e1$sll
##'  e1$score
##'  e1$info.Matrix
##'
##' @export
##' @md
BARC.extract <- function(yt, xreg = NULL, nnew = 0, xnew  = NULL,
                         p, r, coefs = list(), lags = list(),
                         fixed.values = list(), fixed.lags = list(),
                         y.start = NULL, xreg.start = NULL,
                         xregar = TRUE,  error.scale = 0, map = 4,
                         linkg = c("linear","linear"), linkh = "linear",
                         ctt.h = 1, llk = TRUE, sco = FALSE,
                         info = FALSE, debug = FALSE){

  if(is.null(coefs) & is.null(fixed.values))
    stop("Please, provide a list of coefficients")
  if(!is.null(coefs)){
    if(! "list" %in% class(coefs)) stop("coefs must be a list")}

  if(!is.null(fixed.values)){
    if(! "list" %in% class(fixed.values)) stop("fixed.values must be a list")}
  else{ fixed.values <- list()}

  ##----------------------------------------------------
  ## checking if the required parameters are present
  ##----------------------------------------------------
  if(is.null(coefs$u0)){
    if(is.null(fixed.values$u0)) stop("u0 is missing with no default")}

  if(is.null(coefs$theta)) theta <- fixed.values$theta
  else theta <- coefs$theta

  ##--------------------------------------------
  ## checking the map, theta and the link "h"
  ##--------------------------------------------
  cb <- .barc.configs(map = map, theta = theta, linkh = linkh)

  if(missing(p)) p = length(coefs$phi) + length(fixed.values$phi)
  if(missing(r)) r = cb$r

  ##----------------------------------------------------------------------
  ## theoretical score vector and information matrix are not implemented
  ## yet so there is no extra information to be extracted from FORTRAN
  ##----------------------------------------------------------------------
  cf <-  .extract.configs(model = "BARC", yt = yt, y.start = y.start,
                          y.lower = 0, y.upper = 1, openIC = c(TRUE, TRUE),
                          xreg = xreg, xnew = xnew, nnew = nnew,
                          xreg.start = xreg.start, linkg = linkg,
                          p = p, q = r, inf = 0, m = 0, xregar = xregar,
                          error.scale = error.scale, coefs = coefs,
                          lags = lags, fixed.values = fixed.values,
                          fixed.lags = fixed.lags, llk = llk, sco = sco,
                          info = info, extra = FALSE)
  # fixing dummy argument
  if(cf$npar == 0) cf$coefs <- NULL
  #----------------------------------------------------------------
  # merging the information about the map and the linkh with the
  # other configurations
  #----------------------------------------------------------------
  cb <- cb[ "theta" != names(cb)]
  cf[names(cb)] <- cb
  cf$u0 <- .coefs.convert(parname = "u0", fvalues = fixed.values$u0, flags = NULL,
                          coefs = coefs$u0, lags = NULL, npar = 1)
  if(cf$u0$nfix == 0){
    cf$coefs <- c(cf$coefs, u0 = cf$u0$coefs)
    cf$coefsname <- c(cf$coefsname, "u0")
  }
  cf$npar <- as.integer(length(cf$coefs))
  cf$ctt.h <- ctt.h

  # fixing dummy argument
  if(cf$npar == 0) cf$coefs <- 0

  out <- .barc.extract(yt = yt, configs = cf, debug = debug)
  out$model <- "BARC"
  class(out) <- c(class(out), "barc")

  invisible(out)
}


##------------------------------------------------------------------------------------
## internal function:  makes the calculations and reports only the relevant variables
##------------------------------------------------------------------------------------
.barc.extract <- function(yt, configs, debug){

  temp <- .Fortran("barcR",
                   n = configs$n,
                   yt = yt,
                   gyt = numeric(configs$n),
                   ystart = configs$y.start,
                   nreg = configs$nreg,
                   xreg = configs$xreg,
                   xstart = configs$xreg.start,
                   mut = numeric(configs$n),
                   etat = numeric(configs$n),
                   error = numeric(configs$n),
                   escale = configs$error.scale,
                   Ts = numeric(configs$n),
                   nnew = configs$nnew,
                   xnew = configs$xnew,
                   ynew = numeric(max(1,configs$nnew)),
                   Tnew = numeric(max(1,configs$nnew)),
                   link = c(configs$linkg, configs$linkh),
                   ctt.h = configs$ctt.h,
                   map = configs$map,
                   npar = max(1L, configs$npar),
                   coefs = configs$coefs,
                   fixa = configs$alpha$nfix,
                   alpha = configs$alpha$fvalues,
                   fixb = configs$beta$nfix,
                   flagsb = configs$beta$flags,
                   beta = configs$beta$fvalues,
                   p = configs$p,
                   fixphi = configs$phi$nfix,
                   flagsphi = configs$phi$flags,
                   phi = configs$phi$fvalues,
                   xregar = configs$xregar,
                   r = configs$r,
                   fixtheta = configs$theta$nfix,
                   flagstheta = configs$theta$flags,
                   theta = configs$theta$fvalues,
                   fixnu = configs$nu$nfix,
                   nu = configs$nu$fvalues,
                   fixu0 = configs$u0$nfix,
                   u0 = configs$u0$fvalues,
                   llk = configs$llk,
                   sll = 0,
                   sco = configs$sco,
                   U = numeric(max(1, configs$npar*configs$sco)),
                   info = configs$info,
                   K = diag(max(1,configs$npar*configs$info)))

  out <- list(model = "BARC")
  vars <- c("coefs","yt", "xreg", "Ts", "gyt", "mut", "etat", "error")
  out[vars] <- temp[vars]
  if(configs$nreg == 0) out$xreg = NULL

  if(configs$llk == 1) out$sll <- temp$sll
  if(configs$sco == 1){
    out$score <- temp$U
    names(out$score) <- names(configs$coefs)
  }
  if(configs$info == 1){
    out$info.Matrix <- as.matrix(temp$K)
    colnames(out$info.Matrix) <- names(configs$coefs)
    rownames(out$info.Matrix) <- names(configs$coefs)
  }

  if(configs$nnew > 0){
    out$yt.new <- temp$ynew
    out$Ts.new <- temp$Tnew
  }
  if(debug) out$out.Fortran <- temp
  invisible(out)
}

##' @rdname BARC.functions
##' @order 4
##'
##' @details
##' The function \code{BARC.fit} fits a BARC model to a given univariate time
##' series. For now, available optimization algorithms are \code{"L-BFGS-B"} and
##' \code{"Nelder-Mead"}. Both methods accept bounds for the parameters. For
##' \code{"Nelder-Mead"}, bounds are set via parameter transformation.
##'
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
##'   \item \code{model}: string with the text \code{"BARC"}
##'
##'   \item \code{convergence}: An integer code. 0 indicates successful completion.
##'   The error codes depend on the algorithm used.
##'
##'   \item \code{message}: A character string giving any additional information
##'   returned by the optimizer, or NULL.
##'
##'   \item \code{counts}: an integer giving the number of function evaluations.
##'
##'   \item \code{control}: a list of control parameters.
##'
##'   \item \code{start}: the starting values used by the algorithm.
##'
##'   \item \code{coefficients}: 	The best set of parameters found.
##'
##'   \item \code{n}: the sample size used for estimation.
##'
##'   \item \code{series}: the observed time series
##'
##'   \item \code{gyt}: the transformed time series \eqn{g_2(y_t)}
##'
##'   \item \code{fitted.values}: the conditional mean, which corresponds to
##'   the in-sample forecast, also denoted fitted values
##'
##'   \item \code{etat}: the linear predictor \eqn{g_1(\mu_t)}
##'
##'   \item \code{error.scale}: the scale for the error term.
##'
##'   \item \code{error}: the error term \eqn{r_t}
##'
##'   \item \code{residual}: the observed minus the fitted values. The same as
##'   the \code{error} term if \code{error.scale = 0}.
##'
##'   \item \code{forecast}: the out-of-sample forecast for \eqn{y_t} (if requested).
##'
##'   \item \code{Ts.forecas}: the out-of-sample forecast for \eqn{T^t(u_0)}
##'   (if requested).
##'
##'   \item \code{xnew}: the observations of the regressors observed/predicted
##'   values corresponding to the period of out-of-sample forecast.
##'   Only inlcudes if \code{xreg} is not \code{NULL} and \code{nnew > 0}.
##'
##'   \item \code{sll}: the sum of the conditional log-likelihood (if requested)
##'
##'   \item \code{info.Matrix}: the information matrix  (if requested)
##'
##'   \item \code{configs}: a list with the configurations adopted to fit the model.
##'   This information is used by the prediction function.
##'
##'   \item \code{out.Fortran}: FORTRAN output  (if requested)
##'
##'   \item \code{call}: a string with the description of the fitted model.
##'
##' }
##'
##' @seealso
##' \code{\link{btsr.fit}}
##'
##' @examples
##'
##'  #------------------------------------------------------------
##'  # Generating a sample from a BARC model
##'  #------------------------------------------------------------
##'
##'  m1 <- BARC.sim(linkg = "linear", linkh = "linear",
##'                n = 100, seed = 2021, complete = TRUE, ctt.h = 0.6,
##'                coefs = list(nu = 15, theta = 0.85, u0 = pi/4))
##'
##'  #------------------------------------------------------------
##'  #  Fitting a BARC model. Assuming only alpha fixed.
##'  #------------------------------------------------------------
##'   f1 = BARC.fit(yt = m1$yt, map = 4, ctt.h = 0.6,
##'                 start = list(nu = 10, theta = 0.6, u0 = 0.5),
##'                 lower = list(nu = 0, theta = 0, u0 = 0),
##'                 upper = list(theta = 1, u0 = 1),
##'                 fixed.values = list(alpha = 0),
##'                 control = list(iprint = -1, method = "Nelder-Mead"))
##'
##'   coefficients(f1)
##'
##'   plot.ts(m1$yt)
##'   lines(f1$fitted.values, col = "red")
##'
##'  #------------------------------------------------------------
##'  #  Out-of-sample forecast
##'  #------------------------------------------------------------
##'  pred = predict(f1, nnew = 5)
##'  pred$forecast
##'  pred$Ts.forecast
##'
##' @export
##'
##' @md
BARC.fit <- function(yt, xreg = NULL, nnew = 0, xnew = NULL,
                     p = 0, r = 1, start = list(), lags = list(),
                     fixed.values = list(), ignore.start = FALSE,
                     fixed.lags = list(), lower = list(nu = 0, u0 = 0),
                     upper = list(nu = Inf, u0 = 1), map = 4,
                     linkg = c("linear","linear"), linkh = "linear",
                     ctt.h = 1, sco = FALSE, info = FALSE, xregar = TRUE,
                     y.start = NULL, xreg.start = NULL,
                     error.scale = 0, control = list(), report = TRUE,
                     debug = FALSE,...){

  if(report) info = TRUE

  ##--------------------------------------------------------------
  ##  checking if u0, theta, the map and the corresponding link
  ##  Here we need to take into account that both start and
  ##  fixed.values can be empty if start was not initialized yet.
  ##
  ##  if initialization is required the value theta.barc will
  ##  replace the starting value provided by coefs.start
  ##--------------------------------------------------------------
  start0 <- fv <- list()
  start0[names(start)] <- start
  fv[names(fixed.values)] <- fixed.values
  # initilization required for u0?
  if(is.null(c(start0$u0, fv$u0))) start0$u0 = pi/4  # u0 was not provided
  # initilization required for theta?
  if(is.null(start0$theta)) theta.barc <- fv$theta            # theta is fixed
  else theta.barc <- start0$theta                             # checking if theta was given
  if(is.null(theta.barc)) theta.barc <- .theta.start.barc(map) # theta needs initialization


  # default values for nu (merge with user provided values)
  lw <- list(nu = 0, u0 = 0); up <- list(nu = Inf, u0 = 1)
  lw[names(lower)] <- lower; up[names(upper)] <- upper
  lower <- lw; upper <- up
  # fix the lower and upper values for theta (if needed)
  if(is.null(fv$theta)){
    lu <- .theta.lu.fix(map = map, lower = lower$theta, upper = upper$theta)
    lower$theta = lu$lower
    upper$theta = lu$upper
  }

  cb <- .barc.configs(map = map, theta = theta.barc, linkh = linkh)

  ## ------------------------------------------------------------------------
  ## updating start and fixed.values to pass to the configuration function
  ## ------------------------------------------------------------------------
  start <- start0[names(start0) != "u0"]
  fixed.values <- fv[names(fv) != "u0"]
  cf <- .fit.configs(model = "BARC", yt = yt, y.start = y.start,
                     y.lower = 0, y.upper = 1, openIC = c(TRUE, TRUE),
                     xreg = xreg, xnew = xnew, nnew = nnew,
                     xreg.start = xreg.start, linkg = linkg,
                     p = p, d = FALSE, q = r, inf = 0, m = 0,
                     xregar = xregar, error.scale = error.scale,
                     start = start, ignore.start = ignore.start,
                     lags = lags, fixed.values = fixed.values,
                     fixed.lags = fixed.lags, lower = lower,
                     upper = upper, control = control,
                     sco = sco, info = info, extra = FALSE,
                     theta.barc = theta.barc)
  if(!is.null(cf$conv)) return(invisible(out))
  cb <- cb[ "theta" != names(cb)]
  cf[names(cb)] <- cb

  ##--------
  ## u0
  ##--------
  cf$u0 <- .coefs.convert(parname = "u0", fvalues = fv$u0, lags = NULL,
                          flags = NULL, coefs = start0$u0, npar = 1)
  if(cf$u0$nfix == 0){
    cf$coefs <- c(cf$coefs, u0 = cf$u0$coefs)
    cb <- .bounds.convert(npar = 1, lower = lower$u0, upper = upper$u0)
    cf$lower <- c(cf$lower, u0 = cb$lower)
    cf$upper <- c(cf$upper, u0 = cb$upper)
    cf$nbd <- c(cf$nbd, u0 = cb$nbd)
    cf$coefsname = c(cf$coefsname, "u0")
  }

  cf$npar <- length(cf$coefs)
  cf$ctt.h <- ctt.h

  out <- .barc.fit(yt = yt, configs = cf, debug = debug)
  out$call <- .fit.print(model = "BARC", p = cf$p, q = NULL, d = FALSE, nreg = cf$nreg)
  class(out) <- c(class(out), "barc")
  if(report) print(summary(out))
  invisible(out)
}


##------------------------------------------------------------------------------------
## internal function:  makes the calculations and reports only the relevant variables
##------------------------------------------------------------------------------------
.barc.fit <- function(yt, configs, debug){

  if(configs$control$method == "L-BFGS-B"){
    temp <- .Fortran("optimlbfgsbbarcR",
                     npar = max(1L, configs$npar),
                     coefs = configs$coefs,
                     nbd = configs$nbd,
                     lower = configs$lower,
                     upper = configs$upper,
                     n = configs$n,
                     yt = yt,
                     gy = numeric(configs$n),
                     ystart = configs$y.start,
                     nreg = configs$nreg,
                     xreg = configs$xreg,
                     xstart = configs$xreg.start,
                     mut = numeric(configs$n),
                     etat = numeric(configs$n),
                     error = numeric(configs$n),
                     escale = configs$error.scale,
                     Ts = numeric(configs$n),
                     nnew = configs$nnew,
                     xnew = configs$xnew,
                     ynew = numeric(max(1,configs$nnew)),
                     Tnew = numeric(max(1,configs$nnew)),
                     link = c(configs$linkg, configs$linkh),
                     ctt.h = configs$ctt.h,
                     map = configs$map,
                     fixa = configs$alpha$nfix,
                     alpha = configs$alpha$fvalues,
                     fixb = configs$beta$nfix,
                     flagsb = configs$beta$flags,
                     beta = configs$beta$fvalues,
                     p = configs$p,
                     fixphi = configs$phi$nfix,
                     flagsphi = configs$phi$flags,
                     phi = configs$phi$fvalues,
                     xregar = configs$xregar,
                     r = configs$r,
                     fixtheta = configs$theta$nfix,
                     flagstheta = configs$theta$flags,
                     theta = configs$theta$fvalues,
                     fixnu = configs$nu$nfix,
                     nu = configs$nu$fvalues,
                     fixu0 = configs$u0$nfix,
                     u0 = configs$u0$fvalues,
                     sll = 0,
                     U = numeric(max(1, configs$npar)),
                     info = configs$info,
                     K = diag(max(1, configs$npar*configs$info)),
                     iprint = as.integer(configs$control$iprint),
                     factr = configs$control$factr,
                     pgtol = configs$control$pgtol,
                     maxit = as.integer(configs$control$maxit),
                     neval = 0L,
                     conv = 0L)
  }else{
    temp <- .Fortran("optimnelderbarcR",
                     npar = max(1L, configs$npar),
                     coefs = configs$coefs,
                     nbd = configs$nbd,
                     lower = configs$lower,
                     upper = configs$upper,
                     n = configs$n,
                     yt = yt,
                     gy = numeric(configs$n),
                     ystart = configs$y.start,
                     nreg = configs$nreg,
                     xreg = configs$xreg,
                     xstart = configs$xreg.start,
                     mut = numeric(configs$n),
                     etat = numeric(configs$n),
                     error = numeric(configs$n),
                     escale = configs$error.scale,
                     Ts = numeric(configs$n),
                     nnew = configs$nnew,
                     xnew = configs$xnew,
                     ynew = numeric(max(1,configs$nnew)),
                     Tnew = numeric(max(1,configs$nnew)),
                     link = c(configs$linkg, configs$linkh),
                     ctt.h = configs$ctt.h,
                     map = configs$map,
                     fixa = configs$alpha$nfix,
                     alpha = configs$alpha$fvalues,
                     fixb = configs$beta$nfix,
                     flagsb = configs$beta$flags,
                     beta = configs$beta$fvalues,
                     p = configs$p,
                     fixphi = configs$phi$nfix,
                     flagsphi = configs$phi$flags,
                     phi = configs$phi$fvalues,
                     xregar = configs$xregar,
                     r = configs$r,
                     fixtheta = configs$theta$nfix,
                     flagstheta = configs$theta$flags,
                     theta = configs$theta$fvalues,
                     fixnu = configs$nu$nfix,
                     nu = configs$nu$fvalues,
                     fixu0 = configs$u0$nfix,
                     u0 = configs$u0$fvalues,
                     sll = 0,
                     sco = configs$sco,
                     U = numeric(max(1,configs$npar*configs$sco)),
                     info = configs$info,
                     K = diag(max(1, configs$npar*configs$info)),
                     iprint = as.integer(configs$control$iprint),
                     stopcr = configs$control$stopcr,
                     maxit = as.integer(configs$control$maxit),
                     neval = 0L,
                     conv = 0L)
  }

  temp$llk <- 1
  temp$sco <- configs$sco
  # for some reason, sometimes info returns NULL
  # from Fortran. We need to fix this!!
  temp$info <- configs$info

  out <- .fit.get.results(model = "BARC", temp, configs = configs)
  if(debug) out$out.Fortran <- temp
  invisible(out)

}

##------------------------------------------------------------------------------------
## internal function:  makes the calculations and reports only the relevant variables
##------------------------------------------------------------------------------------
.barc.predict <- function(object, debug){

  if(object$model != "BARC") stop("Wrong configurations for BARC models")

  temp <-  .Fortran("predictbarcR",
                    n = object$n,
                    series = object$series,
                    gyt = object$gyt,
                    nreg = object$nreg,
                    xreg = object$xreg,
                    escale = object$error.scale,
                    error = object$error,
                    Ts = object$Ts,
                    nnew = object$nnew,
                    xnew = object$xnew,
                    ynew = numeric(max(1,object$nnew)),
                    Tnew = numeric(max(1,object$nnew)),
                    link = c(object$linkg, object$linkh),
                    ctt.h = object$ctt.h,
                    map = object$map,
                    npar = max(1L, object$npar),
                    coefs = object$coefs,
                    fixa = object$alpha$nfix,
                    alpha = object$alpha$fvalues,
                    fixb = object$beta$nfix,
                    flagsb = object$beta$flags,
                    beta = object$beta$fvalues,
                    p = object$p,
                    fixphi = object$phi$nfix,
                    flagsphi = object$phi$flags,
                    phi = object$phi$fvalues,
                    xregar = object$xregar,
                    r = object$r,
                    fixtheta = object$theta$nfix,
                    flagstheta = object$theta$flags,
                    theta = object$theta$fvalues,
                    fixnu = object$nu$nfix,
                    nu = object$nu$fvalues,
                    fixu0 = object$u0$nfix,
                    u0 = object$u0$fvalues)

  out <- list(model = object$model,
              yt.new = temp$ynew,
              Ts.new = temp$Tnew)
  if(debug) out$out.Fortran <- temp
  invisible(out)
}


##-------------------------------------------------------------------------
## Internal function: Used to check if theta is compatible with the map
## selected by the user
##-------------------------------------------------------------------------
.check.map <- function(map, theta){
  ##------------------------------------------
  ##           Maps
  ##------------------------------------------
  ##  1 = (kx)(mod 1).  k integer
  ##  2 = Rafael's map.  0 <= theta <= 1
  ##  3 = logistic map. 0 <= theta  <= 4
  ##  4 = Manneville-Pomeau. 0 < theta < 1
  ##  5 = Lasota-Mackey's map. No theta
  ##-------------------------------------------

  if((map != 5) & is.null(theta))
    stop("theta is missing with no default")
  maps <- data.frame(map = c(1:5),
                     r = c(1,1,1,1,0),
                     lower = c(1,0,0,0,NA),
                     upper = c(Inf,1,4,1,NA))
  r <- maps[map, "r"]
  if(length(theta) > r & r == 0)
    print(paste(msg, "Theta will be ignored.", sep = " "))
  if(r == 0) return(invisible(list(theta = 0, r = 0L)))

  ## checking the length of theta
  msg <- NULL
  if(length(theta) != r)
    msg <- paste(msg,"Length of theta = ", length(theta),
                 " but the map requires length ", r, ".", sep = "")
  if(length(theta) > r & r > 0)
    msg <- paste(msg, "Only the first ", r, " values will be used.", sep = "")

  if(length(theta) < r)
    stop("Please provide a new theta with length ", r)

  ## checking the range of theta
  ## !!! if theta is a vector this part of the code will need revision !!!!!
  if(theta < maps[map, "lower"] |  theta > maps[map, "upper"])
    stop("Theta is out of range")

  if(!is.null(msg)) print(msg)
  invisible(list(map = as.integer(map), theta = theta[1:r], r = as.integer(r)))
}


##-------------------------------------------------------------------------
## Internal function: Used to check the configurations for simulation
## of BARC models
##-------------------------------------------------------------------------
.barc.configs <- function(map = 4, theta = 0.5, linkh = "linear"){

  out <- c()

  ## checking the map and the corresponding parameter
  out <- .check.map(map = map, theta = theta)

  ## link function
  out$linkh <- .link.convert(linkh)
  if(is.na(out$linkh)) stop(paste("link ", linkh, " not implemented", sep = ""))

  return(out)

}

##-------------------------------------------------------------------------
## Internal function: Used to check if theta is compatible with the map
## selected by the user
##-------------------------------------------------------------------------
.theta.start.barc <- function(map){
  ##------------------------------------------
  ##           Maps
  ##------------------------------------------
  ##  1 = (kx)(mod 1).  k integer
  ##  2 = Rafael's map.  0 <= theta <= 1
  ##  3 = logistic map. 0 <= theta  <= 4
  ##  4 = Manneville-Pomeau. 0 < theta < 1
  ##  5 = Lasota-Mackey's map. No theta
  ##-------------------------------------------
  theta <- c(3, 0.5, 3.5, 0.5, 0)
  theta[map]
}

##-------------------------------------------------------------------------
## Internal function: Used to check if the limits for theta are correct
##-------------------------------------------------------------------------
.theta.lu.fix <- function(map, lower, upper){
  maps <- data.frame(map = c(1:5),
                     lower = c(1,0,0,0,NA),
                     upper = c(Inf,1,4,1,NA))

  # checking if lower and upper limits are provided
  fix.lower <- fix.upper <- FALSE
  if(is.null(lower)) fix.lower <- TRUE
  if(is.null(upper)) fix.upper <- TRUE

  # checking if lower and upper are in the correct range
  if(!is.null(lower))
    if(lower < maps[map, "lower"] | lower > maps[map, "upper"])
      fix.lower <- TRUE
  if(!is.null(upper))
    if(upper < maps[map, "lower"] | upper > maps[map, "upper"])
      fix.upper <- TRUE

  # if needed, fix the wrong values
  if(fix.lower) lower <- maps[map, "lower"]
  if(fix.upper) upper <- maps[map, "upper"]

  if(lower > upper)
    stop("Please, check the lower and upper limits for theta")

  if(lower == upper & map != 5){
    msg = "lower and upper limits for theta are the same.\n "
    msg = paste0(msg, "The range for the selected map is,\n")
    warning(paste0(msg, "lower = ", maps[map, "lower"],"\n",
                   "upper = ", maps[map, "upper"]))
  }
  invisible(list(lower = lower, upper = upper))
}
