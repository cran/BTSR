#---------------------   BARC MODELS    ----------------------------------------------
#' @title
#' Functions to simulate, extract components and fit BARC models
#'
#' @name BARC.functions
#' @order 1
#'
#' @description
#' These functions can be used to simulate, extract components and fit any model
#' of the class `barc`. A model with class `barc` is a special case of a model
#' with class `btsr`. See the Section \sQuote{The BTSR structure} in
#' [btsr-package] for more details on the general structure. See also
#' \sQuote{Details} below.
#'
#' @param n a strictly positive integer. The sample size of `yt` (after
#'   burn-in). Default is `n = 1`.
#'
#' @param burn a non-negative integer. length of "burn-in" period. Default is
#'   `burn = 0`.
#'
#' @param nnew optionally, the number of out-of sample predicted values
#'   required. Default is `nnew = 0`.
#'
#' @param yt a numeric vector with the observed time series. If missing, an
#'   error message is issued.
#'
#' @param y.start optionally, an initial value for  \eqn{Y_t} (to be used in the
#'   recursions). Default is `y.start = NULL`, in which case, the recursion
#'   assumes that \eqn{Y_t = g_{12}^{-1}(0)}, for \eqn{t < 1}. Only relevant if
#'   \eqn{p > 0}.
#'
#' @param xreg optionally, a vector or matrix of external regressors. Default
#'   is `xreg = NULL`. For simulation purposes, the length of `xreg` must be
#'   equal to `n + burn`. For extraction or fitting purposes, the length of
#'   `xreg` must be the same as the length of the observed time series
#'   \eqn{Y_t}.
#'
#' @param xnew  a vector or matrix, with `nnew` observations of the regressors
#'   observed/predicted values corresponding to the period of out-of-sample
#'   forecast. Default is `xreg = NULL`. If `xreg = NULL` or `nnew = 0`, `xnew`
#'   is ignored. If `nnew > 0` and the number of regressors in `xnew` does not
#'   match `xreg` an error message is issued.
#'
#' @param xreg.start optionally, a vector of initial value for
#'   \eqn{\boldsymbol{X}_t} (to be used in the recursions). Default is
#'   `xreg.start = NULL`, in which case, the average of the first \eqn{p} values
#'   is used, that is, the recursion assumes that  \eqn{\boldsymbol{X}_t =
#'   p^{-1}\sum_{k = 1}^p  \boldsymbol{X}_k}, for \eqn{t < 1}. Only relevant if
#'   `xregar = TRUE` and  \eqn{p > 0}.
#'
#' @param xregar logical; indicates whether `xreg` should be included in the AR
#'   recursion of the model. Default is `xregar = TRUE`. Only relevant if `xreg`
#'   is included and \eqn{p > 0}. See the Section \sQuote{The BTSR structure} in
#'   [btsr-package] for details.
#'
#' @param p optionally, a non-negative integer. The order of the AR polynomial.
#'   Default is `p = NULL`, in which case the value of `p` is computed
#'   internally, based on the size of the argument `phi` in the lists of
#'   coefficients (or staring values), fixed lags, and fixed values. For fitting
#'   purposes, if `p` and `start` are both `NULL`, an error message is issued.
#'
#' @param ignore.start logical; indicates whether the argument `start` should be
#'   ignored. If starting values are not provided, the function uses the default
#'   values and `ignore.start` is ignored. In case starting values are provided
#'   and `ignore.start = TRUE`, those starting values are ignored and
#'   recalculated. The default is `ignore.start = FALSE`.
#'
#' @param start a list with the starting values for the non-fixed coefficients
#'   of the model. The default is `start = NULL`, in which case the function
#'   [coefs.start] is used internally to obtain starting values for the
#'   parameters.
#'
#' @param coefs a list with the coefficients of the model. An empty list will
#'   result in an error. The arguments that can be passed through this list are
#'   \itemize{
#'    \item `alpha`: optionally, a numeric value corresponding to the intercept.
#'     If the argument is missing, it will be treated as zero.
#'
#'    \item `beta`: optionally, a vector of coefficients corresponding to the
#'     regressors in `xreg`. For simulation purposes, if `xreg` is provided but
#'     `coefs` does not have a `beta` argument, an error message is issued. The
#'     extracting function also verify the `fixed.values` list before issuing an
#'     error message.
#'
#'    \item`phi`: optionally, for the simulation function this must be a vector
#'    of size, \eqn{p} corresponding to the autoregressive coefficients
#'    (including the ones that are zero), where \eqn{p} is the AR order. For the
#'    extraction and fitting functions, this is a vector with the non-fixed
#'    values in the vector of autoregressive coefficients.
#'
#'    \item `theta` the parameter (or vector of parameters) corresponding to the
#'    map function. If `map = 5` this value is ignored. For simulation,
#'    purposes, the default is `map = 4` and `theta = 0.5`. \strong{Note:} Do
#'    not confuse `theta` from a BARC model with the moving average term in the
#'    general BTSR class of models
#'
#'    \item `nu` the dispersion parameter. If missing, an error message is
#'    issued.
#'
#'    \item `u0` a numeric value in the interval \eqn{(0,1)}, corresponding to
#'    the value of the random variable \eqn{U_0}. For simulation purposes, the
#'    default is `u0 = pi/4`.
#'   }
#'   For simulation purposes, an empty list will result in an error message. For
#'   extraction purposes, an error message will be issued if both `coefs` and
#'   `fixed.values` are empty. The argument `coefs` is not used when fitting a
#'   model. Missing parameters are treated as zero.
#'
#' @param lags optionally, a list with the lags that the values in `coefs`
#'   correspond to. The names of the entries in this list must match the ones in
#'   `coefs` (or `start`). For one dimensional coefficients, the `lag` is
#'   obviously always 1 and can be suppressed. The default is `lags = NULL`, in
#'   which the `lags` are computed from the `fixed.lags` argument (if provided).
#'   If both, `lags` and `fixed.lags` are missing, it is assumed that all lags
#'   must be used. The arguments `lags` and `fixed.lags` are complementary.
#'   Either suffices, or mix them (e.g., `lags` for some parameters,
#'   `fixed.lags` for others).
#'
#' @param fixed.values optionally, a list with the values of the coefficients
#'   that are fixed. The default is `fixed.lags = NULL`. By default, if a given
#'   vector (such as the vector of AR coefficients) has fixed values and the
#'   corresponding entry in this list is empty, the fixed values are set as
#'   zero. The names of the entries in this list must match the ones in `coefs`
#'   (or `start`).
#'
#' @param fixed.lags optionally, a list with the lags that the fixed values in
#'   `fixed.values` correspond to. The names of the entries in this list must
#'   match the ones in `fixed.values`. For one dimensional coefficients, the
#'   `lag` is obviously always 1 and can be suppressed. If an empty list is
#'   provided and the model has fixed lags, the argument `lags` is used as
#'   reference.
#'
#' @param lower  optionally, list with the lower bounds for the parameters. The
#'   names of the entries in these lists must match the ones in `start`. Default
#'   is `lower = NULL`. The default is to assume that the parameters have no
#'   lower bound except for `nu`, for which de default is 0. Only the bounds for
#'   bounded parameters need to be specified.
#'
#' @param upper optionally, list with the upper bounds for the parameters. The
#'   names of the entries in these lists must match the ones in `start`. Default
#'   is `upper = NULL`.  The default is to assume that the parameters have no
#'   upper bound. Only the bounds for bounded parameters need to be specified.
#'
#' @param error.scale the scale for the error term. Default is `error.scale = 0`
#'   (data scale).
#'
#' @param linkg character or a two character vector indicating which links must
#'   be used in the model.  See the Section \sQuote{The BTSR structure} in
#'   [btsr-package] for details and [link.btsr] for valid links. If only one
#'   value is provided, the same link is used for \eqn{\mu_t} and for \eqn{Y_t}
#'   in the AR part of the model.  Default is `linkg = "linear"`.
#'
#' @param configs.linkg a list with two elements, `ctt` and `power`, which
#'   define the constant \eqn{a} and the exponent \eqn{b} in the link function
#'   \eqn{g(x) = a x^b}. Each element can be a single numeric value (applied
#'   uniformly across all linear links), a numeric vector of length 2, or a
#'   named list with entries `g11` and `g12`. This argument is only used when
#'   the link function is `"linear"` or `"polynomial"`. The default is
#'   `configs.linkg = NULL`, in which case the function internally assumes
#'   `configs.linkg = list(ctt = list(g11 = 1, g12 = 1), power = list(g11 = 1,
#'   g12 = 1))`.
#'
#' @param linkh a character indicating which link must be associated to the
#'   chaotic process.  See the Section \sQuote{The BTSR structure} in
#'   [btsr-package] for details and [link.btsr] for valid links. Default is
#'   `linkh = "linear"`.
#'
#' @param configs.linkh a list with extra configurations for the link \eqn{h}.
#'   For now, only used if `linkh = "linear"` or `"polynomial"`. Default is
#'   `configs.linkh = list(ctt = 1, power = 1)`.
#'
#' @param complete logical; if FALSE returns only `yt`, else returns additional
#'   components. Default is `complete = FALSE`.
#'
#' @param llk logical; indicates whether the value of the log-likelihood
#'   function should be returned. Default is `llk = TRUE`.
#'
#' @param sco logical; indicates whether the score vector should be returned.
#'   Default is `sco = FALSE`. For now, the score vector is computed using
#'   numerical derivatives.
#'
#' @param info logical; indicates whether the information matrix should be
#'   returned. Default is `info = FALSE`. For the fitting function, `info` is
#'   automatically set to `TRUE` when `report = TRUE`.  For now, the information
#'   matrix is computed using numerical derivatives.
#'
#' @param control a list with configurations to be passed to the optimization
#'   subroutines. Default is `control = NULL`. Missing arguments will receive
#'   default values. For details, see [fit.control].
#'
#' @param report logical; indicates whether the summary from the fitted model
#'   should be be printed. Default is `report = TRUE`, in which case `info` is
#'   automatically set to `TRUE`.
#'
#' @param debug logical, if `TRUE` the output from Fortran is return (for
#'   debugging purposes). Default is `debug = FALSE`.
#'
#' @param ... further arguments passed to the internal functions. See, for
#'   instance, [summary.btsr] for details.
#'
#' @details
#'
#' ## Sim, Extract and Fit functions
#'
#' The function \code{BARC.sim} generates a random sample from a BARC(\eqn{p})
#' model.
#'
#' The function `BARC.extract` allows the user to extract the components
#' \eqn{Y_t}, \eqn{\mu_t},  \eqn{\eta_t = g(\mu_t)}, \eqn{e_t},  \eqn{T^t(U_0)},
#' the log-likelihood, the score vector and the information matrix associated to
#' a given set of parameters. This function can be used by any user to create an
#' objective function that can be passed to optimization algorithms not
#' available in the BTSR Package.
#'
#' The function `BARC.fit` fits a BARC model to a given univariate time series.
#' For now, available optimization algorithms are `"L-BFGS-B"` and
#' `"Nelder-Mead"`. Both methods accept bounds for the parameters. For
#' `"Nelder-Mead"`, bounds are set via parameter transformation.
#'
#' ## Particular cases
#'
#' Neither the beta regression or an i.i.d. sample from a beta distribution can
#' be obtained as special cases of the BARC model since the term \eqn{h(T(U_0))}
#' is always present.
#'
#' The model from \insertCite{pumi2021;textual}{BTSR} is obtained by setting
#' `xregar = TRUE` (so that the regressors are included in the AR part of the
#' model) and using the same link for \eqn{Y_t} and \eqn{\mu_t}.
#'
#' @template param_map
#'
#' @return
#' By default, the function `BARC.sim` returns the simulated time series `yt`.
#' If `complete = TRUE`, it returns a list with the following components
#' \itemize{
#'  \item `model`: string with the text `"BARC"`
#'
#'  \item `yt`: the simulated time series \eqn{Y_t}
#'
#'  \item `mut`: the conditional mean \eqn{\mu_t}
#'
#'  \item `etat`: the linear predictor \eqn{\eta_t = g_{11}(\mu_t)}
#'
#'  \item `u0`: the starting values of \eqn{U_0}
#'
#'  \item `Ts`: the chaotic process \eqn{T^t(U_0)}
#'
#'  \item `error`: the error term \eqn{e_{1t}}
#'
#'  \item `out.Fortran`: the output from FORTRAN (if requested).
#' }
#'
#' The function `BARC.extract` returns a list with the following components.
#'
#' \itemize{
#'  \item `model`: string with the text `"BARC"`
#'
#'  \item `yt`: the observed time series \eqn{Y_t}
#'
#'  \item `TS`: the chaotic process \eqn{T^t(U_0)}.
#'
#'  \item `mut`: the conditional mean \eqn{\mu_t}
#'
#'  \item `etat`: the linear predictor \eqn{\eta_t = g_{11}(\mu_t)}
#'
#'  \item `error`: the error term \eqn{e_{1t}}
#'
#'  \item `forecast`: the out-of-sample forecast  (if requested)
#'
#'  \item `xnew`: the out-of-sample values of `xreg` provided by the user (only
#'  present if the model includes regressors and forecast is requested)
#'
#'  \item `sll`: the sum of the conditional log-likelihood (if requested)
#'
#'  \item `score`: the score vector  (if requested)
#'
#'  \item `info.Matrix.`: the score vector  (if requested)
#'
#'  \item `out.Fortran`: FORTRAN output  (if requested)
#' }
#'
#' The function `BARC.fit` returns a list with the following components.
#'
#' \itemize{
#'  \item `model`: string with the text `"BARC"`
#'
#'  \item `call`: string with a complete description of the model, including
#'   the AR and MA order.
#'
#'  \item `n`: the sample size used for estimation.
#'
#'  \item `series`: the observed time series \eqn{Y_t}
#'
#'  \item `gyt`: a vector or a matrix with the transformed time series
#'   \eqn{g_{11}(Y_t)} and \eqn{g_{12}(Y_t)}. Only returns a matrix if the links
#'   \eqn{g_{11}} and \eqn{g_{12}} are not the same.
#'
#'  \item `xreg`: a vector or matrix of regressors \eqn{\boldsymbol{X}_t} (if
#'   included in the model).
#'
#'  \item `control`: a list of control parameters.
#'
#'  \item `convergence`: An integer code. 0 indicates successful completion. The
#'   error codes depend on the algorithm used.
#'
#'  \item `message`: A character string giving any additional information
#'   returned by the optimizer (if any), or `NULL`.
#'
#'  \item `counts`: an integer giving the number of function evaluations.
#'
#'  \item `start`: the starting values used by the algorithm.
#'
#'  \item `coefficients`: The best set of parameters found.
#'
#' \item `fitted.values`: the conditional time series \eqn{\mu_t} and the
#'  chaotic process \eqn{T^t(U_0)}, which corresponds to the in-sample forecast,
#'  also denoted fitted values.
#'
#'  \item `etat`: the linear predictor \eqn{\eta_{1t} = g_{11}(\mu_t)}
#'
#'  \item `error`: the error term \eqn{e_{1t}}
#'
#'  \item `residual`: the observed values \eqn{Y_t} minus the fitted values
#'   \eqn{\mu_t}. The same as the `error` term if `error.scale = 0`.
#'
#'  \item `forecast`: a matrix with the out-of-sample forecast (if requested)
#'   for \eqn{\mu_t} and \eqn{\eta_{1t}}
#'
#'  \item `xnew`: the observations of the regressors observed/predicted values
#'   corresponding to the period of out-of-sample forecast. Only included if
#'   `xreg` is not `NULL` and `nnew > 0`.
#'
#'  \item `sll`: the sum of the conditional log-likelihood (if requested)
#'
#'  \item `score`: the score vector (if requested)
#'
#'  \item `info.Matrix`: the information matrix (if requested)
#'
#'  \item `link`: the codes for the link functions (for summary purposes)
#'
#'  \item `configs`: a list with the configurations passed to FORTRAN to fit the
#'   model. This information is used by the prediction function.
#'
#'  \item `out.Fortran`: FORTRAN output  (if requested).
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' [BTSR.functions]: sim, extract and fit functions for BTSR models
#'
#' [BTSR.parent.models]: sim, extract and fit functions for parent models
#'
#' [get.defaults]: Retrieve default arguments for BTSR package functions
#'
#'
#' @examples
#' \donttest{
#'
#' #########################################################################
#' #
#' #   Example of usage of BARC.sim, BARC.extract and BARC.fit
#' #
#' #########################################################################
#'
#' #------------------------------------------------------------
#' # Generating a sample from a BARC model
#' #------------------------------------------------------------
#' set.seed(1234)
#' m1 <- BARC.sim(
#'   coefs = list(nu = 15, theta = 0.85, u0 = pi / 4),
#'   linkg = "linear",
#'   linkh = "linear",
#'   configs.linkh = list(ctt = 0.6),
#'   n = 100,
#'   complete = TRUE
#' )
#'
#' plot.ts(m1$yt)
#' lines(m1$mut, col = "red")
#'
#' #------------------------------------------------------------
#' #  Extracting the conditional time series given yt and
#' #  a set of parameters
#' #------------------------------------------------------------
#'
#' e1 <- BARC.extract(
#'   yt = m1$yt,
#'   map = 4,
#'   coefs = list(nu = 15, theta = 0.85),
#'   fixed.values = list(u0 = pi / 4),
#'   linkg = "linear",
#'   linkh = "linear",
#'   configs.linkh = list(ctt = 0.6),
#'   llk = TRUE,
#'   sco = TRUE,
#'   info = TRUE
#' )
#'
#' #----------------------------------------------------
#' # comparing the simulated and the extracted values
#' #----------------------------------------------------
#' cbind(head(m1$mut), head(e1$mut))
#'
#' #---------------------------------------------------------
#' # the log-likelihood, score vector and information matrix
#' # score vector and information matrix are obtained
#' # numerically.
#' #---------------------------------------------------------
#' e1$sll
#' e1$score
#' e1$info.Matrix
#'
#' #------------------------------------------------------------
#' #  Fitting a BARC model. Assuming only alpha fixed.
#' #------------------------------------------------------------
#' f1 <- BARC.fit(
#'   yt = m1$yt,
#'   map = 4,
#'   configs.linkh = list(ctt = 0.6),
#'   start = list(nu = 10, theta = 0.6, u0 = 0.5),
#'   lower = list(nu = 0, theta = 0, u0 = 0),
#'   upper = list(theta = 1, u0 = 1),
#'   fixed.values = list(alpha = 0),
#'   control = list(iprint = -1, method = "Nelder-Mead")
#' )
#'
#' coefficients(f1)
#'
#' plot.ts(m1$yt)
#' lines(f1$fitted.values[, "mut"], col = "red")
#'
#' #------------------------------------------------------------
#' #  Out-of-sample forecast
#' #------------------------------------------------------------
#' pred <- predict(f1, nnew = 5)
#' pred$forecast
#' }
#'
NULL
#> NULL


#' @rdname BARC.functions
#' @order 2
#' @export
BARC.sim <- function(
    n = 1, burn = 0, y.start = NULL,
    xreg = NULL, xreg.start = NULL, xregar = TRUE,
    coefs = NULL, map = 4, error.scale = 0,
    linkg = "linear", configs.linkg = NULL,
    linkh = "linear", configs.linkh = list(ctt = 1, power = 1),
    complete = FALSE, debug = FALSE) {
  invisible(btsr.sim(
    model = "BARC",
    n = n, burn = burn, xreg = xreg, map = map, coefs = coefs, y.start = y.start,
    xreg.start = xreg.start, xregar = xregar, error.scale = error.scale,
    linkg = linkg, configs.linkg = configs.linkg, linkh = linkh,
    configs.linkh = configs.linkh, complete = complete, debug = debug
  ))
}


#' @rdname BARC.functions
#' @order 3
#' @export
BARC.extract <- function(
    yt, y.start = NULL, xreg = NULL, xreg.start = NULL,
    xnew = NULL, xregar = TRUE, nnew = 0, p = NULL,
    coefs = NULL, lags = NULL, fixed.values = NULL, fixed.lags = NULL,
    error.scale = 0, map = 4,
    linkg = "linear", configs.linkg = NULL,
    linkh = "linear", configs.linkh = list(ctt = 1, power = 1),
    llk = TRUE, sco = FALSE, info = FALSE, debug = FALSE) {
  invisible(btsr.extract(
    model = "BARC",
    yt = yt, xreg = xreg, nnew = nnew, xnew = xnew, p = p, coefs = coefs,
    lags = lags, fixed.values = fixed.values, fixed.lags = fixed.lags,
    y.start = y.start, xreg.start = xreg.start, xregar = xregar,
    error.scale = error.scale, map = map, linkg = linkg,
    configs.linkg = configs.linkg, linkh = linkh, configs.linkh = configs.linkh,
    llk = llk, sco = sco, info = info, debug = debug
  ))
}

#' @rdname BARC.functions
#' @order 4
#' @export
BARC.fit <- function(
    yt, y.start = NULL, xreg = NULL, xreg.start = NULL, xregar = TRUE,
    xnew = NULL, nnew = 0, p = NULL, ignore.start = FALSE,
    start = NULL, lags = NULL, fixed.values = NULL,
    fixed.lags = NULL, lower = NULL, upper = NULL,
    map = 4,
    linkg = "linear", configs.linkg = NULL,
    linkh = "linear", configs.linkh = list(ctt = 1, power = 1),
    sco = FALSE, info = FALSE, error.scale = 0,
    control = NULL, report = TRUE, debug = FALSE, ...) {
  invisible(btsr.fit(
    model = "BARC",
    yt = yt, xreg = xreg, nnew = nnew, xnew = xnew, p = p, d = FALSE,
    m = 0, inf = 0, start = start, ignore.start = ignore.start,
    lags = lags, fixed.values = fixed.values, fixed.lags = fixed.lags,
    lower = lower, upper = upper, y.start = y.start, xreg.start = xreg.start,
    xregar = xregar, error.scale = error.scale, map = map, linkg = linkg,
    configs.linkg = configs.linkg, linkh = linkh, configs.linkh = configs.linkh,
    llk = 1, sco = sco, info = info, extra = FALSE, control = control,
    report = report, debug = debug, ...
  ))
}
