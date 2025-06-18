#---------------  Generic functions: sim, extract and fit ----------------------
#' @title
#' Generic functions to simulate, extract components and fit BTSR models
#'
#' @name BTSR.functions
#' @order 1
#'
#' @description
#' These generic functions can be used to simulate, extract components and fit
#' any model of the class `btsr`. See \sQuote{Details} below. \cr
#' @template section_structure_description
#'
#' @inheritParams arguments
#'
#' @details
#'
#' A detailed description of the general structure (mathematical formulation) of
#' BTSR models, associated to the `btsr` class, is presented in Section
#' \sQuote{The BTSR Structure} of [btsr-package]. Particular models are
#' discussed in [arguments.model].\cr
#'
#' All functions are compatible with the new format for the arguments,
#' introduced in version 1.0.0. and the previous format.
#' \itemize{
#'  \item The function `btsr.sim` is used to generate random samples from any
#'   BTSR models.
#'
#'  \item The function `btsr.extract` allows the user to extract all conditional
#'   time series, the log-likelihood, and the vectors and matrices used to
#'   calculate the score vector and the information matrix associated to a given
#'   set of parameters. This function can be used by any user to create an
#'   objective function that can be passed to optimization functions not available
#'   in BTSR Package. At this point, there is no other use for which this function
#'   was intended.
#'
#'  \item The function `btsr.fit` fits a BTSR model to a given univariate time
#'   series. For now, available optimization algorithms are `"L-BFGS-B"` and
#'   `"Nelder-Mead"`. Both methods accept bounds for the parameters. For
#'   `"Nelder-Mead"`, bounds are set via parameter transformation.
#' }
#'
#' For compatibility with previous versions of the package, all functions
#' associated to parent models (e.g. BARFIMA) are still available (see
#' [BTSR.parent.models]). Also, analogous functions are available for parent
#' models with time varying \eqn{\nu} (e.g. BARFIMAV). The list of arguments and
#' default values for these specific functions can be accessed using the
#' function [get.defaults].\cr
#'
#' Particular models (e.g. BETA, BARMA) share the same arguments as the parent
#' model, however, some arguments can have different default values (see the
#' documentation for [shared arguments][arguments] for details). Information on
#' the parent model can be obtained using the function [BTSR.model.defaults].
#'
#' @template section_return_sim
#' @template section_return_extract
#' @template section_return_fit
#'
#' @seealso
#' [BARC.functions]: sim, extract and fit functions for BARC models
#'
#' [BTSR.parent.models]: sim, extract and fit functions for parent models
#'
#' [get.defaults]: Retrieve default arguments for BTSR package functions
#'
#' @examples
#' \donttest{
#' #########################################################################
#' #
#' #   Examples of usage of btsr.sim
#' #
#' #########################################################################
#' #------------------------------------------------------------------------
#' # Generating a Beta model were both mu and nu do not vary with time
#' # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
#' #------------------------------------------------------------------------
#'
#' # CASE 1: using the legacy format for coefs
#' set.seed(1234)
#' y1 <- btsr.sim(
#'   model = "BETA", n = 1000,
#'   coefs = list(alpha = 0.2, nu = 20)
#' )
#' hist(y1)
#'
#' # CASE 2: using the new layout for coefs
#' set.seed(1234)
#' y2 <- btsr.sim(
#'   model = "BETA", n = 1000,
#'   coefs = list(part1 = list(alpha = 0.2), part2 = list(alpha = 20))
#' )
#' hist(y2)
#'
#' # CASE 3: function for the parent model plus legacy format for coefs.
#' # - requires setting linkg = "linear", otherwhise the default "logit"
#' #   link is used.
#' set.seed(1234)
#' y3 <- BARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.2, nu = 20)
#' )
#' hist(y3)
#'
#' # CASE 4: function for the parent model plus new format for coefs.
#' # - requires setting linkg = "linear", otherwhise the default "logit"
#' #   link is used.
#' set.seed(1234)
#' y4 <- BARFIMA.sim(
#'   n = 1000, linkg = "linear",
#'   coefs = list(part1 = list(alpha = 0.2), part2 = list(alpha = 20))
#' )
#' hist(y4)
#'
#' # comparing the results:
#' range(abs(y2 - y1))
#' range(abs(y3 - y1))
#' range(abs(y3 - y4))
#'
#' #------------------------------------------------------------------------
#' # Generating a sample from a Beta regression model
#' #------------------------------------------------------------------------
#' burn <- 100
#' n <- 500
#' N <- n + burn
#' covar <- cbind(sin(2 * pi * (1:N) / 50), 1:N)
#'
#' set.seed(1234)
#' y1 <- btsr.sim(
#'   model = "BREG", linkg = "logit",
#'   n = n, burn = burn, xreg = covar,
#'   coefs = list(alpha = -1.3, beta = c(0.6, 0.002), nu = 30),
#'   complete = TRUE
#' )
#'
#' # The regressors: X1 = sin(2*pi*t/50) and X2 = t
#' plot.ts(
#'   covar,
#'   main = "Regressors" ~ X[1][t] == sin(2 * pi * t / 50) ~ "and" ~ X[2][t] == t
#' )
#'
#' # Conditional time series:
#' plot.ts(y1$etat, main = "Linear predictor" ~ eta[t] == g[11](mu[t]))
#' plot.ts(y1$mut, main = "Conditional mean" ~ mu[t])
#' plot.ts(y1$yt, main = "Time series" ~ Y[t])
#'
#' #########################################################################
#' #
#' #   Examples of usage of btsr.extract
#' #
#' #########################################################################
#' #------------------------------------------------------------------------
#' # Generating a sample from a BARMAX(1,1) model (BARMA with covariates)
#' #------------------------------------------------------------------------
#' burn <- 100
#' n <- 500
#' N <- n + burn
#' covar <- cbind(sin(2 * pi * (1:N) / 50), 1:N)
#'
#' set.seed(1234)
#' m1 <- btsr.sim(
#'   model = "BARMA", linkg = "logit",
#'   n = n, burn = burn, xreg = covar,
#'   coefs = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ),
#'   error.scale = 1, complete = TRUE
#' )
#'
#' # Extracting components assuming that all coefficients are non-fixed
#' e1 <- btsr.extract(
#'   model = "BARMA", yt = m1$yt,
#'   xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   coefs = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' # Extracting components assuming that all coefficients are fixed
#' #  - no need to provide fixed.lags in this case.
#' e2 <- btsr.extract(
#'   model = "BARMA", yt = m1$yt,
#'   xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   fixed.values = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' # Extracting components assuming that some coefficients are fixed
#' #  - e3 and e4 give the same result
#' #  - e3 uses legacy format for all arguments
#' #  - e4 uses the new format for all arguments (not optimal here)
#' e3 <- btsr.extract(
#'   model = "BARMA", yt = m1$yt,
#'   xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   coefs = list(
#'     phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002)
#'   ),
#'   fixed.values = list(alpha = 0, nu = 20),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' e4 <- btsr.extract(
#'   model = "BARMA", yt = m1$yt,
#'   xreg = list(part1 = covar[(burn + 1):N, ]),
#'   linkg = list(g11 = "logit"),
#'   coefs = list(part1 = list(
#'     phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002)
#'   )),
#'   fixed.values = list(
#'     part1 = list(alpha = 0),
#'     part2 = list(alpha = 20)
#'   ),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' #----------------------------------------------------
#' # comparing the simulated and the extracted values
#' #----------------------------------------------------
#' cbind(head(m1$mut), head(e1$mut), head(e2$mut), head(e3$mut), head(e4$mut))
#'
#' #----------------------------------------------------
#' # comparing the log-likelihood values obtained (must be the all equal)
#' #----------------------------------------------------
#' c(e1$sll, e2$sll, e3$sll, e4$sll)
#'
#' #----------------------------------------------------
#' # comparing the score vectors:
#' #----------------------------------------------------
#' # - e1 must have 6 values: dl/dro values and dl/dlambda values
#' # - e2 must be empty
#' # - e3 and e4 must have only the values corresponding
#' #    to the non-fixed coefficient. The value sof
#' #    dl/dlambda are the same as in e1, but dl/drho are
#' #    differente since the mixed derivatives are not present.
#' #----------------------------------------------------
#' round(e1$score, 4)
#' e2$score
#' e3$score
#' e4$score
#'
#' #----------------------------------------------------
#' # comparing the information matrices.
#' #----------------------------------------------------
#' # - e1 must be a 6x6 matrix
#' # - e2 must be empty
#' # - e3 and e4 must have only the value corresponding
#' #    to the non-fixed coefficient
#' #----------------------------------------------------
#' round(e1$info.Matrix, 4)
#' e2$info.Matrix
#' e3$info.Matrix
#' e4$info.Matrix
#'
#' #------------------------------------------------------------------------
#' # Generating a sample from a BARFIMAVX(1,d1,1)x(1,0,1) with d1 = 0.35
#' # (BARFIMAV with covariates)
#' # Here using the nre format for coefficients and xreg is required.
#' #------------------------------------------------------------------------
#' burn <- 100
#' n <- 500
#' N <- n + burn
#' covar1 <- cbind(sin(2 * pi * (1:N) / 50), 1:N)
#' covar2 <- sin(2 * pi * (1:N) / 25)
#'
#' set.seed(1234)
#' m1 <- btsr.sim(
#'   model = "BARFIMAV",
#'   linkg = list(g11 = "logit", g2 = "default", g21 = "logit"),
#'   n = n, burn = burn, xreg = list(part1 = covar1, part2 = covar2),
#'   coefs = list(
#'     part1 = list(
#'       alpha = 0, d = 0.35, phi = -0.65,
#'       theta = -0.25, beta = c(0.6, -0.002)
#'     ),
#'     part2 = list(
#'       alpha = -3, phi = 0.25,
#'       theta = -0.2, beta = -0.15
#'     )
#'   ),
#'   error.scale = 1, complete = TRUE, vt.start = 0.02
#' )
#'
#'
#' # Extracting components assuming that all coefficients are non-fixed
#' e1 <- btsr.extract(
#'   model = "BARFIMAV", yt = m1$yt,
#'   xreg = list(part1 = covar1[(burn + 1):N, ], part2 = covar2[(burn + 1):N]),
#'   linkg = list(g11 = "logit", g2 = "default", g21 = "logit"),
#'   coefs = list(
#'     part1 = list(
#'       alpha = 0, d = 0.35, phi = -0.65,
#'       theta = -0.25, beta = c(0.6, -0.002)
#'     ),
#'     part2 = list(
#'       alpha = -3, phi = 0.25,
#'       theta = -0.2, beta = -0.15
#'     )
#'   ),
#'   vt.start = 0.02,
#'   llk = TRUE, sco = TRUE, info = TRUE, extra = TRUE, debug = TRUE
#' )
#'
#' # score vector
#' e1$score
#'
#' # information matrix
#' e1$info.Matrix
#'
#' #########################################################################
#' #
#' #   Examples of usage of btsr.fit
#' #
#' #########################################################################
#' #------------------------------------------------------------------------
#' #  Generating a sample from a BARMAVX(1,0)x(0,1) (BARMAV with covariates)
#' #------------------------------------------------------------------------
#' burn <- 100
#' n <- 500
#' N <- n + burn
#' covar1 <- cbind(sin(2 * pi * (1:N) / 50), 1:N)
#' covar2 <- sin(2 * pi * (1:N) / 25)
#'
#' set.seed(1234)
#' m1 <- btsr.sim(
#'   model = "BARMAV",
#'   linkg = list(g11 = "logit", g2 = "default", g21 = "logit"),
#'   n = n, burn = burn, xreg = list(part1 = covar1),
#'   coefs = list(
#'     part1 = list(
#'       alpha = 0, phi = -0.3,
#'       beta = c(0.6, -0.002)
#'     ),
#'     part2 = list(alpha = -2.5, theta = -0.4)
#'   ),
#'   error.scale = 1, complete = TRUE
#' )
#'
#' # fitting the model
#' f1 <- btsr.fit(
#'   model = "BARMAV", yt = m1$yt, report = FALSE, info = TRUE,
#'   xreg = list(part1 = covar1[(burn + 1):N, ]),
#'   linkg = list(g11 = "logit", g2 = "default", g21 = "logit"),
#'   p = c(1, 0), q = c(0, 1)
#' )
#'
#' # fitting the model using the name string for the parent model
#' #   - same result
#' f2 <- btsr.fit(
#'   model = "BARFIMAV", yt = m1$yt, report = FALSE, info = TRUE,
#'   xreg = list(part1 = covar1[(burn + 1):N, ]),
#'   linkg = list(g11 = "logit", g2 = "default", g21 = "logit"),
#'   p = c(1, 0), q = c(0, 1), d = FALSE
#' )
#'
#' summary(f1, full.report = TRUE) # default
#' summary(f2, full.report = TRUE)
#'
#' summary(f1, full.report = FALSE) # simplified output
#' summary(f2, full.report = FALSE)
#' }
#'
NULL
#> NULL


#' @rdname BTSR.functions
#' @order 2
#' @export
btsr.sim <- function(model, ...) {
  # ~~~~~~~~~~~~~~~~~~~~~~~   Models  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # *ARFIMAV: General model with time varying parameters.
  # *ARFIMA: Classical long memory models with nu fixed.
  # *ARMAV: Short memory models (d = 0) with time varying parameters.
  # *ARMA: Classical short memory models (d = 0) with nu fixed.
  # *REGV: Classical regression model with time varying parameters.
  # *REG: Classical regression model with nu fixed.
  #  BETA, GAMMA, KUMA or UW: iid samples.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # - for Compatibility with Previous Versions the format of the
  #   coefficients and link are checked and converted (if necessary)
  # - coefficients and links are set to default values for particular
  #   models
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  model <- toupper(model[1])
  MODEL <- .get.base.model(model)
  if (!(MODEL %in% c("BARC", .current.models(type = "base")))) {
    .stop.with.message("The requested model is not available")
  }

  if (MODEL == "BARC") {
    return(.barc.sim(...))
  }

  return(.btsr.sim(model = model, ...))
}


#' @rdname BTSR.functions
#' @order 3
#' @export
btsr.extract <- function(model, ...) {
  # ~~~~~~~~~~~~~~~~~~~~~~~   Models  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # *ARFIMAV: General model with time varying parameters.
  # *ARFIMA: Classical long memory models with nu fixed.
  # *ARMAV: Short memory models (d = 0) with time varying parameters.
  # *ARMA: Classical short memory models (d = 0) with nu fixed.
  # *REGV: Classical regression model with time varying parameters.
  # *REG: Classical regression model with nu fixed.
  #  BETA, GAMMA, KUMA or UW: iid samples.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # - for Compatibility with Previous Versions the format of the
  #   coefficients and link are checked and converted (if necessary)
  # - coefficients and links are set to default values for particular
  #   models
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  model <- toupper(model[1])
  MODEL <- .get.base.model(model)
  if (!(MODEL %in% c("BARC", .current.models(type = "base")))) {
    .stop.with.message("The requested model is not available")
  }

  if (MODEL == "BARC") {
    return(.barc.extract(...))
  }

  return(.btsr.extract(model = model, ...))
}


#' @rdname BTSR.functions
#' @order 4
#' @export
btsr.fit <- function(model, ...) {
  # ~~~~~~~~~~~~~~~~~~~~~~~   Models  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # *ARFIMAV: General model with time varying parameters.
  # *ARFIMA: Classical long memory models with nu fixed.
  # *ARMAV: Short memory models (d = 0) with time varying parameters.
  # *ARMA: Classical short memory models (d = 0) with nu fixed.
  # *REGV: Classical regression model with time varying parameters.
  # *REG: Classical regression model with nu fixed.
  #  BETA, GAMMA, KUMA or UW: iid samples.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # - for Compatibility with Previous Versions the format of the
  #   coefficients and link are checked and converted (if necessary)
  # - coefficients and links are set to default values for particular
  #   models
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  model <- toupper(model[1])
  MODEL <- .get.base.model(model)
  if (!(MODEL %in% c("BARC", .current.models(type = "base")))) {
    .stop.with.message("The requested model is not available")
  }

  if (MODEL == "BARC") {
    return(.barc.fit(...))
  }

  return(.btsr.fit(model = model, ...))
}


#' @title Predict method for BTSR
#'
#' @description Predicted values based on btsr object.
#'
#' @param object Object of class inheriting from `"btsr"`
#'
#' @param newdata A matrix with new values for the regressors. If omitted and
#'   `"xreg"` is present in the model, the fitted values are returned. If the
#'   model does not include regressors, the functions will use the value of
#'   `nnew`.
#'
#' @param nnew number of out-of-sample forecasts required. If `newdata` is
#'   provided, `nnew` is ignored.
#'
#' @param debug logical, if `TRUE` the output from Fortran is return (for
#'   debugging purposes). Default is `debug = FALSE`.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' `predict.btsr` produces predicted values, obtained by evaluating the
#' regression function in the frame `newdata`.
#'
#' If `newdata` is omitted the predictions are based on the data used for the
#' fit.
#'
#' For now, prediction intervals are not provided.
#'
#' @return
#' If `nnew = 0`, returns a list with in-sample predictions (`fitted.values`,
#' `etat` and `error`), otherwise, returns a list with the following arguments
#'
#' \itemize{
#'  \item `fitted.values`: in-sample forecast.\cr
#'   If \eqn{\nu_t} is fixed: a vector with the in-sample value of \eqn{\mu_t}.\cr
#'   If \eqn{\nu_t} is time varying: a matrix with the in-sample values of
#'   \eqn{\mu_t}, \eqn{\nu_t} and \eqn{\vartheta_t}.
#'
#'  \item `etat`: the linear predictor(s)\cr For models with \eqn{\nu} fixed,
#'   returns \eqn{\eta_{1t} = g_{11}(\mu_t)}\cr For models with time varying
#'   \eqn{\nu}, returns a matrix whose columns are \eqn{\eta_{1t} = g_{11}(\mu_t)}
#'   and \eqn{\eta_{2t} = g_{21}(\vartheta_t)}.
#'
#'  \item `error`: the error term \eqn{e_{1t}} (depends on the argument
#'   \code{error.scale})
#'
#'  \item `residual`: The (in-sample) residuals, that is, the observed values
#'   \eqn{Y_t} minus the fitted values \eqn{\mu_t}. The same as the `error` term
#'   if `error.scale = 0`.
#'
#'  \item `forecast`: the out-of-sample forecast.\cr
#'  If \eqn{\nu_t} is fixed: a vector with the predicted values for \eqn{\mu_t}
#'  and \eqn{\eta_{1t}}\cr
#'  If \eqn{\nu_t} is time varying: a matrix the predicted values for \eqn{\mu_t}
#'  and \eqn{\eta_{1t}}, \eqn{\nu_t}, \eqn{\vartheta_t} and \eqn{\eta_{2t}}.\cr
#'  For BARC models also returs a column with predicted values for the iterated
#'  map.
#'
#'  \item `TS`: only for `"BARC"` models. The iterated map.
#'
#'  \item `xnew`: out-of-sample values for `xreg` (if presented). These are the
#'   values passed through `newdata`.
#' }
#'
#' @examples
#' #------------------------------------------------------------
#' # Generating a Beta model were mut does not vary with time
#' # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
#' #------------------------------------------------------------
#'
#' y <- btsr.sim(
#'   model = "BARFIMA", linkg = "linear",
#'   n = 100, coefs = list(alpha = 0.2, nu = 20)
#' )
#'
#' # fitting the model
#' f <- btsr.fit(
#'   model = "BARFIMA", yt = y, report = TRUE,
#'   start = list(alpha = 0.5, nu = 10),
#'   linkg = "linear", d = FALSE
#' )
#'
#' pred <- predict(f, nnew = 5)
#' pred$forecast
#'
#' @export
#'
predict.btsr <- function(object, newdata, nnew = 0, debug = FALSE, ...) {
  # check if it is a valid model
  if (!.get.base.model(object$model) %in% c("BARC", .current.models("base"))) {
    .stop.with.message(paste0(
      " Model ", object$model, " not implemented yet "
    ))
  }

  # check if newdata is provided
  if (missing(newdata)) newdata <- NULL

  # If wew data was not provided
  nms.out <- .predict.names
  nms.out <- nms.out[nms.out %in% names(object)]

  # if new values are nor required, return the default components
  if (nnew <= 0) {
    return(object[nms.out])
  }

  # if nreg > 0 return an error message
  if (is.null(newdata) && sum(object$configs$order[, "nreg"]) > 0) {
    .stop.with.message(" Please provide the new values for the regressors")
  }

  # get predictions
  if (object$model == "BARC") {
    temp <- .barc.predict(
      obj = object, newdata = newdata, nnew = nnew, debug = debug
    )
  } else {
    temp <- .btsr.predict(
      obj = object, newdata = newdata, nnew = nnew, debug = debug
    )
  }

  # update the object from input
  object <- object[nms.out]
  object[names(temp)] <- temp
  return(object)
}
