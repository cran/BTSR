#----------------------- BTSR MODELS DOCUMENTATION --------------------------
#' @title
#' BTSR models with \eqn{\nu} constant over time
#'
#' @name BTSR.parent.models
#' @order 1
#'
#' @description
#' Function to simulate, extract components, and fit BTSR parent models
#' \itemize{
#' \item \eqn{\nu} constant over time:\cr
#'  BARFIMA, GARFIMA, KARFIMA, MARFIMA, ULARFIMA, and UWARFIMA
#'
#' \item \eqn{\nu} varying over time: \cr
#' BARFIMAV, GARFIMAV, KARFIMAV and UWARFIMAV
#' }
#' all of which are special cases of the general BTSR structure. See the Section
#' \sQuote{The BTSR Structure} in [btsr-package] for details. These functions
#' are maintained for backward compatibility.
#'
#' All models share core arguments with
#' - \code{BARFIMA.sim()} for simulation
#' - \code{BARFIMA.extract()} for extraction
#' - \code{BARFIMA.fit()} for fitting.
#'
#' @inheritParams arguments
#'
#' @details
#' All functions implemented in the current version of the package are
#' compatible with the new format for the arguments, introduced in version
#' 1.0.0. and the previous format.
#' \itemize{
#'  \item The \emph{simulatio functions} (e.g. `BARFIMA.sim`) are used to
#'   generate random samples from the corresponding model.
#'
#'  \item The \emph{extraction functions} (e.g. `BARFIMA.extract`) allow the
#'   user to extract all conditional time series, the log-likelihood, and the
#'   vectors and matrices used to calculate the score vector and the information
#'   matrix associated to a given set of parameters. This function can be used by
#'   any user to create an objective function that can be passed to optimization
#'   functions not available in BTSR Package. At this point, there is no other
#'   use for which this function was intended.
#'
#'  \item The \emph{fitting functions} (e.g. `BARFIMA.fit`) fit the corresponding
#'   model to a given univariate time series. For now, available optimization
#'   algorithms are `"L-BFGS-B"` and `"Nelder-Mead"`. Both methods accept bounds
#'   for the parameters. For `"Nelder-Mead"`, bounds are set via parameter
#'   transformation.
#' }
#'
#' @returns
#' These functions return the same ouptuts as [btsr.sim],
#' [btsr.extract] and [btsr.fit].
#'
#' @seealso
#'  [BTSR.functions], [BARC.functions], [link.btsr], [get.defaults]
#'
#' @examples
#' \donttest{
#' #########################################################################
#' #
#' #   Examples of usage of "MODEL.sim" type of functions
#' #
#' #########################################################################
#' #------------------------------------------------------------------------
#' # iid samples
#' #------------------------------------------------------------------------
#' # BARFIMA: yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
#' # CASE 1: using coefs as in the previous version of the package
#' set.seed(1234)
#' yb1 <- BARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.5, nu = 20)
#' )
#' hist(yb1)
#'
#' # CASE 2: using the new layout
#' set.seed(1234)
#' yb2 <- BARFIMA.sim(
#'   n = 1000,
#'   linkg = list(g11 = "linear", g2 = "linear", g21 = "linear"),
#'   coefs = list(part1 = list(alpha = 0.5), part2 = list(alpha = 20))
#' )
#' hist(yb2)
#'
#' # comparing the results
#' range(abs(yb2 - yb1))
#'
#' # samples from other models in the package
#' yg <- GARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.5, nu = 20)
#' )
#' yk <- KARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.5, nu = 20)
#' )
#' ym <- MARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.5)
#' )
#' yul <- ULARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.5)
#' )
#' yuw <- UWARFIMA.sim(
#'   linkg = "linear", n = 1000,
#'   coefs = list(alpha = 0.5, nu = 20)
#' )
#'
#' # comparing the different distributions
#' par(mfrow = c(2, 3))
#' hist(yb1, xlim = c(0, 1))
#' hist(yk, xlim = c(0, 1))
#' hist(yg, xlim = c(0, 1))
#' hist(ym, xlim = c(0, 1))
#' hist(yul, xlim = c(0, 1))
#' hist(yuw, xlim = c(0, 1))
#'
#' #------------------------------------------------------------------------
#' #  BARFIMA(1,d,1) with d = 0.25 and no regressors
#' #------------------------------------------------------------------------
#'
#' # CASE 1: using coefs as in the previous version of the package
#' set.seed(1234)
#' y1 <- BARFIMA.sim(
#'   n = 1000,
#'   linkg = "logit",
#'   coefs = list(alpha = 0.2, phi = 0.2, theta = 0.4, d = 0.25, nu = 20)
#' )
#'
#' # CASE 2: using the new layout
#' set.seed(1234)
#' y2 <- BARFIMA.sim(
#'   n = 1000,
#'   linkg = list(g11 = "logit", g2 = "linear", g21 = "linear"),
#'   coefs = list(
#'     part1 = list(alpha = 0.2, phi = 0.2, theta = 0.4, d = 0.25),
#'     part2 = list(alpha = 20)
#'   )
#' )
#'
#' # comparing the results
#' range(abs(y1 - y2))
#'
#' #########################################################################
#' #
#' #   Examples of usage of "MODEL.extract" type of functions
#' #
#' #########################################################################
#'
#' #------------------------------------------------------------------------
#' #  code to simulate and extract components of a BARMA(1,1) model
#' #------------------------------------------------------------------------
#' burn <- 100
#' n <- 500
#' N <- n + burn
#' covar <- cbind(sin(2 * pi * (1:N) / 50), 1:N)
#'
#' set.seed(1234)
#' m1 <- BARFIMA.sim(
#'   linkg = "logit", n = n, burn = burn, xreg = covar,
#'   coefs = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ), complete = TRUE
#' )
#'
#' # Extracting assuming that all coefficients are non-fixed
#' e1 <- BARFIMA.extract(
#'   yt = m1$yt, xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   coefs = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' # Extracting assuming that all coefficients are fixed
#' e2 <- BARFIMA.extract(
#'   yt = m1$yt, xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   fixed.values = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' # Extracting using a mix of fixed and non-fixed parameters
#' e3 <- BARFIMA.extract(
#'   yt = m1$yt, xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   coefs = list(
#'     phi = -0.65, theta = -0.25,
#'     beta = c(0.6)
#'   ),
#'   fixed.values = list(alpha = 0, nu = 20, beta = -0.002),
#'   fixed.lags = list(beta = 2),
#'   llk = TRUE, sco = TRUE, info = TRUE
#' )
#'
#' # comparing the simulated and the extracted values of mut
#' cbind(head(m1$mut), head(e1$mut), head(e2$mut), head(e3$mut))
#'
#' # comparing the log-likelihood values obtained (must be the all equal)
#' c(e1$sll, e2$sll, e3$sll)
#'
#' # comparing the score vectors:
#' # - e1 must have 6 values: dl/dro values and dl/dlambda values
#' # - e2 must be empty (all parameters are fixed)
#' # - e3 must have only the values corresponding to the non-fixed coefficients.
#' round(e1$score, 4)
#' e2$score
#' round(e3$score, 4)
#'
#' # comparing the information matrices.
#' # - e1 must be a 6x6 matrix
#' # - e2 must be empty
#' # - e3 must have only the value corresponding to the non-fixed coefficient
#' round(e1$info.Matrix, 4)
#' e2$info.Matrix
#' round(e3$info.Matrix, 4)
#'
#' #########################################################################
#' #
#' #   Examples of usage of "MODEL.fit" type of functions
#' #
#' #########################################################################
#'
#' #------------------------------------------------------------------------
#' #  code to simulate and fit a BARMA(1,1) model
#' #------------------------------------------------------------------------
#' burn <- 100
#' n <- 500
#' N <- n + burn
#' covar <- cbind(sin(2 * pi * (1:N) / 50), 1:N)
#'
#' set.seed(1234)
#' m1 <- BARFIMA.sim(
#'   linkg = "logit", n = n, burn = burn, xreg = covar,
#'   coefs = list(
#'     alpha = 0, phi = -0.65, theta = -0.25,
#'     beta = c(0.6, -0.002), nu = 20
#'   ),
#'   complete = TRUE
#' )
#'
#' plot.ts(m1$yt)
#'
#' # Fit a model assuming that all coefficients are non-fixed
#' # for a more simplified summary, set full.report = FALSE
#' f1 <- BARFIMA.fit(
#'   yt = m1$yt, xreg = covar[(burn + 1):N, ],
#'   linkg = "logit", p = 1, q = 1, report = TRUE
#' )
#'
#' # the fitted coefficients (full model, including d)
#' coefficients(f1)
#'
#' # if you try to use `robust` or `outer` without setting `extra = TRUE`, the
#' # code issues a message and uses the information matrix
#' summary(f1, robust = TRUE)
#' summary(f1, outer = TRUE)
#'
#' # Fit a model assuming alpha and d are fixed
#' f2 <- BARFIMA.fit(
#'   yt = m1$yt, xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   p = 1, q = 1, fixed.values = list(alpha = 0, d = 0)
#' )
#' # Alternatively, set `d = FALSE`
#' f2 <- BARFIMA.fit(
#'   yt = m1$yt, xreg = covar[(burn + 1):N, ], linkg = "logit",
#'   p = 1, q = 1, fixed.values = list(alpha = 0), d = FALSE
#' )
#'
#' # comparing the results
#' true <- c(
#'   alpha = 0, beta = c(0.6, -0.002),
#'   phi = -0.65, theta = -0.25,
#'   d = 0, nu = 20
#' )
#' cf1 <- coefficients(f1)
#' cf2 <- c(NA, coefficients(f2)[1:4], NA, coefficients(f2)[5])
#' round(cbind(true, cf1, cf2), 3)
#' }
NULL
#>NULL

#' @rdname BTSR.parent.models
#' @order 2
#' @export
BARFIMA.sim <- function(
    n = 1, burn = 0, y.start = NULL, xreg = NULL, xreg.start = NULL,
    xregar = TRUE, coefs = NULL, error.scale = 1, linkg = "logit",
    configs.linkg = NULL, inf = 1000, complete = FALSE, debug = FALSE) {
  args <- list(
    model = "BARFIMA", n = n, burn = burn, y.start = y.start,
    xreg = xreg, xreg.start = xreg.start, xregar = xregar,
    coefs = coefs, error.scale = error.scale, linkg = linkg,
    configs.linkg = configs.linkg, inf = inf, complete = complete,
    debug = debug
  )
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
BARFIMAV.sim <- function(
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "BARFIMAV", vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
GARFIMA.sim <- function(linkg = "log", ...) {
  args <- c(list(model = "GARFIMA", linkg = linkg), list(...))
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
GARFIMAV.sim <- function(
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "log", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "GARFIMAV", vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
KARFIMA.sim <- function(rho = 0.5, y.lower = 0, y.upper = 1, ...) {
  args <- c(
    list(model = "KARFIMA", rho = rho, y.lower = y.lower, y.upper = y.upper),
    list(...)
  )
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
KARFIMAV.sim <- function(
    rho = 0.5, y.lower = 0, y.upper = 1,
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "logit"), ...) {
  args <- c(
    list(
      model = "KARFIMAV",
      rho = rho, y.lower = y.lower, y.upper = y.upper,
      vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
MARFIMA.sim <- function(linkg = "cloglog", ...) {
  args <- c(list(model = "MARFIMA", linkg = linkg), list(...))
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
ULARFIMA.sim <- function(...) {
  args <- c(list(model = "ULARFIMA"), list(...))
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
UWARFIMA.sim <- function(rho = 0.5, ...) {
  args <- c(list(model = "UWARFIMA", rho = rho), list(...))
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 2
#' @export
UWARFIMAV.sim <- function(
    rho = 0.5, vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "UWARFIMAV", rho = rho,
      vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.sim, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
BARFIMA.extract <- function(
    yt, xreg = NULL, nnew = 0, xnew = NULL, y.start = NULL, xreg.start = NULL,
    p = NULL, q = NULL, coefs = NULL, lags = NULL, fixed.values = NULL,
    fixed.lags = NULL, xregar = TRUE, error.scale = 1, inf = 1000,
    linkg = "logit", configs.linkg = NULL, m = 0, llk = TRUE,
    sco = FALSE, info = FALSE, extra = FALSE, debug = FALSE) {
  args <- list(
    model = "BARFIMA",
    yt = yt, xreg = xreg, y.start = y.start, xreg.start = xreg.start,
    nnew = nnew, xnew = xnew, p = p, q = q, inf = inf, xregar = xregar,
    coefs = coefs, lags = lags, fixed.values = fixed.values,
    fixed.lags = fixed.lags, error.scale = error.scale, linkg = linkg,
    configs.linkg = configs.linkg, m = m, llk = llk, sco = sco, info = info,
    extra = extra, debug = debug
  )
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
BARFIMAV.extract <- function(
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "BARFIMAV", vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
GARFIMA.extract <- function(linkg = "log", ...) {
  args <- c(list(model = "GARFIMA", linkg = linkg), list(...))
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
GARFIMAV.extract <- function(
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "log", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "GARFIMAV", vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
KARFIMA.extract <- function(rho = 0.5, y.lower = 0, y.upper = 1, ...) {
  args <- c(
    list(model = "KARFIMA", rho = rho, y.lower = y.lower, y.upper = y.upper),
    list(...)
  )
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
KARFIMAV.extract <- function(
    rho = 0.5, y.lower = 0, y.upper = 1,
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "logit"), ...) {
  args <- c(
    list(
      model = "KARFIMAV",
      rho = rho, y.lower = y.lower, y.upper = y.upper,
      vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
MARFIMA.extract <- function(linkg = "cloglog", ...) {
  args <- c(list(model = "MARFIMA", linkg = linkg), list(...))
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
ULARFIMA.extract <- function(...) {
  args <- c(list(model = "ULARFIMA"), list(...))
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
UWARFIMA.extract <- function(rho = 0.5, ...) {
  args <- c(list(model = "UWARFIMA", rho = rho), list(...))
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 3
#' @export
UWARFIMAV.extract <- function(
    rho = 0.5, vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "UWARFIMAV", rho = rho,
      vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.extract, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
BARFIMA.fit <- function(
    yt, xreg = NULL, nnew = 0, xnew = NULL, y.start = NULL, xreg.start = NULL,
    p = NULL, d = FALSE, q = NULL, xregar = TRUE, inf = 1000, start = NULL,
    ignore.start = FALSE, lags = NULL, fixed.values = NULL, fixed.lags = NULL,
    lower = NULL, upper = NULL, error.scale = 1, linkg = "logit",
    configs.linkg = NULL, m = 0, llk = TRUE, sco = FALSE, info = FALSE,
    extra = FALSE, control = NULL, report = TRUE, debug = FALSE, ...) {
  args <- c(
    list(
      model = "BARFIMA",
      yt = yt, xreg = xreg, nnew = nnew, xnew = xnew, y.start = y.start,
      xreg.start = xreg.start, p = p, d = d, q = q, inf = inf, xregar = xregar,
      m = m, start = start, ignore.start = ignore.start, lags = lags,
      fixed.values = fixed.values, fixed.lags = fixed.lags, lower = lower,
      upper = upper, linkg = linkg, configs.linkg = configs.linkg, sco = sco,
      error.scale = error.scale, info = info, extra = extra, control = control,
      report = report, debug = debug
    ),
    list(...)
  )
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
BARFIMAV.fit <- function(
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "BARFIMAV", vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
GARFIMA.fit <- function(linkg = "log", ...) {
  args <- c(list(model = "GARFIMA", linkg = linkg), list(...))
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
GARFIMAV.fit <- function(
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "log", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "GARFIMAV", vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
KARFIMA.fit <- function(rho = 0.5, y.lower = 0, y.upper = 1, ...) {
  args <- c(
    list(model = "KARFIMA", rho = rho, y.lower = y.lower, y.upper = y.upper),
    list(...)
  )
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
KARFIMAV.fit <- function(
    rho = 0.5, y.lower = 0, y.upper = 1,
    vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "logit"), ...) {
  args <- c(
    list(
      model = "KARFIMAV",
      rho = rho, y.lower = y.lower, y.upper = y.upper,
      vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
MARFIMA.fit <- function(linkg = "cloglog", ...) {
  args <- c(list(model = "MARFIMA", linkg = linkg), list(...))
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
ULARFIMA.fit <- function(...) {
  args <- c(list(model = "ULARFIMA"), list(...))
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
UWARFIMA.fit <- function(rho = 0.5, ...) {
  args <- c(list(model = "UWARFIMA", rho = rho), list(...))
  do.call(btsr.fit, args)
}

#' @rdname BTSR.parent.models
#' @order 4
#' @export
UWARFIMAV.fit <- function(
    rho = 0.5, vt.start = NULL, e2.start = NULL,
    linkg = list(g11 = "logit", g2 = "default", g21 = "log"), ...) {
  args <- c(
    list(
      model = "UWARFIMAV", rho = rho,
      vt.start = vt.start, e2.start = e2.start, linkg = linkg
    ),
    list(...)
  )
  do.call(btsr.fit, args)
}
