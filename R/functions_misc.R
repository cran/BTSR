# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Similar to make.link()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Create a Link for BTSR models
#'
#' @description
#' Given the name of a link, this function returns a link function, an inverse
#' link function, the derivative   \eqn{d\eta/d\mu} and the derivative
#' \eqn{d\mu/d\eta}.
#'
#' @param link  character; one of  `"linear"`,`"polynomial"`, `"logit"`,
#'   `"log"`, `"loglog"`, `"cloglog"`, `"SIP"`. See \sQuote{Details} below.
#'   Default is `link = "linear"`.
#'
#' @return An object of class `"link-btsr"`, a list with components
#'  \item{linkfun}{Link function `function(mu)` \eqn{g(\mu)}}
#'  \item{linkinv}{Inverse link function `function(eta)` \eqn{g^{-1}(eta)}}
#'  \item{linkdif}{Derivative `function(mu)` \eqn{d\eta/d\mu}}
#'  \item{mu.eta}{Derivative `function(eta)` \eqn{d\mu/d\eta}}
#'  \item{name}{a name to be used for the link}
#'
#' @details The available links are:
#' \describe{
#'  \item{`linear`:}{\eqn{g(x) = ax}, for \eqn{a} real.}
#'
#'  \item{`polynomial`:}{\eqn{g(x) = ax^b}, for \eqn{a} and \eqn{b} real.}
#'
#'  \item{`logit`:}{\eqn{g(x) = \log((x-l)/(u-x))}}
#'
#'  \item{`log`:}{\eqn{g(x) = \log(x - l)}}
#'
#'  \item{`loglog`:}{\eqn{g(x) = \log(-\log((x-l)/(u - l)))}}
#'
#'  \item{`cloglog`:}{\eqn{g(x) = \log(-\log(1-(x-l)/(u - l)))}}
#'
#'  \item{`SIP` (Shifted Inverse Power):}{\eqn{g(x) = 1/(a+x)^b}, with
#'   \eqn{a \in \{0,1\}} and \eqn{b} real.}
#' }
#' Here \eqn{l} and \eqn{u} denote the lower and upper bounds of \eqn{x}. The
#' linear link is a particular case of the polynomial link. It was kept for
#' compatibility with older versions of the package.
#'
#' The parameters \eqn{a}, \eqn{b}, \eqn{l}, and \eqn{u} are configured using
#' the `configs.linkg` list, which can include the following elements
#'  \itemize{
#'   \item `ctt`: A constant multiplier for the link function (default: 1).
#'   \item `power`: The power parameter for polynomial and SIP links (default: 1).
#'   \item `lower`: The lower bound for `mu` (default: 0).
#'   \item `upper`: The upper bound for `mu` (default: 1).
#'  }
#' This list must be  passed to the functions created by `link.btsr`, when
#' invoking them.
#'
#' @examples
#' \donttest{
#' #---------------------------------------------
#' #  0 < y < 1 and linear link
#' #---------------------------------------------
#' mylink <- link.btsr("linear")
#' y <- 0.8
#' a <- 3.4
#' gy <- a * y
#'
#' # comparing the expected result with the output of the function:
#' mylink$linkfun(mu = y, configs.linkg = list(ctt = a))
#' gy
#'
#' mylink$linkinv(eta = gy, configs.linkg = list(ctt = a))
#' y
#'
#' mylink$diflink(mu = y, configs.linkg = list(ctt = a))
#' a
#'
#' mylink$mu.eta(eta = gy, configs.linkg = list(ctt = a))
#' 1 / a
#'
#' # same function, different ctt:
#' mylink$linkfun(mu = y, configs.linkg = list(ctt = a + 1))
#'
#' #---------------------------------------------
#' # For linear link bounds have no effect
#' #---------------------------------------------
#' mylink <- link.btsr("linear")
#' y <- 0.8
#' a <- 3.4
#' gy <- a * y
#'
#' mylink$linkfun(mu = y, configs.linkg = list(ctt = a, lower = 1, upper = 2))
#' mylink$linkfun(mu = y, configs.linkg = list(ctt = a)) # same result
#' gy
#'
#' mylink$linkinv(eta = gy, configs.linkg = list(ctt = a, lower = 1, upper = 2))
#' y
#'
#' mylink$diflink(mu = y, configs.linkg = list(ctt = a, lower = 1, upper = 2))
#' a
#'
#' mylink$mu.eta(eta = gy, configs.linkg = list(ctt = a, lower = 1, upper = 2))
#' 1 / a
#'
#' #---------------------------------------------
#' # 0 < y < 1 and logit link
#' #---------------------------------------------
#' mylink <- link.btsr("logit")
#' y <- 0.8
#' gy <- log(y / (1 - y))
#' ginv <- 1 / (1 + exp(-gy))
#'
#' mylink$linkfun(mu = y)
#' gy
#'
#' mylink$linkinv(eta = gy)
#' ginv
#'
#' mylink$diflink(mu = y)
#' 1 / (y - y**2)
#'
#' mylink$mu.eta(eta = gy)
#' ginv - ginv**2
#'
#' #---------------------------------------------
#' # 1 < y < 2 and logit link
#' #---------------------------------------------
#' mylink <- link.btsr("logit")
#' a <- 1
#' b <- 2
#' y <- 1.8
#' gy <- log((y - a) / (b - y))
#' ginv <- b / (1 + exp(-gy)) + a / (1 + exp(gy))
#'
#' mylink$linkfun(mu = y, configs.linkg = list(lower = 1, upper = 2))
#' gy
#'
#' mylink$linkinv(eta = gy, configs.linkg = list(lower = 1, upper = 2))
#' ginv
#'
#' mylink$diflink(mu = y, configs.linkg = list(lower = 1, upper = 2))
#' (b - a) / ((a + b) * y - y**2 - a * b)
#'
#' mylink$mu.eta(eta = gy, configs.linkg = list(lower = 1, upper = 2))
#' ((a + b) * ginv - ginv**2 - a * b) / (b - a)
#' }
#' @export
link.btsr <- function(link = "linear") {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # linkfun:  Link function function(mu)
  # linkinv:  Inverse link function function(eta)
  # mu.eta:   Derivative function(eta) dmu/deta
  # diflink:  Derivative function(mu) deta/dmu
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # convert character to number
  linktemp <- .convert.link(link)

  # helper function:
  # set common parameters and calls FORTRAN
  .call.linkR <- function(linktemp, x, gx, dlink, id, configs.linkg) {
    # configurations for the link
    # update if necessary
    cf <- list(ctt = 1, power = 1, lower = 0, upper = 1)
    if (!is.null(configs.linkg)) cf[names(configs.linkg)] <- configs.linkg
    # length of the input argument (y, mu or eta)
    n <- length(x)
    .Fortran("linkr",
      link = linktemp, par = c(cf$lower, cf$upper, cf$ctt, cf$power),
      n = n, ind = as.integer(id), x = x, gx = gx, dlink = dlink, NAOK = TRUE, PACKAGE = "BTSR"
    )
  }

  # defines g(mu)
  linkfun <- function(mu, ...) {
    args <- list(...)
    .call.linkR(
      linktemp = linktemp, x = mu, gx = numeric(length(mu)), dlink = 0,
      id = c(1, 0, 0), configs.linkg = args$configs.linkg
    )$gx
  }


  # defines g^{-1}(eta)
  linkinv <- function(eta, ...) {
    args <- list(...)
    .call.linkR(
      linktemp = linktemp, x = numeric(length(eta)), gx = eta, dlink = 0,
      id = c(0, 1, 0), configs.linkg = args$configs.linkg
    )$x
  }

  # defines dg/dmu
  diflink <- function(mu, ...) {
    args <- list(...)
    .call.linkR(
      linktemp = linktemp, x = mu, gx = 0, dlink = numeric(length(mu)),
      id = c(0, 0, 1), configs.linkg = args$configs.linkg
    )$dlink
  }

  # defines dmu/deta = 1/g'(ginv(eta))
  mu.eta <- function(eta, ...) {
    1 / diflink(mu = linkinv(eta = eta, ...), ...)
  }

  structure(list(
    linkfun = linkfun,
    linkinv = linkinv,
    diflink = diflink,
    mu.eta = mu.eta,
    name = link
  ), class = "link-btsr")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Default control list
#'
#' @description
#' Sets default values for constants used by the optimization functions in
#' FORTRAN.
#'
#' @param control a list with configurations to be passed to the optimization
#'   subroutines. Missing arguments will receive default values. See
#'   \sQuote{Details}.
#'
#' @details
#' The `control` argument is a list that can supply any of the following
#' components
#'  \describe{
#'   \item{`method`}{The optimization method. Current available options are
#'   `"L-BFGS-B"` and `"Nelder-Mead"`. Default is `"L-BFGS-B"`.}
#'
#'   \item{`maxit`}{The maximum number of iterations. Defaults is `1000`.}
#'
#'   \item{`iprint`}{The frequency of reports if `control$trace` is positive.
#'   Defaults is -1 (no report).
#'   \itemize{
#'    \item For `"L-BFGS-B"` method:
#'    \itemize{
#'     \item `iprint < 0`  no output is generated;
#'     \item `iprint = 0`  print only one line at the last iteration;
#'     \item `0 < iprint < 99` print also f and |proj g| every iprint iterations;
#'     \item `iprint = 99` print details of every iteration except n-vectors;
#'     \item `iprint = 100`  print also the changes of active set and final x;
#'     \item `iprint > 100`  print details of every iteration including x and g;
#'    }
#'
#'   \item For `"Nelder-Mead"` method:
#'    \itemize{
#'    \item `iprint < 0` No printing
#'    \item `iprint = 0` Printing of parameter values and the function value
#'    after initial evidence of convergence.
#'    \item `iprint > 0` As for `iprint = 0` plus progress reports after every
#'    `iprint` evaluations, plus printing for the initial simplex.
#'     }
#'  }}
#'
#'   \item{`factr`}{controls the convergence of the `"L-BFGS-B"`  method.
#'   Convergence occurs when the reduction in the objective is within this
#'   factor of the machine tolerance. The iteration will stop when
#'   \deqn{\dfrac{(f^k - f^{k+1})}{\max\{|f^k|,|f^{k+1}|,1\}} \le
#'   factr \times epsmch}
#'   where \eqn{epsmch} is the machine precision, which is automatically
#'   generated by the code. Typical values for `factr`: 1.e+12 for low accuracy;
#'   1.e+7 for moderate accuracy; 1.e+1 for extremely high accuracy. Default is
#'   `1e7`, that is a  tolerance of about `1e-8`.}
#'
#'   \item{`pgtol`}{helps control the convergence of the `"L-BFGS-B"` method.
#'   It is a tolerance on the projected gradient in the current search
#'   direction. the iteration will stop when
#'   \deqn{\max\{|\text{proj }g_i |, i = 1, ..., n\} \le pgtol}
#'   where \eqn{\text{proj }g_i} is the ith component of the projected gradient.
#'   Default is `1e-12`.}
#'
#'   \item{`stopcr`}{The criterion applied to the standard deviation of the
#'   values of objective function at the points of the simplex, for
#'   `"Nelder-Mead"` method.}
#' }
#'
#' @return
#' Returns a list with all arguments in \sQuote{Details}.
#'
#' @examples
#' BTSR::fit.control()
#'
#' @export
fit.control <- function(control = list()) {
  con <- list(
    method = "L-BFGS-B",
    maxit = 1000,
    iprint = -1,
    factr = 1e+7,
    pgtol = 1e-12,
    stopcr = 1e-4
  )

  con[names(control)] <- control

  return(con)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Initial values for coefficients
#'
#' @description
#' Generates initial values for coefficients in BTSR models.
#'
#' @details
#' Parameter initialization is done as follows.
#'
#' \enumerate{
#'  \item  Legacy flat lists are converted to nested `part1`/`part2` format.
#'   Link functions and density bounds are validated.
#'
#'  \item \strong{Part 1:} \eqn{\mu_t} related parameters.
#'
#'   A linear model is used to estimate \eqn{\alpha}, \eqn{\boldsymbol
#'   \beta} and \eqn{\boldsymbol \phi} by setting
#'   \deqn{
#'   \boldsymbol{Y} = \begin{pmatrix}
#'   Y_1 \\ \vdots \\ Y_n
#'   \end{pmatrix}
#'   \quad \text{and} \quad
#'   D = \begin{pmatrix}
#'   1 & X_{11} & \cdots & X_{1s} & g_{12}(Y_0) & \cdots & g_{12}(Y_{1-p}) \\
#'   \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
#'   1 & X_{n1} & \cdots & X_{ns} & g_{12}(Y_{n-1}) & \cdots & g_{12}(Y_{n-p})
#'   \end{pmatrix}
#'  }
#'  where \eqn{s} is the number of regressors and \eqn{p} is the AR order, and
#'  solving the linear regression problem \eqn{\boldsymbol{Y} =
#'  D\boldsymbol{\gamma} + \boldsymbol{\epsilon}} via `lm.fit`.
#'
#'  MA coefficients \eqn{\boldsymbol{\theta}} are initialized to zero (though
#'  small non-zero values may help with optimization stability)
#'
#'  The long-memory parameter \eqn{d} starts at 0.01 when estimated
#'
#'  For BARC models:
#'  \itemize{
#'   \item Map parameters use:\cr
#' ```
#'    +-------+-------+-------+-------+-------+-------+
#'    |  map  |   1   |   2   |   3   |   4   |   5   |
#'    +-------+-------+-------+-------+-------+-------+
#'    | theta |   3   |  0.5  |  3.5  |  0.5  |  NA   |
#'    +-------+-------+-------+-------+-------+-------+
#' ```
#'   \item \eqn{u_0} defaults to \eqn{\pi/4} when not fixed
#'  }
#'
#'   \item \strong{Part 2:} \eqn{\nu_t} related parameters.
#'
#'   If presented and not time varying, \eqn{\nu} is initialized as follows:
#'   \itemize{
#'    \item \eqn{\nu = 5}, for the Kumaraswamy and the Unit-Weibull distributions,
#'    \item \eqn{
#'     \nu = \displaystyle\frac{1}{n}\sum_{t=1}^n\dfrac{\hat \mu_t (1 - \hat \mu_t)}{\sigma^2} - 1},
#'     for the Beta distribution,
#'    \item \eqn{\nu = \displaystyle\frac{1}{n}\sum_{t=1}^n\frac{\hat \mu_t^2}{\sigma^2}},
#'     for the Gamma distribution
#'   }
#'   where \eqn{(\hat\mu_1, \cdots, \hat\mu_n)' = D\hat{\boldsymbol{\gamma}}}
#'   are the fitted values from the regression model and \eqn{\sigma^2} is the
#'   estimated variance of the residuals.
#'
#'   If \eqn{\nu} is time varying,
#'   \itemize{
#'    \item set \eqn{\alpha} as \eqn{g_{12}(g_2(\nu))}, with \eqn{\nu}
#'     estimated as in the case where the parameter does not vary on time.
#'
#'    \item set \eqn{\boldsymbol{\beta}},  \eqn{\boldsymbol{\phi}} and
#'     \eqn{\boldsymbol{\theta}} to zero.
#'
#'    \item The long-memory parameter \eqn{d} starts at 0.01 when estimated.
#'   }
#' }
#'
#' @inheritParams arguments
#'
#' @param yt numeric vector with the observed time series. Missing values (NA's)
#'   are not allowed.
#'
#' @param p the AR order. Default is `p = 0`. Can be provided as either a single
#'   integer (legacy format) or a length 2 integer vector (new format)
#'   specifying orders for `part1`/`part2`. If \eqn{\nu} is time-varying and a
#'   single value is provided it is assumed that \eqn{p_1 = p_2 = p}.
#'
#' @param q the MA order. Default is `q = 0`. Can be provided as either a single
#'   integer (legacy format) or a length 2 integer vector (new format)
#'   specifying orders for `part1`/`part2`. If \eqn{\nu} is time-varying and a
#'   single value is provided it is assumed that \eqn{q_1 = q_2 = q}.
#'
#' @param d a length 1 (legacy format) or 2 (new format) logical vector
#'   indicating whether the long memory parameter \eqn{d} is presented in the
#'   model (either as a fixed or non-fixed parameter). In the new format, if
#'   only one value is provided the code assumes that the same option is valid
#'   for both parts of the model. Default is `d = FALSE`.
#'
#' @param lags optional; specification of which lags to include in the model.
#'   For one dimensional coefficients, the lag is obviously always 1 and can
#'   be suppressed. Can be specified in one of two ways
#'   \itemize{
#'    \item a list with components `beta`, `phi`, and `theta` (legacy format)
#'     specifying which lags to include for each parameter type.
#'
#'    \item a list with elements `part1` and `part2` (new format), each being a
#'    list with components `beta`, `phi`, and `theta` specifying which lags to
#'    include for each parameter type.
#'   }
#'   Default is `lags = NULL`, in which the `lags` are computed from the
#'   `fixed.lags` argument (if provided). When components are missing or empty
#'   in both, `lags` and `fixed.lags`, the default behavior is to include all
#'   lags based on `nreg = ncol(xreg)`, `p`, and `q`. The arguments `lags` and
#'   `fixed.lags` are complementary. Either suffices, or mix them (e.g., `lags`
#'   for some parameters, `fixed.lags` for others).
#'
#' @param fixed.values optional; specification of fixed parameter values.
#'   Can be specified in one of two ways
#'   \itemize{
#'    \item a list with optional components `alpha`, `beta`, `phi`, `theta`, `d`
#'    and `nu` (legacy format) containing fixed values for each parameter type.
#'
#'    \item a list with elements `part1` and `part2` (new format), each being a
#'    list with optional components `alpha`, `beta`, `phi`, `theta` and `d`
#'    containing fixed values for each parameter type.
#'   }
#'   If fixed values are provided and there exists more than one possible lag,
#'   either `lags` or `fixed.lags` must also be provided. The default is
#'   `fixed.lags = NULL`. By default, if a given vector (such as the vector of
#'   AR coefficients) has fixed values and the corresponding entry in this list
#'   is empty, the fixed values are set as zero.
#'
#' @param fixed.lags optional; specification of which lags should be fixed.
#'   For one dimensional coefficients, the lag is obviously always 1 and can
#'   be suppressed. Can be specified in one of two ways
#'   \itemize{
#'    \item a list with components `beta`, `phi`, and `theta` (legacy format)
#'     specifying which lags should be fixed.
#'
#'    \item a list with elements `part1` and `part2` (new format), each being a
#'     list with components `beta`, `phi`, and `theta` specifying which lags
#'     should be fixed.
#'   }
#'   For missing components, fixed values will are set based on `lags`.
#'
#' @param linkg specification of link functions. Can be specified in one of two
#'   ways
#'   \itemize{
#'    \item A character or two-character vector (legacy format). If only one
#'     string is provided, the same link name is used for \eqn{\mu_t}
#'     (\eqn{g_{11}}) and for \eqn{Y_t} (\eqn{g_{12}}).
#'
#'    \item A named list (new format) with elements `g11`, `g12`, `g2`, `g21`,
#'    and `g22` (order does not matter). Particular models (see
#'    \sQuote{Particular Models} in [BTSR.functions]) have default values for
#'    some links. Missing links follow these rules
#'    \itemize{
#'     \item If either `g11` or `g12` is missing (but not both), assumes `g12 = g11`
#'     \item If `phi = NULL` for part 1, `g12` is not required
#'     \item If `phi = NULL` for part 2, `g22` is not required
#'     \item If either `g21` or `g22` is missing (but not both), assumes `g22 = g21`
#'    }
#'    Special case: `g2 = "default"` uses the model's default link. The default
#'    depends on the model.
#'   }
#'   Default is `linkg = "linear"`, which is equivalent (done internally) to set
#'   all links as `"linear"` in the new format. See [link.btsr] for valid
#'   links. For details, see the Section \sQuote{The BTSR structure} in
#'   [btsr-package].
#'
#' @return
#' For models where \eqn{\nu_t} is not time-varying, returns a list (legacy
#' format) with starting values for the parameters of the selected model.
#' Possible outputs are
#'
#' \item{alpha}{ the intercept.}
#' \item{beta}{ the coefficients for the regressors.}
#' \item{phi}{ the AR coefficients.}
#' \item{theta}{ for BARC models, the parameter associate to the map function.
#'  For any other model, the MA coefficients.}
#' \item{d}{ the long memory parameter.}
#' \item{nu}{ distribution related parameter (usually, the precision).}
#' \item{u0}{ for BARC models, the initial value of the iterated map.}
#'
#' For models where \eqn{\nu_t} is time-varying, returns a list whose elements
#' are `part1` and `part2`. Each element is a list with starting values for the
#' parameters corresponding to each part o the selected model. Possible outputs
#' for each part are the same as for the legacy format.
#'
#' @examples
#' mu <- 0.5
#' nu <- 20
#'
#' yt <- rbeta(100, shape1 = mu * nu, shape2 = (1 - mu) * nu)
#' # using the general model BARFIMA
#' coefs.start(model = "BARFIMA", yt = yt, linkg = "linear")
#' # same output as the specific model BETA
#' coefs.start(model = "BETA", yt = yt, linkg = "linear")
#'
#' yt <- rgamma(100, shape = nu, scale = mu / nu)
#' coefs.start(model = "GARFIMA", yt = yt, linkg = "linear")
#'
#' @importFrom  stats lm.fit fitted residuals
#'
#' @export
coefs.start <- function(
    model, yt, y.start = NULL, y.lower = 0, y.upper = 1,
    xreg = NULL, p = 0, q = 0, d = FALSE, map = .default.map.barc,
    lags = NULL, fixed.values = NULL, fixed.lags = NULL,
    linkg = "linear", configs.linkg = NULL) {
  # check BARC defaults
  is.barc <- model == "BARC"
  if (is.barc && is.null(map)) {
    map <- .default.map.barc
    .warning.with.message(paste0(
      "'map' is missing. Using the default value map = ", map
    ))
  }

  n <- length(yt)
  # If y.lower or y.upper is missing, use the default values
  defaults <- .get.model.settings(model = model)
  if (defaults$model != "KARFIMA" || is.null(y.lower)) {
    y.lower <- defaults$lconfig["g12", 1]
  }
  if (defaults$model != "KARFIMA" || is.null(y.upper)) {
    y.upper <- defaults$lconfig["g12", 2]
  }

  # for compatibility with the old structure:
  #   - checking the structure and converting to the new format
  lags <- .check.coefs.defaults(
    model = model, object = lags, type = "lags"
  )
  fixed.values <- .check.coefs.defaults(
    model = model, object = fixed.values, type = "fixed.values"
  )
  fixed.lags <- .check.coefs.defaults(
    model = model, object = fixed.lags, type = "fixed.values"
  )
  temp <- .convert.xreg(
    model = model, n = n, nnew = 0, xreg = xreg, xnew = NULL, skip.forecast = TRUE
  )
  xreg <- temp$xreg
  nreg <- temp$nreg
  if (length(p) == 1) p <- c(p, p)
  if (length(q) == 1) q <- c(q, q)
  if (length(d) == 1) d <- c(d, d)

  # checking the the format of the link.
  # Setting default values if necessary
  if (is.barc) {
    allow.any <- TRUE
  } else {
    allow.any <- apply(
      cbind(nreg = nreg, p = p, q = q), 1, function(x) sum(x) == 0
    )
  }
  temp <- .check.link(
    model = model,
    link = list(
      linkg = linkg,
      configs.linkg = configs.linkg, escale = NULL
    ),
    lk.default = list(
      linkg = defaults$linkg,
      configs.linkg = list(
        ctt = defaults$lconfig[, "ctt"], power = defaults$lconfig[, "power"]
      )
    ),
    allow.any = allow.any
  )
  linkg <- temp$link
  configs.linkg <- temp$configs

  if (is.null(y.start)) y.start <- .get.link.starting.values(linkg$g12)

  #
  #  Starting values for part 1:
  #

  # helper function:
  #  updates the values in configs list to pass to the link function
  .update.cf <- function(key) {
    return(list(
      ctt = 1, power = configs.linkg$power[[key]],
      lower = y.lower, upper = y.upper
    ))
  }

  # creates the link functions g11 and g12
  # and computing the values g11(yt) and g12(yt)
  linktemp1 <- link.btsr(link = linkg$g11)
  g11 <- linktemp1$linkfun
  g1y <- g11(yt, configs.linkg = .update.cf("g11"))
  check <- .compare.link(
    link = linkg, configs.linkg = configs.linkg, key1 = "g11", key2 = "g12"
  )
  # if g11 = g12 skip creating g12
  if (check) {
    g2y <- g1y
    g12 <- g11
  } else {
    linktemp2 <- link.btsr(link = linkg$g12)
    g12 <- linktemp2$linkfun
    g2y <- g12(yt, configs.linkg = .update.cf("g12"))
  }

  # starting values for alpha, phi and beta.
  # Using the coefficients of the regression
  #  g(Y) = a + Xb + Gc
  # where G is the matrix with lagged values of g(Y)
  X <- matrix(1, nrow = n)
  nreg1 <- nreg[1]
  if (nreg1 > 0) {
    lag <- 1:nreg[1]
    if (!is.null(lags$part1$beta)) {
      lag <- lags$part1$beta
    } else {
      if (!is.null(fixed.lags$part1$beta)) {
        fl <- fixed.lags$part1$beta
        lag <- lag[-fl]
      }
    }
    nreg1 <- length(lag)
    if (nreg1 > 0) X <- cbind(X, as.matrix(xreg$part1)[, lag])
  }
  p1 <- p[1]
  if (p1 > 0) {
    gystart <- g12(y.start, configs.linkg = .update.cf("g12"))
    lag <- 1:p[1]
    if (!is.null(lags$part1$phi)) {
      lag <- lags$part1$phi
    } else {
      if (!is.null(fixed.lags$part1$phi)) {
        fl <- fixed.lags$part1$phi
        lag <- lag[-fl]
      }
    }
    p1 <- length(lag)
    if (p1 > 0) {
      P <- matrix(gystart, nrow = n, ncol = p1)
      for (i in 1:p1) P[-c(1:lag[i]), i] <- g2y[1:(n - lag[i])]
      X <- cbind(X, P)
    }
  }
  w <- sum(is.na(X[, ncol(X)]))
  if (w > 0) {
    X <- X[-c(1:w), ]
    g1y <- g1y[-c(1:w)]
  }

  fit <- lm.fit(x = X, y = g1y)
  mqo <- c(fit$coefficients, use.names = FALSE)
  mqo[is.na(mqo)] <- 0
  k <- length(mqo)

  # initializing the parameter values
  start <- list()

  a <- as.integer(is.null(fixed.values$part1$alpha))
  if (a == 1) {
    start$alpha <- mqo[1]
  } else {
    mqo <- mqo[-1]
    k <- k - 1
  }

  if (nreg1 > 0) start$beta <- mqo[(a + 1):(a + nreg1)]

  if (p1 > 0) start$phi <- mqo[(a + nreg1 + 1):k]

  q1 <- max(
    q[1] - max(length(fixed.values$part1$theta), length(fixed.lags$part1$theta)),
    length(lags$part1$theta)
  )
  if (q1 > 0) {
    if (is.barc) {
      start$theta <- .get.barc.starting.theta(map)
    } else {
      start$theta <- rep(0, q1)
    }
  }

  if (d[1] == TRUE) {
    if (is.null(fixed.values$part1$d)) start$d <- 0.01
  }

  if (is.barc && is.null(fixed.values$part1$u0)) {
    start$u0 <- .default.u0.barc
  }

  if (endsWith(model, "V")) start <- list(part1 = start)

  #
  #  Starting values for part 2:
  #

  nu <- NULL
  npar2 <- nreg[2] + p[2] + q[2] + is.null(fixed.values$part2$alpha)

  # if part2 is fixed and nu is provided return
  if (npar2 == 0) {
    return(start)
  }

  # part2 is not fixed or nu is not provided
  if (is.null(fixed.values$part2$alpha)) {
    n1 <- length(g1y)
    mu <- fitted(fit)
    mu <- linktemp1$linkinv(mu, configs.linkg = .update.cf("g11"))
    dlink <- linktemp1$diflink(mu, configs.linkg = .update.cf("g11"))
    er <- residuals(fit)
    sigma2 <- sum(er^2 / ((n1 - k) * (dlink)^2)) # estimated variance
    nu.type <- switch(EXPR = .get.base.model(model[1]),
      BARC = "T1",
      BARFIMA = "T1",
      GARFIMA = "T2",
      KARFIMA = "T3",
      MARFIMA = "T4",
      ULARFIMA = "T4",
      UWARFIMA = "T3"
    )
    nu <- NULL
    nu <- switch(EXPR = nu.type,
      T1 = mean(mu * (1 - mu) / sigma2) - 1,
      T2 = mean(mu^2 / sigma2),
      T3 = 5
    )

    # if nu is fixed, add the initial value to the list and return
    if (!endsWith(model, "V")) {
      start$nu <- nu
      return(start)
    }

    # g2(nut)
    linktemp3 <- link.btsr(link = linkg$g2)
    g2 <- linktemp3$linkfun
    cf <- list(
      ctt = 1, power = configs.linkg$power[["g2"]],
      lower = defaults$lconfig["g2", 1], upper = defaults$lconfig["g2", 2]
    )
    gnu <- g2(nu, configs.linkg = cf)
    check <- .compare.link(
      link = linkg, configs.linkg = configs.linkg, key1 = "g2", key2 = "g21"
    )
    if (check) {
      gvt <- gnu
    } else {
      linktemp4 <- link.btsr(link = linkg$g21)
      g21 <- linktemp4$linkfun
      aa <- g2(defaults$lconfig["g2", 1], configs.linkg = cf)
      bb <- g2(defaults$lconfig["g2", 2], configs.linkg = cf)
      cf$lower <- min(aa, bb)
      cf$upper <- max(aa, bb)
      cf$power <- configs.linkg$power[["g21"]]
      gvt <- g21(gnu, configs.linkg = cf)
    }
    start$part2$alpha <- gvt
  }

  nreg2 <- max(
    nreg[2] - max(length(fixed.values$part2$beta), length(fixed.lags$part2$beta)),
    length(lags$part2$beta)
  )
  if (nreg2 > 0) start$part2$beta <- rep(0, nreg2)

  p2 <- max(
    p[2] - max(length(fixed.values$part2$phi), length(fixed.lags$part2$phi)),
    length(lags$part2$phi)
  )
  if (p2 > 0) start$part2$phi <- rep(0, p2)

  q2 <- max(
    q[2] - max(length(fixed.values$part2$theta), length(fixed.lags$part2$theta)),
    length(lags$part2$theta)
  )
  if (q2 > 0) start$part2$theta <- rep(0, q2)

  if (d[2] == TRUE) {
    if (is.null(fixed.values$part2$d)) start$part2$d <- 0.01
  }

  return(start)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Summary Method of class BTSR
#'
#' @description \code{summary} method for class \code{"btsr"}.
#'
#' @name summary
#'
#' @aliases summary.btsr
#' @aliases print.summary.btsr
#'
#' @param object object of class `"btsr"`.
#'
#' @param robust logical. If \code{TRUE} the robust covariance matrix is
#'   computed
#'
#' @param outer logical. If \code{TRUE} uses the outer product of the gradient
#'   to compute the covariance matrix. If \code{robust = TRUE}, \code{outer} is
#'   used as a second option (in case of error computing the robust version)
#'
#' @param full.report logical. If \code{TRUE} prints a more detailed report.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' The function `summary.btsr` computes and returns a list of summary statistics
#' of the fitted model given in `object`. Returns a list of class
#' `summary.btsr`, which contains the following components
#'
#' \itemize{
#'  \item `model`: the corresponding model.
#'
#'  \item `call`: the matched call.
#'
#'  \item `residuals`: The (in-sample) residuals, that is, the observed values
#'  \eqn{Y_t} minus the fitted values \eqn{\mu_t}. The same as the `error` term
#'  if `error.scale = 0`.
#'
#'  \item `coefficients`: a \eqn{k \times 4}{k x 4} matrix with columns for the
#'  estimated coefficient, its standard error, z-statistic and corresponding
#'  (two-sided) p-value.
#'
#'  \item `sigma.res`: the square root of the estimated variance of the error
#'  term in `residuals`
#'  \deqn{
#'  \hat\sigma^2 = \displaystyle\frac{1}{n-k}\sum_{i=1}^{n-k}{e_i^2},
#'  }
#'  where \eqn{e_i} is the \eqn{i}-th residual.
#'
#'  \item `df`: degrees of freedom, a 2-vector \eqn{(k, n-k)}, the first being
#'  the number of estimated coefficients.
#'
#'  \item `vcov`: a \eqn{k \times k}{k \times k} matrix of (unscaled)
#'  covariances. The inverse ov the information matrix.
#'
#'  \item `loglik`: the sum of the log-likelihood values (\eqn{L})
#'
#'  \item `aic`: the AIC value. \eqn{AIC = -2L + 2k}.
#'
#'  \item `bic`: the BIC value. \eqn{BIC = -2L + k\log(n)}.
#'
#'  \item `hqc`: the HQC value. \eqn{HQC = -2L + k\log(\log(n))}.
#' }
#'
#' @importFrom stats pnorm
#'
#' @export
summary.btsr <- function(
    object, robust = FALSE, outer = FALSE, full.report = TRUE, ...) {
  if (!"btsr" %in% class(object)) {
    .stop.with.message(" The argument 'object' must be a 'btsr' object")
  }

  ans <- list(
    model = object$model,
    call = object$call,
    full.report = full.report,
    algorithm = list(
      method = object$control$method,
      neval = object$counts,
      convergence = ifelse(object$convergence == 0, TRUE, FALSE)
    ),
    residuals = object$residuals
  )

  npar <- length(object$coefficients)
  n <- length(object$residuals)
  rdf <- ans$df.residuals <- n - npar
  ans$sigma.res <- sqrt(sum(ans$residuals^2) / rdf)
  class(ans) <- "summary.btsr"

  if (npar == 0) {
    ans$df <- c(0L, n)
    ans$coefficients <- matrix(NA_real_, 0L, 4L,
      dimnames = list(NULL, c(
        "Estimate", "Std. Error",
        "z value", "Pr(>|t|)"
      ))
    )
  } else {
    ans$df <- c(npar, rdf)

    gradient <- NULL
    if (object$configs$extra == 1) gradient <- .gradient(object)
    vcov.temp <- .cov.matrix(robust, outer, object$info.Matrix, gradient)
    ans$vcov <- vcov.temp$cov
    ans$vcov.type <- vcov.temp$type
    stderror <- sqrt(diag(abs(ans$vcov)))
    zstat <- abs(object$coefficients / stderror)
    ans$coefficients <- cbind(
      Estimate = object$coefficients,
      `Std. Error` = stderror,
      `z value` = zstat,
      `Pr(>|t|)` = 2 * (1 - pnorm(zstat))
    )
  }

  ans$loglik <- object$sll
  ans$aic <- -2 * ans$loglik + 2 * npar
  ans$bic <- -2 * ans$loglik + log(n) * npar
  ans$hqc <- -2 * ans$loglik + log(log(n)) * npar

  nms <- .valid.links("code.base")
  ans$escale <- object$link["g13"]

  # helper function
  # prints the name of the link and the constants
  # for SIP and polynomial links
  f <- function(x, lk) {
    lname <- nms[nms["code"] == lk[x, "link"], "link"]
    out <- paste0(
      names(lk[, "link"])[x], ": ", lname
    )
    if (lname %in% c("SIP", "polynomial")) {
      out <- paste0(
        out, " (a = ", lk[x, "ctt"], ", b = ", lk[x, "power"], ")"
      )
    }
    out
  }
  link <- c("g11", "g12")
  link <- object$link[link, ]
  es <- names(.current.error.scale)
  ans$link <- c(
    sapply(1:nrow(link), f, lk = link),
    paste0("g13: ", es[object$link["g13", "link"] == .current.error.scale], " scale")
  )
  if (endsWith(object$model, "V")) {
    link <- object$link[c("g2", "g21", "g22", "g23"), ]
    ans$link <- cbind(c(ans$link, ""), sapply(1:nrow(link), f, lk = link))
  }

  return(ans)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Print Method of class BTSR
#'
#' @description
#' Print method for objects of class `btsr`.
#'
#' @param x object of class `btsr`.
#'
#' @param digits  minimal number of significant digits, see [print.default].
#'
#' @param ... further arguments to be passed to or from other methods. They are
#'   ignored in this function
#'
#' @details
#' Users are not encouraged to call these internal functions directly.
#' Internal functions for package BTSR.
#'
#' @return Invisibly returns its argument, `x`.
#'
#' @importFrom stats coef
#'
#' @export
print.btsr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
      print.gap = 2L, quote = FALSE
    )
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function for printing the summary
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @rdname summary
#'
#' @importFrom stats quantile printCoefmat
#'
#' @param x an object of class `"summary.btsr"`, usually, a result of a call to
#'   `summary.btsr`.
#'
#' @param digits  minimal number of significant digits, see [print.default].
#'
#' @param signif.stars logical. If `TRUE`, \sQuote{significance stars} are
#'   printed for each coefficient.
#'
#' @details
#' `print.summary.btsr` tries to be smart about formatting the coefficients,
#' standard errors, etc. and additionally provides \sQuote{significance stars}.
#'
#' @export
print.summary.btsr <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"), ...) {
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat("\n-----------------------------------------------")
  cat("\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n",
    sep = ""
  )

  if (x$full.report) {
    cat("\nLink Functions:\n")
    .print.link(x$link)
    cat("-----------------------------------------------")
  }

  cat("\nResiduals:\n")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L) {
      structure(apply(t(resid), 1L, quantile),
        dimnames = list(nam, dimnames(resid)[[2L]])
      )
    } else {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  } else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  } else { # rdf == 0 : perfect fit!
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }
  if (x$df[1] == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    cat("\nCoefficients:\n")
    coefs <- x$coefficients

    printCoefmat(coefs,
      digits = digits,
      signif.stars = signif.stars,
      cs.ind = 1:2, # which columns are coefficients and standard errors
      tst.ind = 3, # which column is the test statistic
      P.values = TRUE,
      na.print = "NA", ...
    )
  }

  cat(
    "\nResidual standard error:",
    format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom"
  )

  if (x$full.report) {
    cat("\n-----------------------------------------------")
    cat("\n Covariance matrix:", x$vcov.type)
    cat("\n Optimization algorithm:", x$algorithm$method)
    cat("\n Number of function Evaluations:", x$algorithm$neval)
    cat("\n Convergence achieved:", x$algorithm$convergence)
    cat("\n-----------------------------------------------")
    cat("\n Log-likelihood:", format(x$loglik, digits = max(4L, digits + 1L)))
    cat("\n AIC:", format(x$aic, digits = max(4L, digits + 1L)))
    cat("\n BIC:", format(x$bic, digits = max(4L, digits + 1L)))
    cat("\n HQC:", format(x$hqc, digits = max(4L, digits + 1L)))
  }
  cat("\n-----------------------------------------------\n")
  invisible(x)
}


#' @title
#' Print Model Default Settings
#'
#' @description
#' Displays the default settings for a specified model in the BTSR package,
#' including link functions and their configurations.
#'
#' @param model Character string specifying the model name (e.g., "KREGV",
#'   "MARMA").
#'
#' @return Invisibly returns a list of data frames containing:
#' \itemize{
#'   \item \code{basic_info} - Model name and parent model (if different)
#'   \item \code{link_functions} - Link functions with their `ctt` and `power`
#'   parameters (for `"polynomial"` links)
#' }
#'
#' @examples
#' \dontrun{
#' # Print settings for KREGV model
#' BTSR.model.defaults("KREGV")
#'
#' # Print settings for MARMA model
#' BTSR.model.defaults("MARMA")
#' }
#'
#' @export
BTSR.model.defaults <- function(model) {
  model.settings <- .get.model.settings(model)

  # Create a list to hold all the table data
  tables <- list()

  if (endsWith(model, "V")) {
    model.settings$model <- paste0(model.settings$model, "V")
  }
  # Basic info table
  tables$basic_info <- data.frame(
    Model = model
  )
  if (model.settings$model != model) {
    tables$basic_info$Parent <- model.settings$model
  }

  link <- c("g11", "g12", "g2", "g21", "g22")
  model.settings$linkg <- model.settings$linkg[link]

  # Link function settings
  tables$link_functions <- data.frame(
    Link = names(model.settings$linkg),
    Function = unlist(model.settings$linkg),
    ctt = model.settings$lconfig[link, "ctt"],
    power = model.settings$lconfig[link, "power"]
  )
  linear <- tables$link_functions$Function != "polynomial"
  tables$link_functions[linear, c("ctt", "power")] <- "-"

  # Print all tables
  for (i in seq_along(tables)) {
    cat("\n", names(tables)[i], ":\n", sep = "")
    print(tables[[i]], row.names = FALSE)
    cat("\n")
  }

  invisible(tables)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Extracts the exact arguments (args) from "..." (mc)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.build.list <- function(mc, args) {
  configs <- lapply(args, function(name) {
    if (exists(name, where = mc)) {
      return(mc[[name, exact = TRUE]])
    } else {
      return(NULL)
    }
  })
  names(configs) <- args
  return(configs)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Check if two links are the equal
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.compare.link <- function(link, configs.linkg, key1, key2) {
  check <- link[[key1]] == link[[key2]]
  keys <- names(configs.linkg)
  if (is.null(keys)) {
    return(check)
  }
  for (key in keys) {
    check <- check && configs.linkg[[key]][[key1]] == configs.linkg[[key]][[key2]]
  }
  return(check)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Convert some values to pass to FORTRAN.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.update.control <- function(control) {
  # set default values
  control <- fit.control(control)

  # add extra information to pass to FORTRAN
  method <- substring(tolower(control$method), first = 1, last = 1)
  control$method.code <- switch(method,
    l = 0L,
    n = 1L,
    .stop.with.message(paste0(
      "\n Method ", control$method, " is not implemented.",
      "\n Available options are: `L-BFGS-B` or `Nelder-Mead`"
    ))
  )

  # vector of constants
  if (method == "l") {
    control$control2 <- c(factr = control$factr, pgtol = control$pgtol)
  } else {
    control$control2 <- c(stopcr = control$stopcr)
  }

  return(control)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Used to print information about the selected model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.print.call <- function(model, p, q, d, nreg) {
  if (model == "BARC") {
    if (sum(nreg) > 0) {
      model <- "BARCX"
    }
    return(paste0(model, "(", p[1], ") model"))
  }

  MODEL <- .get.model.type(model)

  if (MODEL == "iid") {
    return(paste0("i.i.d sample from a ", model, " distribution"))
  }

  if (MODEL == "REG") {
    dist <- .get.dist.name(model)
    return(paste0(dist, " regression"))
  }

  # ARMA-type model, ignore d
  if (any(endsWith(model, c("ARMA", "ARMAV")))) {
    dname <- NULL
  } else {
    dname <- sapply(d, function(x) ifelse(x, ",d", ",0"))
  }

  msg <- model
  # ARMAX/ARFIMAX-type model
  if (sum(nreg) > 0) msg <- paste0(msg, "X")

  msg <- paste0(msg, "(", p[1], dname[1], ",", q[1], ")")
  if (endsWith(model, "V")) msg <- paste0(msg, "x(", p[2], dname[2], ",", q[2], ")")
  msg <- paste0(msg, " model")
  return(msg)
}

.print.link <- function(link) {
  if (!is.matrix(link)) {
    # For single part
    cat(paste(link, collapse = ",   "), "\n")
  } else {
    # For two parts
    # Calculate column widths
    col1_width <- max(nchar(c("Part 1", link[, 1])))
    col2_width <- max(nchar(c("Part 2", link[, 2])))
    total_width <- col1_width + col2_width + 3

    # Create separator lines
    sep_line1 <- paste0(
      strrep("-", col1_width),
      "   ",
      strrep("-", col2_width)
    )

    # Print header
    cat(sep_line1, "\n", sep = "")
    cat(sprintf(
      paste0("%-", col1_width, "s   %-", col2_width, "s\n"),
      "Part 1 (P1)", "Part 2 (P2)"
    ))
    cat(sep_line1, "\n", sep = "")

    # Print each row
    for (i in 1:nrow(link)) {
      cat(sprintf(
        paste0("%-", col1_width, "s   %-", col2_width, "s\n"),
        link[i, 1], link[i, 2]
      ))
    }
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Compute the covariance matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.cov.matrix <- function(robust, outer, infoMatrix, gradient) {
  # check for null arguments
  minfo <- is.null(infoMatrix)
  mgrad <- is.null(gradient)

  # Nothing can be done if both are missing
  if (minfo && mgrad) {
    .stop.with.message(
      paste0(
        " Unable to proceed:",
        "\n Both the information matrix and the gradient are missing.",
        "\n Cannot calculate the covariance matrix.",
        "\n---------------------------------------------------------------",
        "\n Please fit the model again",
        "\n Set info = TRUE, to calculate the information matrix.",
        "\n To calculate the outer product (outer = TRUE) or",
        "\n the robust covariance matrix (robust = TRUE), set extra = TRUE"
      )
    )
  }

  # Use Hinv if outer = FALSE and robust = FALSE
  # back up the use options for future reference
  Hinv <- (!robust) && (!outer)
  cov.type <- c(Hinv = Hinv, robust = robust, outer = outer)

  if (minfo) {
    cat(
      "\n-------------------------------------------------",
      if (Hinv || cov.type["robust"]) "\n The information matrix is missing.",
      if (cov.type["robust"]) "\n Unable to compute the robust covariance matrix.",
      if (!cov.type["outer"]) "\n Using outer product instead.",
      "\n To calculate the information matrix, please",
      "\n fit the model again and set info = TRUE.",
      "\n-------------------------------------------------"
    )

    G <- t(gradient) %*% gradient
    cov.matrix <- try(solve(G))
    if (!("try.error" %in% class(cov.matrix))) {
      return(list(
        cov = cov.matrix,
        type = "Inverse of the outer product of the gradient"
      ))
    }
    .stop.with.message(
      paste0(
        " Failed inverting of outer product of the gradient",
        "\n No more options available.",
        "\n Try fitting the model again setting info = TRUE",
        "\n and use the robust matrix (robust = TRUE) or",
        "\n inverse of the information matrix (outer = FALSE, robust = FALSE)"
      )
    )
  }

  if (mgrad) {
    if (cov.type["robust"] || cov.type["outer"]) {
      cat(
        "\n-------------------------------------------------",
        if (cov.type["robust"] || cov.type["outer"]) "\n The gradient is missing.",
        if (cov.type["robust"]) "\n Unable to compute the robust covariance matrix.",
        if (cov.type["outer"]) "\n Unable to compute the outer product.",
        if (cov.type["robust"] || cov.type["outer"]) "\n Using the inverse of the information matrix instead.",
        "\n To compute the gradient, please",
        "\n fit the model again and set extra  = TRUE.",
        "\n-------------------------------------------------\n"
      )
    }
    cov.matrix <- try(solve(infoMatrix))
    if (!("try.error" %in% class(cov.matrix))) {
      return(list(
        cov = cov.matrix,
        type = "Inverse of the information matrix"
      ))
    }
    .stop.with.message(
      paste0(
        " Fail to calculate the inverse of the information matrix",
        "\n No more options available.",
        "\n Try fitting the model again setting extra = TRUE and",
        "\n use the outer product of the gradient (robust = TRUE)"
      )
    )
  }

  # information matrix and gradient are not missing
  # If robust = TRUE and outer = TRUE, first try robust
  if (robust && outer) outer <- FALSE

  if (robust || Hinv) {
    # find the inverse of the information matrix
    H.inv <- try(solve(infoMatrix))
    if (!("try-error" %in% class(H.inv))) {
      if (robust) {
        G <- t(gradient) %*% gradient
        return(list(
          cov = H.inv %*% G %*% H.inv,
          type = "Robust"
        ))
      }
      return(list(
        cov = H.inv,
        type = "Inverse of the information matrix"
      ))
    }
    # In case of error, use the outer product
    # If Hinv = TRUE, robust = FALSE and outer = FALSE
    cat(
      "\n------------------------------------------------------",
      "\n Fail to calculate the inverse of the information matrix",
      if (robust) "\n Unable to compute the robust covariance matrix.",
      "\n Using the outer product instead",
      "\n------------------------------------------------------\n"
    )
    Hinv <- FALSE
    outer <- TRUE
  }

  # find the inverse of the outer product
  cov.matrix <- try(solve(t(gradient) %*% gradient))
  if (!("try-error" %in% class(cov.matrix))) {
    return(list(
      cov = cov.matrix,
      type = "Inverse of the outer product of the gradient"
    ))
  }

  # CASE 1: outer = FALSE or robust = TRUE on input.
  # In this case the code gets here because Hinv failed
  if (!cov.type["outer"]) {
    .stop.with.message(
      " Failed inverting of outer product of the gradient. \n No more options available."
    )
  }

  # CASE 2: outer = TRUE and robust = FALSE on input.
  # In this Hinv was never evaluated
  cat(
    "\n---------------------------------------------------------",
    "\n Failed inverting of outer product of the gradient",
    "\n Using the inverse of the Hessin instead.",
    "\n---------------------------------------------------------\n"
  )
  cov.matrix <- try(solve(infoMatrix))
  if (!("try-error" %in% class(cov.matrix))) {
    return(list(
      cov = cov.matrix,
      type = "Inverse of the information matrix"
    ))
  }
  .stop.with.message(
    " Fail to calculate the inverse of the information matrix, \n No more options available."
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Compute the gradient
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.gradient <- function(object) {
  gradient <- NULL

  nd <- as.integer(sum(c(1 + (object$configs$npar[2] > 0), 1) * object$configs$npar))
  n <- length(object$series)

  D <- object$D
  T <- object$T
  h <- object$h

  # check if the gradient can be computed
  if (any(is.null(D), is.null(T), is.null(h)) ||
    any(is.nan(range(cbind(D, T, h)))) ||
    any(is.infinite(range(cbind(D, T, h))))) {
    return(NULL)
  }

  if (object$configs$extra == 1) {
    gradient <- try(.Fortran("gradient",
      n = n,
      npar = object$configs$npar,
      nd = nd,
      D = D,
      T = T,
      h = h,
      grad = matrix(0, nrow = n, ncol = max(1, sum(object$configs$npar))),
      PACKAGE = "BTSR"
    )$grad)
    if ("try-error" %in% class(gradient)) {
      return(NULL)
    }
    names(gradient) <- object$configs$coefsnames
  }

  return(gradient)
}
