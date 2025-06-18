# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               BTSR MODELS
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function:
# Interface between R and FORTRAN
# Also used to summarize the results of the simulation and return
# only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.btsr.sim <- function(model, ...) {
  # match the input arguments
  mc <- list(...)

  # create the configs object
  configs <- .get.model.configs(model, mc = mc, type = "sim")

  # call the FORTRAN subroutine
  N <- configs$n + configs$burn
  temp <- .Fortran("simbtsr",
    code = configs$code,
    length = c(burn = configs$burn, n = configs$n),
    order = configs$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = numeric(N),
      mut = numeric(N),
      eta1t = numeric(N),
      error1 = numeric(N),
      nut = numeric(N),
      varthetat = numeric(N),
      eta2t = numeric(N),
      error2 = numeric(N)
    ),
    xreg1 = configs$xreg$part1,
    xreg2 = configs$xreg$part2,
    tstart = c(
      ystart = configs$y.start,
      vtstart = configs$vt.start,
      e2start = configs$e2.start
    ),
    xstart = configs$xreg.start,
    link = configs$link[.current.link.names],
    lconfig = configs$lconfig[.current.link.names, ],
    np = length(configs$pdist),
    pdist = configs$pdist,
    xregar = configs$xregar,
    alpha = configs$alpha,
    beta = configs$beta,
    phi = configs$phi,
    theta = configs$theta,
    d = configs$d,
    conv = 0L, NAOK = TRUE, PACKAGE = "BTSR"
  )

  configs$complete <- ifelse(is.null(mc$complete), FALSE, mc$complete)
  invisible(
    .get.model.output(
      model = model, obj = temp, type = "sim",
      configs = configs, debug = mc$debug
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function:
# Interface between R and FORTRAN
# Also used to summarize the results of the extraction and return
# only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.btsr.extract <- function(model, ...) {
  # match the input arguments
  mc <- list(...)

  # create the configs object
  configs <- .get.model.configs(model, mc = mc, type = "extract")

  # call the FORTRAN subroutine
  nd <- as.integer(sum(c(1 + (configs$npar[2] > 0), 1) * configs$npar))
  temp <- .Fortran("extractbtsr",
    code = configs$code,
    length = c(n = configs$n, nnew = configs$nnew),
    order = configs$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = mc$yt,
      g11yt = numeric(configs$n),
      g12yt = numeric(configs$n),
      mut = numeric(configs$n),
      eta1t = numeric(configs$n),
      error1 = numeric(configs$n),
      nut = numeric(configs$n),
      varthetat = numeric(configs$n),
      eta2t = numeric(configs$n),
      g22varthetat = numeric(configs$n),
      error2 = numeric(configs$n)
    ),
    xreg1 = configs$xreg$part1,
    xreg2 = configs$xreg$part2,
    ts.start = c(
      y.start = configs$y.start,
      vt.start = configs$vt.start,
      e2.start = configs$e2.start
    ),
    xstart = configs$xreg.start,
    xnew1 = configs$xnew$part1,
    xnew2 = configs$xnew$part2,
    forecast = cbind(
      mut = numeric(max(1, configs$nnew)),
      eta1t = numeric(max(1, configs$nnew)),
      nut = numeric(max(1, configs$nnew)),
      varthetat = numeric(max(1, configs$nnew)),
      eta2t = numeric(max(1, configs$nnew))
    ),
    link = configs$linkg[.current.link.names],
    lconfig = configs$lconfig[.current.link.names, ],
    npar = configs$npar,
    coefs = configs$coefs,
    xregar = configs$xregar,
    nfix = configs$nfix[, c("alpha", "beta", "phi", "theta", "d")],
    alpha = configs$alpha$fvalues,
    flagsb = configs$beta$flags,
    beta = configs$beta$fvalues,
    flagsphi = configs$phi$flags,
    phi = configs$phi$fvalues,
    flagstheta = configs$theta$flags,
    theta = configs$theta$fvalues,
    d = configs$d$fvalues,
    np = length(configs$pdist),
    pdist = configs$pdist,
    extras = c(
      m = configs$m,
      llk = configs$llk,
      sco = configs$sco,
      info = configs$info,
      extra = configs$extra
    ),
    sll = 0,
    U = numeric(max(1, sum(configs$npar) * configs$sco)),
    K = diag(max(1, sum(configs$npar) * configs$info)),
    nd = nd,
    D = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, nd * configs$extra)
    ),
    T = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, 2 * configs$extra)
    ),
    E = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, 3 * configs$extra)
    ),
    h = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, 2 * configs$extra)
    ),
    conv = 0L, NAOK = TRUE, PACKAGE = "BTSR"
  )

  invisible(
    .get.model.output(
      model = model, obj = temp, type = "extract",
      configs = configs, debug = mc$debug
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function:
# Interface between R and FORTRAN
# Also used to summarize the results of the optimization
# Returns only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.btsr.fit <- function(model, ...) {
  # match the input arguments
  mc <- list(...)

  report <- mc$report
  if (is.null(report)) report <- TRUE
  if (report) mc$info <- TRUE
  mc$llk <- TRUE
  full.report <- mc$full.report
  if (is.null(full.report)) full.report <- report

  # create the configs object
  configs <- .get.model.configs(model, mc = mc, type = "fit")

  # call the FORTRAN subroutine
  nd <- as.integer(sum(c(1 + (configs$npar[2] > 0), 1) * configs$npar))
  temp <- .Fortran("optimbtsr",
    method = configs$control$method.code,
    code = configs$code,
    length = c(n = configs$n, nnew = configs$nnew),
    order = configs$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = mc$yt,
      g11yt = numeric(configs$n),
      g12yt = numeric(configs$n),
      mut = numeric(configs$n),
      eta1t = numeric(configs$n),
      error1 = numeric(configs$n),
      nut = numeric(configs$n),
      varthetat = numeric(configs$n),
      eta2t = numeric(configs$n),
      g22varthetat = numeric(configs$n),
      error2 = numeric(configs$n)
    ),
    xreg1 = configs$xreg$part1,
    xreg2 = configs$xreg$part2,
    ts.start = c(
      y.start = configs$y.start,
      vt.start = configs$vt.start,
      e2.start = configs$e2.start
    ),
    xstart = configs$xreg.start,
    xnew1 = configs$xnew$part1,
    xnew2 = configs$xnew$part2,
    forecast = cbind(
      mut = numeric(max(1, configs$nnew)),
      eta1t = numeric(max(1, configs$nnew)),
      nut = numeric(max(1, configs$nnew)),
      varthetat = numeric(max(1, configs$nnew)),
      eta2t = numeric(max(1, configs$nnew))
    ),
    link = configs$linkg[.current.link.names],
    lconfig = configs$lconfig[.current.link.names, ],
    npar = configs$npar,
    coefs = configs$coefs,
    nbd = configs$nbd,
    bounds = cbind(lower = configs$lower, upper = configs$upper),
    xregar = configs$xregar,
    nfix = configs$nfix[, c("alpha", "beta", "phi", "theta", "d")],
    alpha = configs$alpha$fvalues,
    flagsb = configs$beta$flags,
    beta = configs$beta$fvalues,
    flagsphi = configs$phi$flags,
    phi = configs$phi$fvalues,
    flagstheta = configs$theta$flags,
    theta = configs$theta$fvalues,
    d = configs$d$fvalues,
    np = length(configs$pdist),
    pdist = configs$pdist,
    extras = c(
      m = configs$m,
      llk = configs$llk,
      sco = configs$sco,
      info = configs$info,
      extra = configs$extra
    ),
    sll = 0,
    U = numeric(max(1, sum(configs$npar) * configs$sco)),
    K = diag(max(1, sum(configs$npar) * configs$info)),
    nd = nd,
    D = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, nd * configs$extra)
    ),
    T = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, 2 * configs$extra)
    ),
    E = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, 3 * configs$extra)
    ),
    h = matrix(0,
      nrow = max(1, configs$n * configs$extra),
      ncol = max(1, 2 * configs$extra)
    ),
    control1 = c(
      iprint = as.integer(configs$control$iprint),
      maxit = as.integer(configs$control$maxit)
    ),
    nc2 = length(configs$control$control2),
    control2 = configs$control$control2,
    neval = 0L,
    conv = 0L, NAOK = TRUE, PACKAGE = "BTSR"
  )

  out <- .get.model.output(
    model = model, obj = temp, type = "fit",
    configs = configs, debug = mc$debug
  )
  if (report) print(summary(out, ...))

  invisible(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function:
# Interface between R and FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.btsr.predict <- function(obj, newdata, nnew, debug) {
  # get configurations for prediction
  object <- .get.predict.configs(obj = obj, newdata = newdata, nnew = nnew)

  # call FORTRAN
  temp <- .Fortran("predictbtsr",
    length = c(n = object$n, nnew = object$nnew),
    order = object$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = object$series,
      g11yt = object$gyt[, "g11yt"],
      g12yt = object$gyt[, "g12yt"],
      error1 = object$error[, "error1"],
      varthetat = object$fitted.values[, "varthetat"],
      eta2t = object$etat[, "eta2t"],
      g22varthetat = object$g22varthetat,
      error2 = object$error[, "error2"]
    ),
    xreg1 = object$xreg$part1,
    xreg2 = object$xreg$part2,
    xnew1 = object$xnew$part1,
    xnew2 = object$xnew$part2,
    forecast = cbind(
      mut = numeric(max(1, object$nnew)),
      eta1t = numeric(max(1, object$nnew)),
      nut = numeric(max(1, object$nnew)),
      varthetat = numeric(max(1, object$nnew)),
      eta2t = numeric(max(1, object$nnew))
    ),
    link = object$linkg[.current.link.names],
    lconfig = object$lconfig[.current.link.names, ],
    npar = object$npar,
    coefs = object$coefs,
    xregar = object$xregar,
    nfix = object$nfix,
    alpha = object$alpha$fvalues,
    flagsb = object$beta$flags,
    beta = object$beta$fvalues,
    flagsphi = object$phi$flags,
    phi = object$phi$fvalues,
    flagstheta = object$theta$flags,
    theta = object$theta$fvalues,
    d = object$d$fvalues, NAOK = TRUE, PACKAGE = "BTSR"
  )

  nm <- .get.output.names(model = object$model)
  # predicted values
  out <- list(forecast = temp$forecast[, nm$forecast])
  # regressors
  if (sum(object$order[, "nreg"]) > 0) {
    if (length(nm$xnew) == 1) {
      out$xnew <- obj$xnew1
    } else {
      w <- object$order[, "nreg"] > 0
      out$xnew[.parts[w]] <- as.list(temp[nm$xnew[w]])
    }
  }

  if (debug) out$out.Fortran <- temp
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               BARC MODELS
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function
# Simulates BARC models
# Makes the calculations and reports only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.barc.sim <- function(...) {
  mc <- list(...)

  # create the configs object
  configs <- .get.model.configs(model = "BARC", mc = mc, type = "sim")

  temp <- .Fortran("simbarc",
    length = c(n = configs$n, burn = configs$burn),
    order = configs$order[1, c("nreg", "p", "q")],
    ts = cbind(
      yt = numeric(configs$n + configs$burn),
      mut = numeric(configs$n + configs$burn),
      eta1t = numeric(configs$n + configs$burn),
      error1 = numeric(configs$n + configs$burn),
      Tt = numeric(configs$n + configs$burn)
    ),
    xreg = configs$xreg$part1,
    ystart = configs$y.start,
    xstart = configs$xreg.start[1, ],
    linkg = configs$linkg[c("g11", "g12", "g13", "h")],
    lconfig = configs$lconfig[c("g11", "g12", "g13", "h"), ],
    map = configs$map,
    xregar = configs$xregar[1],
    alpha = configs$alpha,
    beta = configs$beta[1, ],
    phi = configs$phi[1, ],
    theta = configs$theta[1, ],
    u0 = configs$u0[1],
    conv = 0L, NAOK = TRUE, PACKAGE = "BTSR"
  )

  configs$complete <- ifelse(is.null(mc$complete), FALSE, mc$complete)
  invisible(
    .get.model.output(
      model = "BARC", obj = temp, type = "sim",
      configs = configs, debug = mc$debug
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function
# Extracts components from BARC models
# Makes the calculations and reports only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.barc.extract <- function(...) {
  mc <- list(...)
  mc$inf <- 0
  mc$extra <- FALSE

  # create the configs object
  configs <- .get.model.configs(model = "BARC", mc = mc, type = "extract")

  temp <- .Fortran("extractbarc",
    length = c(n = configs$n, nnew = configs$nnew),
    order = configs$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = mc$yt,
      g11yt = numeric(configs$n),
      g12yt = numeric(configs$n),
      mut = numeric(configs$n),
      eta1t = numeric(configs$n),
      error1 = numeric(configs$n),
      Tt = numeric(configs$n)
    ),
    xreg1 = configs$xreg$part1,
    ystart = configs$y.start,
    xstart = configs$xreg.start,
    xnew1 = configs$xnew$part1,
    forecast = cbind(
      mut = numeric(max(1, configs$nnew)),
      eta1t = numeric(max(1, configs$nnew)),
      Tt = numeric(max(1, configs$nnew))
    ),
    link = configs$linkg[c(.current.link.names, "h")],
    lconfig = configs$lconfig[c(.current.link.names, "h"), ],
    map = configs$map,
    npar = configs$npar,
    coefs = configs$coefs,
    xregar = configs$xregar,
    nfix = configs$nfix[, c("alpha", "beta", "phi", "theta", "u0")],
    alpha = configs$alpha$fvalues,
    flagsb = configs$beta$flags,
    beta = configs$beta$fvalues,
    flagsphi = configs$phi$flags,
    phi = configs$phi$fvalues,
    flagstheta = configs$theta$flags,
    theta = configs$theta$fvalues,
    u0 = configs$u0$fvalues[1],
    extras = c(
      m = 0L,
      llk = configs$llk,
      sco = configs$sco,
      info = configs$info,
      extra = 0L
    ),
    sll = 0,
    U = numeric(max(1, sum(configs$npar) * configs$sco)),
    K = diag(max(1, sum(configs$npar) * configs$info)),
    conv = 0L, NAOK = TRUE, PACKAGE = "BTSR"
  )

  invisible(
    .get.model.output(
      model = "BARC", obj = temp, type = "extract",
      configs = configs, debug = mc$debug
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function
# Fit a BARC model
# Makes the calculations and reports only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.barc.fit <- function(...) {
  # match the input arguments
  mc <- list(...)
  mc$inf <- 0

  report <- mc$report
  if (is.null(report)) report <- TRUE
  if (report) mc$info <- TRUE
  mc$llk <- TRUE
  full.report <- mc$full.report
  if (is.null(full.report)) full.report <- report

  # create the configs object
  configs <- .get.model.configs("BARC", mc = mc, type = "fit")

  # call the FORTRAN subroutine
  nd <- as.integer(sum(c(1 + (configs$npar[2] > 0), 1) * configs$npar))
  temp <- .Fortran("optimbarc",
    method = configs$control$method.code,
    length = c(n = configs$n, nnew = configs$nnew),
    order = configs$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = mc$yt,
      g11yt = numeric(configs$n),
      g12yt = numeric(configs$n),
      mut = numeric(configs$n),
      eta1t = numeric(configs$n),
      error1 = numeric(configs$n),
      Tt = numeric(configs$n)
    ),
    xreg1 = configs$xreg$part1,
    ystart = configs$y.start,
    xstart = configs$xreg.start,
    xnew1 = configs$xnew$part1,
    forecast = cbind(
      mut = numeric(max(1, configs$nnew)),
      eta1t = numeric(max(1, configs$nnew)),
      Tt = numeric(max(1, configs$nnew))
    ),
    link = configs$linkg[c(.current.link.names, "h")],
    lconfig = configs$lconfig[c(.current.link.names, "h"), ],
    map = configs$map,
    npar = configs$npar,
    coefs = configs$coefs,
    nbd = configs$nbd,
    bounds = cbind(lower = configs$lower, upper = configs$upper),
    xregar = configs$xregar,
    nfix = configs$nfix[, c("alpha", "beta", "phi", "theta", "u0")],
    alpha = configs$alpha$fvalues,
    flagsb = configs$beta$flags,
    beta = configs$beta$fvalues,
    flagsphi = configs$phi$flags,
    phi = configs$phi$fvalues,
    flagstheta = configs$theta$flags,
    theta = configs$theta$fvalues,
    u0 = configs$u0$fvalues[1],
    extras = c(
      m = 0L,
      llk = configs$llk,
      sco = configs$sco,
      info = configs$info,
      extra = 0L
    ),
    sll = 0,
    U = numeric(max(1, sum(configs$npar) * configs$sco)),
    K = diag(max(1, sum(configs$npar) * configs$info)),
    control1 = c(
      iprint = as.integer(configs$control$iprint),
      maxit = as.integer(configs$control$maxit)
    ),
    nc2 = length(configs$control$control2),
    control2 = configs$control$control2,
    neval = 0L,
    conv = 0L, NAOK = TRUE, PACKAGE = "BTSR"
  )

  out <- .get.model.output(
    model = "BARC", obj = temp, type = "fit",
    configs = configs, debug = mc$debug
  )
  if (report) print(summary(out, ...))
  invisible(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function
# Prediction for BARC models
# Makes the calculations and reports only the relevant variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.barc.predict <- function(obj, newdata, nnew, debug) {
  if (!"barc" %in% class(obj)) {
    .stop.with.message("The object is not a 'BARC' model")
  }

  # get configurations for prediction
  object <- .get.predict.configs(obj = obj, newdata = newdata, nnew = nnew)

  # call FORTRAN
  temp <- .Fortran("predictbarc",
    length = c(n = object$n, nnew = object$nnew),
    order = object$order[, c("nreg", "p", "q", "inf")],
    ts = cbind(
      yt = object$series,
      g11yt = object$gyt[, "g11yt"],
      g12yt = object$gyt[, "g12yt"],
      error1 = object$error,
      Tt = object$fitted.values[, "Tt"]
    ),
    xreg1 = object$xreg$part1,
    xnew1 = object$xnew$part1,
    forecast = cbind(
      mut = numeric(max(1, object$nnew)),
      eta1t = numeric(max(1, object$nnew)),
      Tt = numeric(max(1, object$nnew))
    ),
    link = object$linkg[c(.current.link.names, "h")],
    lconfig = object$lconfig[c(.current.link.names, "h"), ],
    map = object$map,
    npar = object$npar,
    coefs = object$coefs,
    xregar = object$xregar,
    nfix = object$nfix,
    alpha = object$alpha$fvalues,
    flagsb = object$beta$flags,
    beta = object$beta$fvalues,
    flagsphi = object$phi$flags,
    phi = object$phi$fvalues,
    flagstheta = object$theta$flags,
    theta = object$theta$fvalues,
    u0 = object$u0$fvalues, NAOK = TRUE, PACKAGE = "BTSR"
  )

  nm <- .get.output.names(model = object$model)
  # predicted values
  out <- list(forecast = temp$forecast[, nm$forecast])
  # regressors
  if (sum(object$order[, "nreg"]) > 0) {
    if (length(nm$xnew) == 1) {
      out$xnew <- obj$xnew1
    } else {
      w <- object$order[, "nreg"] > 0
      out$xnew[.parts[w]] <- as.list(temp[nm$xnew[w]])
    }
  }

  if (debug) out$out.Fortran <- temp
  return(out)
}
