# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Returns the default values for m, llk, sco, info, extra
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.extra.configs <- function(extra.configs) {
  # default values
  default <- .default.extra.configs

  # update the default values using user provided values
  for (key in names(default)) {
    if (!is.null(extra.configs[[key]])) {
      default[[key]] <- extra.configs[[key]]
    }
  }
  # convert to FORTRAN format
  default <- lapply(default, as.integer)
  return(default)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function checks and converts the error.code to pass to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.escale <- function(escale) {
  check <- escale == .current.error.scale
  if (any(check)) {
    return(as.integer(escale))
  }
  .stop.with.message(.scale.message)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Initializes y.start and xreg.start to pass to FORTRAN.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.data.start <- function(linkg, data.start, nreg) {
  out <- c()

  # starting values are provided, uses this values
  # otherwise, the initial values are set as the default value
  out$y.start <- ifelse(
    is.null(data.start$y.start),
    .get.link.starting.values(linkg$g12),
    data.start$y.start
  )
  out$vt.start <- ifelse(
    is.null(data.start$vt.start),
    .get.link.starting.values(linkg$g22),
    data.start$vt.start
  )
  out$e2.start <- ifelse(
    is.null(data.start$e2.start),
    0,
    data.start$e2.start
  )
  out$xreg.start <- matrix(0, nrow = 2, ncol = max(1, nreg))

  # update only if necessary
  if (nreg[1] > 0) {
    out$xreg.start[1, 1:nreg[1]] <- data.start$xreg.start$part1
  }
  if (nreg[2] > 0) {
    out$xreg.start[2, 1:nreg[2]] <- data.start$xreg.start$part2
  }

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function convert xreg to the format required by FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.xreg.part <- function(xreg, n, part) {
  if (is.null(xreg)) {
    #  nreg = 0
    xreg <- matrix(0, ncol = 1, nrow = max(1, n))
    nreg <- 0L
  } else {
    #  nreg > 0
    xreg <- as.matrix(xreg)
    nreg <- ncol(xreg)
    if (nrow(xreg) != n) {
      .stop.with.message(
        paste0(
          " In part ", part, " of the model:",
          "\n xreg and y do not have the same number of observations"
        )
      )
    }
  }
  return(list(xreg = xreg, nreg = nreg))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function convert xnew to the format required by FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.xnew.part <- function(xnew, nnew, nreg, part, skip.forecast) {
  nnew <- as.integer(max(0, nnew))

  # set default values if
  #  - if nnew = 0 => ncol = max(1, nreg)
  #  - forecast is not required (for simulation purposes) => ncol = 1
  if (skip.forecast) {
    xnew <- matrix(0, nrow = 1, ncol = max(1, nreg))
    nnew <- 0L
    return(list(xnew = xnew, nnew = nnew))
  }

  # no forecast required or nreg = 0
  #  - create a dummy variable
  if (nnew * nreg == 0) {
    xnew <- matrix(0, nrow = max(1, nnew), ncol = max(1, nreg))
    return(list(xnew = xnew, nnew = nnew))
  }

  # if forecast is required
  #  - if the code gets here nreg > 0
  if (is.null(xnew)) {
    .stop.with.message(" New values for `xreg` are missing with no default. ")
  } else {
    xnew <- as.matrix(xnew)
    nnew <- min(nnew, nrow(xnew))
    if (nnew < nrow(xnew)) xnew <- xnew[1:nnew, ]
    if (nnew == 1) xnew <- matrix(xnew, nrow = 1)

    # chec if xnew is compatible with xreg
    if (ncol(xnew) != nreg) {
      .stop.with.message(
        paste0(
          " In part ", part, " of the model:",
          "\n The number of columns in xnew and xreg are not the same"
        )
      )
    }
  }
  return(list(xnew = xnew, nnew = nnew))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Initializes xreg and xnew to pass to FORTRAN.
# This function creates dummy matrices of size 1 x 1 when
# the model does not have regressors or forecast is not required.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.xreg <- function(model, xreg, xnew, n, nnew, skip.forecast) {
  # for iid samples, ignore the regressors
  if (.get.model.type(model) == "iid") xreg <- xnew <- NULL

  # set default values
  if (is.null(nnew)) nnew <- 0
  out <- list(
    nnew = as.integer(nnew),
    nreg <- c(0L, 0L)
  )

  # For compatibility with previous versions
  out$xreg <- .check.xreg.format(xreg)
  out$xnew <- .check.xreg.format(xnew)
  if (!endsWith(model, "V")) {
    out$xreg$part2 <- NULL
    out$xnew$part2 <- NULL
  }

  # convert xreg an xnew to FORTRAN format
  for (i in 1:length(.parts)) {
    key <- .parts[i]
    temp <- .convert.xreg.part(xreg = out$xreg[[key]], n = n, part = i)
    out$xreg[[key]] <- temp$xreg
    out$nreg[i] <- temp$nreg

    temp <- .convert.xnew.part(
      xnew = out$xnew[[key]], nnew = out$nnew, nreg = out$nreg[i],
      part = i, skip.forecast = skip.forecast
    )
    out$xnew[[key]] <- temp$xnew
    out$nnew <- temp$nnew
  }
  if (nrow(out$xnew$part1) > out$nnew) out$xnew$part1 <- out$xnew$part1[1:out$nnew, ]

  return(out)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function updates values, if provided or set to zero (if NULL)
# Also, converts to the required FORTRAN format
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.coefs.sim <- function(coefs, keys, npars, npar.key, scalar) {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # helper function
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # replaces NULL by zero
  .update.coefs.each <- function(par1, par2, npar, scalar) {
    if (scalar) {
      # par1 and par2 are scalar return a vector
      if (is.null(par1)) par1 <- 0
      if (is.null(par2)) par2 <- 0
      par <- c(par1, par2)
    } else {
      # par1 and par2 are vectors return a matrix
      par <- matrix(0, nrow = 2, ncol = max(1, npar))
      if (npar[1] > 0) par[1, 1:npar[1]] <- par1
      if (npar[2] > 0) par[2, 1:npar[2]] <- par2
    }
    return(par)
  }

  out <- list(npar = c(0, 0))
  # loop over the coefficients that need update
  for (i in 1:length(keys)) {
    # which parameter needs update
    key <- keys[i]
    if (scalar[i]) {
      npar <- c(length(coefs$part1[[key]]), length(coefs$part2[[key]]))
    } else {
      npar <- npars[, npar.key[i]]
    }
    #  convert to a vector with the values from part 1 and part 2
    out[[key]] <- .update.coefs.each(
      par1 = coefs$part1[[key]],
      par2 = coefs$part2[[key]],
      npar = npar,
      scalar = scalar[i]
    )
    out$npar <- out$npar + npar
  }
  out$npar <- as.integer(out$npar)

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Not supposed to be called by the user.
# Initializes the variables with fixed values and fixed lags to pass to FORTRAN.
# Used only for extraction and to fit a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.coefs.each <- function(parname, fvalues, flags, coefs, lags, npar, part) {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # helper functions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # set fit and fixed values or only fit
  .set.fit <- function(nop, fvl, flg, nf, cf, lg, pname) {
    cfnames <- NULL

    # check the total number of parameters
    if (nop == 0) {
      cf <- NULL
      nf <- 0
    } else {
      # check if there exist any non-fixed parameter
      if (nf < nop) {
        if (!(pname %in% c("alpha", "d", "u0"))) {
          cfnames <- paste(pname, "(", lg, ")", sep = "")
        } else {
          cfnames <- pname
        }
        if (is.null(cf)) cf <- rep(0, nop)
      }
    }

    # check fixed values and lags
    if (nf == 0) {
      fvl <- 0
      flg <- 0L
      nf <- 0L
    } else {
      if (is.null(fvl)) fvl <- rep(0, max(nf))
      if (nf == nop && is.null(flg)) flg <- 1:nop
    }

    return(list(
      fvalues = fvl,
      flags = as.integer(flg),
      nfix = as.integer(nf),
      coefs = setNames(cf, cfnames),
      coefsnames = cfnames
    ))
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # helper function. Error message
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # print an error message to help debug
  .input.error <- function(msg, inpar = NULL,
                           pcoefs = FALSE, pnpar = FALSE,
                           plags = FALSE, pfvalues = FALSE, pflags = FALSE,
                           icoefs = NULL, ilc = NULL, ilags = NULL, ill = NULL,
                           ifvalues = NULL, ilfv = NULL, iflags = NULL, ilfl = NULL) {
    is.beta <- parname == "beta"
    paste0(
      " For parameter ", parname, " in part ", part, ":",
      "\n  - ", msg,
      "\n\n Based on the input:",
      if (pnpar && !is.beta) paste0("\n  - number of parameters: ", inpar),
      if (pnpar && is.beta) paste0("\n  - number of regressors: ", inpar),
      if (pcoefs) {
        paste0(
          "\n  - length(coefs) = ", ilc,
          ifelse(ilc > 0, ". Coefficients provided: ", ".")
        )
      },
      if (pcoefs && ilc > 0) paste0("\n    ", paste0(icoefs, collapse = ", ")),
      if (plags) {
        paste0(
          "\n  - length(lags) = ", ill,
          ifelse(ill > 0, ". Lags provided: ", ".")
        )
      },
      if (plags && ill > 0) paste0("\n    ", paste0(ilags, collapse = ", ")),
      if (pfvalues) {
        paste0(
          "\n  - length(fixed.values) = ", ilfv,
          ifelse(ilfv > 0, ". Fixed values provided: ", ".")
        )
      },
      if (pfvalues && ilfv > 0) paste0("\n    ", paste0(ifvalues, collapse = ", ")),
      if (pflags) {
        paste0(
          "\n  - length(fixed.lags) = ", ilfl,
          ifelse(ilfl > 0, ". Fixed lags provided: ", ".")
        )
      },
      if (pflags && ilfl > 0) paste0("\n    ", paste0(iflags, collapse = ", ")),
      "\n\n Please check if all fixed values/lags and",
      "\n all non-fixed values/lags were correctly informed.",
      if (is.beta) "\n Also check if the regressors have the correct dimension."
    )
  }

  # check inputs
  coefs.null <- is.null(coefs)
  lags.null <- is.null(lags)
  fvalues.null <- is.null(fvalues)
  flags.null <- is.null(flags)


  # SPECIAL CASE:
  # CASE 1: all variables are NULL
  #   Set all parameters as fixed and equal to zero.
  if (coefs.null && lags.null && fvalues.null && flags.null) {
    return(.set.fit(
      nop = npar, fvl = NULL, flg = NULL, nf = npar,
      cf = NULL, lg = NULL, pname = parname
    ))
  }

  lc <- length(coefs)
  ll <- length(lags)
  lfv <- length(fvalues)
  lfl <- length(flags)

  if (npar > 0) {
    # set number of fixed and non-fixed values
    nfit <- max(lc, ll)
    nfix <- max(lfv, lfl)
    if (nfit == npar && lags.null) lags <- 1:npar
    if (nfix == npar && flags.null) flags <- 1:npar
  }

  # Basic checks:
  # lags and fixed lags cannot have any intersection
  if (any(lags %in% flags)) {
    .stop.with.message(
      .input.error(
        msg = "lags and fixed.lags have non-empty intersection. ",
        plags = TRUE, pflags = TRUE, ilags = lags, ill = ll,
        iflags = flags, ilfl = lfl
      )
    )
  }

  # setting some default values
  # lags/flags missing
  if (npar == 1) {
    if (lags.null && !coefs.null) {
      lags <- 1
      ll <- 1
    }
    if (flags.null && !fvalues.null) {
      flags <- 1
      lfl <- 1
    }
  } else {
    # if lags are not provided, they must complement fixed lags
    if (lags.null && !flags.null && lfl < npar) {
      lags <- c(1:npar)[-flags]
      ll <- npar - lfl
      lags.null <- FALSE
    }
    # if flags are not provided, they must complement lags
    if (!lags.null && flags.null && ll < npar) {
      flags <- c(1:npar)[-lags]
      lfl <- npar - ll
      flags.null <- FALSE
    }
  }
  # coefs/fvalues missing
  if (coefs.null && !lags.null) {
    coefs <- rep(0, ll)
    lc <- ll
    coefs.null <- FALSE
  }
  if (fvalues.null && !flags.null) {
    fvalues <- rep(0, lfl)
    lfv <- lfl
    fvalues.null <- FALSE
  }

  # set number of fixed and non-fixed values
  nfit <- max(lc, ll)
  nfix <- max(lfv, lfl)


  # Basic check: size of the vectors
  if (nfit + nfix != npar) {
    msg <- "values provided are not compatible. "
    if (ll + lfl != npar) {
      msg <- paste0(
        msg, "\n  - unable to decide which lags must be fixed. "
      )
    }
    if (lc + lfv != npar) {
      ppname <- switch(parname,
        "beta" = "\n    'nreg' (the number of regressors).",
        "  `npar`. "
      )
      msg <- paste0(
        msg, "\n  - 'coefs' and 'fixed.values' size do not match",
        ppname
      )
    }
    .stop.with.message(
      .input.error(msg,
        pnpar = TRUE, plags = TRUE, pflags = TRUE,
        pcoefs = TRUE, pfvalues = TRUE, inpar = npar,
        ilags = lags, ill = ll, iflags = flags, ilfl = lfl,
        icoefs = coefs, ilc = lc, ifvalues = fvalues, ilfv = lfv
      )
    )
  }

  # Basic check: coefs + fvalues and lags and/or flags
  if (!coefs.null && !fvalues.null && lags.null && flags.null) {
    .stop.with.message(
      .input.error(
        paste0(
          "`lags` and `fixed.lags` are missing",
          "\n  - unable to decide which lags must be fixed"
        ),
        pcoefs = TRUE, pfvalues = TRUE, icoefs = coefs, ilc = lc,
        ifvalues = fvalues, ilfv = lfv
      )
    )
  }

  # if non-fixed values and lags are provided, they must have the same size
  if ((!coefs.null && !lags.null) && (lc != ll)) {
    .stop.with.message(
      .input.error("'coefs' and 'lags' are not compatible. ",
        plags = TRUE, pcoefs = TRUE, icoefs = coefs, ilc = lc,
        ilags = lags, ill = ll
      )
    )
  }

  # if fixed values and flags are provided, they must have the same size
  if ((!fvalues.null && !flags.null) && (lfv != lfl)) {
    .stop.with.message(
      .input.error("'fixed.values' and 'fixed.lags' are not compatible. ",
        pfvalues = TRUE, pflags = TRUE, ifvalues = fvalues, ilfv = lfv,
        iflags = flags, ilfl = lfl
      )
    )
  }

  # SPECIAL CASES:
  # CASE 2: all variables are fixed
  if (nfix == npar) {
    return(.set.fit(
      nop = npar, fvl = fvalues, flg = NULL, nf = npar,
      cf = NULL, lg = NULL, pname = parname
    ))
  }
  # CASE 3: all variables are non-fixed
  if (nfit == npar) {
    return(.set.fit(
      nop = npar, fvl = 0, flg = 0, nf = 0, cf = coefs, lg = lags, pname = parname
    ))
  }

  # fixed values
  if (is.null(fvalues)) fvalues <- rep(0, lfl)
  return(.set.fit(
    nop = npar, fvl = fvalues, flg = flags, nf = lfl, cf = coefs,
    lg = lags, pname = parname
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Not supposed to be called by the user.
# Initializes the variables with fixed values and fixed lags for each
# part of the model.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.coefs.part <- function(model, coefs, lags, fixed.values,
                                fixed.lags, order, part) {
  out <- c()
  par <- NULL
  nm <- NULL

  cf.names <- .names.coefs(type = model)
  npar <- c(
    alpha = 1,
    beta = order[["nreg"]],
    phi = order[["p"]],
    theta = order[["q"]],
    d = 1,
    u0 = 1
  )[cf.names]

  for (i in 1:length(cf.names)) {
    key <- cf.names[i]
    out[[key]] <- .convert.coefs.each(
      parname = key, fvalues = fixed.values[[key]],
      flags = fixed.lags[[key]], coefs = coefs[[key]],
      lags = lags[[key]], npar = npar[i], part = part
    )
    par <- c(par, out[[key]]$coefs)
    nm <- c(nm, out[[key]]$coefsname)
  }

  if (!is.null(par)) names(par) <- nm
  out$coefs <- par
  out$coefsnames <- nm

  if (model == "BARC") {
    out$coefsnames[out$coefsnames == "theta(1)"] <- "theta"
    names(out$coefs) <- out$coefsnames
  }
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Not supposed to be called by the user.
# Initializes the variables with fixed values and fixed lags
# to pass to FORTRAN. Merges the information from part 1 and part 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.coefs <- function(model, coefs, lags, fixed.values, fixed.lags,
                           order, is.sim = FALSE) {
  # for simulation, only coefs are required. Ignore all other
  # for extraction and fitting a model, make a full conversion.

  if (is.sim) {
    # updating the coefficients values
    temp <- .need.update(model)
    return(.convert.coefs.sim(
      coefs = coefs, keys = temp$need.update, npars = order,
      npar.key = temp$npar.key, scalar = temp$scalar
    ))
  }

  # - extraction or fitting -
  # Converting each part to FORTRAN format
  parts <- c()
  nparts <- length(.parts)
  MODEL <- ifelse(model == "BARC", model, "arfima")
  for (i in 1:nparts) {
    key <- .parts[i]
    parts[[key]] <- .convert.coefs.part(
      model = MODEL,
      coefs = coefs[[key]],
      lags = lags[[key]],
      fixed.values = fixed.values[[key]],
      fixed.lags = fixed.lags[[key]],
      order = order[i, ], part = i
    )
  }
  out <- list(npar = c(length(parts$part1$coefs), length(parts$part2$coefs)))

  # check if part 2 is ok
  onlymu <- .get.base.model(model) %in% .current.models("nonu")
  if (!onlymu && out$npar[2] == 0) {
    if (is.null(fixed.values$part2)) {
      msg <- ifelse(endsWith(model, "V"), " Parameters in part 2 are ", " 'nu' is ")
      .stop.with.message(paste0(msg, "missing with no default"))
    }
  }

  # for compatibility with previous versions
  if (!endsWith(model, "V")) {
    if (out$npar[2] > 1) {
      .stop.with.message(
        paste0(
          " Unknown parameters in part 2",
          "\n For model ", model, " part 2 must include only 'alpha'"
        )
      )
    }
    if (out$npar[2] == 1) parts$part2$coefsnames <- "nu"
  }

  if (model %in% .current.models("iid")) {
    if (out$npar[1] > 1) {
      .stop.with.message(
        paste0(
          " Unknown parameters in part 1",
          "\n For model ", model, "(iid sample) part 1 must include only `alpha`"
        )
      )
    }
  }

  if (endsWith(model, "V")) {
    for (i in 1:nparts) {
      if (out$npar[i] > 0) {
        key <- .parts[i]
        parts[[key]]$coefsnames <- paste0("P", i, ": ", parts[[key]]$coefsnames)
      }
    }
  }

  out$coefs <- NULL
  out$coefsnames <- NULL
  for (i in 1:nparts) {
    key <- .parts[i]
    out$coefs <- c(out$coefs, parts[[key]]$coefs)
    out$coefsnames <- c(out$coefsnames, parts[[key]]$coefsnames)
  }
  if (!is.null(out$coefs)) names(out$coefs) <- out$coefsnames
  cf.names <- .names.coefs(MODEL)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # helper function.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Merges the outputs from part 1 and part 2 for each parameter
  .merge.coefs <- function(part1, part2) {
    out <- c()

    # number of fixed parameters
    out$nfix <- c(part1$nfix, part2$nfix)
    nfix <- c(max(1, out$nfix[1]), max(1, out$nfix[2]))

    # fixed lags
    out$flags <- matrix(0L, nrow = 2, ncol = max(nfix))
    out$flags[1, 1:nfix[1]] <- part1$flags
    out$flags[2, 1:nfix[2]] <- part2$flags

    # fixed values
    out$fvalues <- matrix(0, nrow = 2, ncol = max(nfix))
    out$fvalues[1, 1:nfix[1]] <- part1$fvalues
    out$fvalues[2, 1:nfix[2]] <- part2$fvalues

    return(out)
  }

  out$nfix <- NULL
  for (key in cf.names) {
    out[[key]] <- .merge.coefs(parts$part1[[key]], parts$part2[[key]])
    out$nfix <- cbind(out$nfix, out[[key]]$nfix)
  }
  colnames(out$nfix) <- cf.names

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# For each parameter, convert the bounds to pass to FORTRAN.
# bounds:
#   0 = no bounds
#   1 = lower bound only
#   2 = lower and upper bounds
#   3 = upper bound only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.bounds.each <- function(npar, lower, upper) {
  out <- c()

  if (is.null(lower)) lower <- rep(-Inf, npar)
  if (is.null(upper)) upper <- rep(Inf, npar)

  out$nbd <- rep(0L, npar)
  out$lower <- lower
  out$upper <- upper

  w1 <- (lower > -Inf) & (upper == Inf)
  w2 <- (lower > -Inf) & (upper < Inf)
  w3 <- (lower == -Inf) & (upper < Inf)

  # lower bound only
  if (sum(w1) > 0) out$nbd[w1] <- 1L

  # upper and lower bound
  if (sum(w2) > 0) out$nbd[w2] <- 2L

  # upper bound only
  if (sum(w3) > 0) out$nbd[w3] <- 3L

  invisible(out)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# For each part, convert the bounds to pass to FORTRAN.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @importFrom stats setNames
.convert.bounds.part <- function(model, nfix, lower, upper, order) {
  #  setting the bounds
  lwr <- NULL
  upr <- NULL
  nbd <- NULL

  # d - not implemented for BARC models
  nfix.comp <- c(
    alpha = 1,
    beta = order[["nreg"]],
    phi = order[["p"]],
    theta = order[["q"]],
    d = 1,
    u0 = 1
  )
  keys <- .names.coefs(type = model)
  pnames <- NULL

  for (key in keys) {
    npar <- nfix.comp[key] - nfix[key]
    if (npar > 0) {
      cb <- .convert.bounds.each(
        npar = npar, lower = lower[[key]], upper = upper[[key]]
      )
      name <- key
      if (length(cb$lower) > 1) {
        name <- paste0(key, seq_along(cb$lower))
      }
      lwr <- c(lwr, setNames(cb$lower, name))
      upr <- c(upr, setNames(cb$upper, name))
      nbd <- c(nbd, setNames(cb$nbd, name))
    }
  }

  return(list(lower = lwr, upper = upr, nbd = nbd))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Convert the bounds to pass to FORTRAN.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.bounds <- function(model, nfix, lower, upper, order, map) {
  # in case the lower and upper limits for nu are missing, use
  # default values for nu (merge with user provided values)
  lower <- .check.coefs.defaults(
    model = model, object = lower, type = "lower bounds"
  )
  upper <- .check.coefs.defaults(
    model = model, object = upper, type = "upper bounds"
  )
  if (!endsWith(model, "V")) {
    if (is.null(lower$part2$alpha)) lower$part2$alpha <- 0
    if (is.null(upper$part2$alpha)) upper$part2$alpha <- Inf
  }

  if (model == "BARC") {
    # fix the lower and upper values for u0 (if needed)
    lower$part1$u0 <- max(0, lower$part1$u0)
    upper$part1$u0 <- min(1, upper$part1$u0)
    if (lower$part1$u0 > upper$part1$u0) {
      .stop.with.message("Parameter u0: \n  - lower limit > upper limit")
    }

    # fix the lower and upper values for theta (if needed)
    lu.theta <- NULL
    if (order[1, "q"] > 0) {
      lu.theta <- .check.lu.theta(
        map = map, lower = lower$part1$theta, upper = upper$part1$theta
      )
    }
    lower$part1$theta <- lu.theta$lower
    upper$part1$theta <- lu.theta$upper
  }

  #  setting the bounds
  out <- list()
  for (i in 1:length(.parts)) {
    temp <- .convert.bounds.part(
      model = model,
      nfix = nfix[i, ], order = order[i, ],
      lower = lower[[.parts[i]]], upper = upper[[.parts[i]]]
    )
    out$lower <- c(out$lower, temp$lower)
    out$upper <- c(out$upper, temp$upper)
    out$nbd <- c(out$nbd, temp$nbd)
  }
  out$nbd <- as.integer(out$nbd)
  names(out$nbd) <- names(out$lower)

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function converts the link to the corresponding integer to be
# passed to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.convert.link <- function(link) {
  lk <- .get.link.info(link, only.code = TRUE)$code
  if (any(is.na(lk))) {
    .stop.with.message(
      paste0(
        " Invalid link(s) passed to the function:",
        "\n  - ", paste0("'", link[is.na(lk)], "'", collapse = ", "),
        "\n Valid links are:",
        "\n  - ", paste0(.valid.links("names"), collapse = ", ")
      )
    )
  }
  return(structure(as.integer(lk), names = names(link)))
}
