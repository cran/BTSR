# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal constants.
#  - The valid error scales currently implemented
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0: added this constants
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# list of all models
.btsr.models <- structure(
  rbind(
    # Beta
    c("Beta", "BETA", "BREG", "BARMA", "BARFIMA", "BARC"),
    c("", "", "BREGV", "BARMAV", "BARFIMAV", ""),
    # Gamma
    c("Gamma", "GAMMA", "GREG", "GARMA", "GARFIMA", ""),
    c("", "", "GREGV", "GARMAV", "GARFIMAV", ""),
    # Kuma
    c("Kumaraswamy", "KUMA", "KREG", "KARMA", "KARFIMA", ""),
    c("", "", "KREGV", "KARMAV", "KARFIMAV", ""),
    # Matsu
    c("Matsuoka", "MATSU", "MREG", "MARMA", "MARFIMA", ""),
    # UL
    c("Unit-Lindley", "UL", "ULREG", "ULARMA", "ULARFIMA", ""),
    # UW
    c("Unit-Weibull", "UW", "UWREG", "UWARMA", "UWARFIMA", ""),
    c("", "", "UWREGV", "UWARMAV", "UWARFIMAV", "")
  ),
  dimnames = list(
    NULL,
    c(
      "Distribution", "i.i.d.", "Regression",
      "Short-Memory", "Long-Memory", "Chaotic"
    )
  )
)

# default map for BARC models
.default.map.barc <- 4
.default.u0.barc <- pi / 4

# error scale
.current.error.scale <- c(data = 0, predictive = 1)
.scale.message <- paste0(
  " Invalid error.scale. Valid options are",
  "\n  - 0: data scale",
  "\n  - 1: predictive scale"
)

# link names
.current.link.names <- c("g1", "g11", "g12", "g13", "g2", "g21", "g22", "g23")

# extra configurations
.default.extra.configs <- list(
  m = 0, llk = TRUE, sco = FALSE, info = FALSE, extra = FALSE
)

# parts in the model
.parts <- c("part1", "part2")

# default names for prediction
.predict.names <- c(
  "gy", "fitted.values", "etat", "g22varthetat", "error", "forecast"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Maps for BARC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  1 = (kx)(mod 1).  k integer
#  2 = Rafael's map.  0 <= theta <= 1
#  3 = logistic map. 0 <= theta  <= 4
#  4 = Manneville-Pomeau. 0 < theta < 1
#  5 = Lasota-Mackey's map. No theta
.barc.maps <- data.frame(
  map = c(1:5),
  r = c(1, 1, 1, 1, 0),
  lower = c(1, 0, 0, 0, NA),
  upper = c(Inf, 1, 4, 1, NA)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Prints the current implemented models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @importFrom graphics lines mtext par rect text
.BTSR.model.table <- function(do.plot = interactive()) {
  # Create text table
  header <- colnames(.btsr.models)
  table_data <- .btsr.models
  nc <- ncol(table_data)
  nr <- nrow(table_data)

  if (do.plot) {
    # Set up empty plot
    par(mar = c(0.5, 0.5, 2, 0.5))
    plot(0, 0,
      type = "n", xlim = c(0, nc), ylim = c(0, nr + 1),
      axes = FALSE, xlab = "", ylab = "",
      main = "The BTSR package implements the following models:"
    )
    # Subtitle
    mtext("For more details, please refer to the documentation",
      side = 1, line = -1, cex = 1, font = 3, adj = 0
    )

    # Column positions
    col_x <- 1:nc - 0.5
    col_width <- 0.9
    row_height <- 0.8

    # Header with proper formatting
    header_bg <- c(0, nr + 0.5, nc, nr + 1)
    rect(header_bg[1], header_bg[2], header_bg[3], header_bg[4],
      col = "gray70", border = NA
    )
    text(col_x, rep(nr + 0.75, nc), header, font = 2, cex = 0.9)

    n <- nrow(table_data)
    opt <- c("gray88", "white")
    color_list <- opt[1]
    for (i in 2:nr) {
      if (table_data[i, 1] != "") opt <- rev(opt)
      color_list <- c(color_list, opt[1])
    }

    # Draw table cells with proper alignment
    for (i in 1:n) {
      y_pos <- n + 1 - i
      fill_col <- color_list[i]

      # Cell background
      rect(0, y_pos - row_height, nc, y_pos + row_height / 2,
        col = fill_col, border = NA
      )

      # Cell text centered
      for (j in 1:nc) {
        text(col_x[j], y_pos, table_data[i, j], cex = 0.8)
      }
    }

    # Grid lines (thinner and lighter)
    for (x in 0:nc) {
      lines(rep(x, 2), c(n + 0.5, 0.5), col = "gray20", lwd = 0.5)
    }
    lines(c(0, nc), rep(0.5, 2), col = "gray20", lwd = 0.5)
  } else {
    # Text output for terminals/non-interactive use
    cat("\nThe BTSR package implements the following models:\n")
    cat("For more details, please refer to the documentation\n\n")

    # Calculate column widths
    all_rows <- rbind(header, table_data)
    max_widths <- apply(all_rows, 2, function(x) max(nchar(x)))

    # Print header
    cat(paste(
      sapply(1:nc, function(i) {
        sprintf("%-*s", max_widths[i], header[i])
      }),
      collapse = " | "
    ), "\n")

    # Print divider
    cat(paste(rep("-", sum(max_widths) + 12), collapse = ""), "\n")

    # Print table rows
    for (i in 1:nr) {
      cat(paste(
        sapply(1:nc, function(j) {
          sprintf("%-*s", max_widths[j], table_data[i, j])
        }),
        collapse = " | "
      ), "\n")
    }
    cat("\n")
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Set the expected arguments that each function must return
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.output.names <- function(model) {
  is.v <- endsWith(model, "V")

  # fitted.values
  fitted <- "mut"
  eta <- "eta1t"
  xreg <- "xreg1"
  forecast <- c("mut", "eta1t")
  xnew <- c("xnew1")
  error <- c("error1")
  if (is.v) {
    fitted <- c(fitted, "nut", "varthetat")
    eta <- c(eta, "eta2t")
    xreg <- c(xreg, "xreg2")
    forecast <- c(forecast, "nut", "varthetat", "eta2t")
    xnew <- c(xnew, "xnew2")
    error <- c(error, "error2")
  }
  if (model == "BARC") {
    fitted <- c("Tt", fitted)
    forecast <- c("Tt", forecast)
  }

  return(list(
    fitted = fitted, eta = eta, xreg = xreg, error = error,
    forecast = forecast, xnew = xnew
  ))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Set the expected arguments for each function: sim, extract and fit
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.arguments <- function(type) {
  fit <- type == "fit"
  sim <- type == "sim"

  coefs.configs <- c("xregar") # all
  if (!fit) {
    # sim and extract
    coefs.configs <- c(coefs.configs, "coefs")
  }
  if (!sim) {
    # fit and extract
    coefs.configs <- c(coefs.configs, "lags", "fixed.values", "fixed.lags")
    if (fit) {
      # only fit
      coefs.configs <- c(coefs.configs, "start", "ignore.start", "lower", "upper")
    }
  }

  series <- "xreg"
  if (sim) {
    series <- c(series, "n", "burn") # sim
  } else {
    series <- c(series, "yt", "xnew", "nnew") # fit and extract
  }

  order.configs <- c("p", "q", "inf")
  if (fit) order.configs <- c(order.configs, "d")

  extra.configs <- NULL
  if (!sim) extra.configs <- c("m", "llk", "sco", "info", "extra")


  # returning the list of arguments to be passed to extract and fit functions
  list(
    pdist = c("rho", "y.lower", "y.upper"),
    series = series,
    order.configs = order.configs,
    data.start = c("y.start", "xreg.start", "vt.start", "e2.start"),
    link.configs = c(
      "linkg", "linkh", "error.scale", "configs.linkg", "configs.linkh"
    ),
    coefs.configs = coefs.configs,
    extra.configs = extra.configs
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function returns the FORTRAN code for the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.name.to.code <- function(model) {
  # Here model is the base model
  switch(model,
    BARC = 0L,
    BARFIMA = 1L,
    GARFIMA = 2L,
    KARFIMA = 3L,
    MARFIMA = 4L,
    ULARFIMA = 5L,
    UWARFIMA = 6L
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# List of coefficients names
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - Added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.names.coefs <- function(type) {
  # get the type of model: iid, REG, ARMA
  if (type %in% c("BARC", "iid", "reg", "arma", "arfima", "all")) {
    model <- type
  } else {
    model <- tolower(.get.model.type(type))
    if (model == "arma") {
      temp <- grep("ima", type, ignore.case = TRUE)
      if (length(temp) > 0) model <- "arfima"
    }
  }

  # return the expected names
  switch(model,
    BARC = c("alpha", "beta", "phi", "theta", "u0"),
    iid = c("alpha"),
    reg = c("alpha", "beta"),
    arma = c("alpha", "beta", "phi", "theta"),
    arfima = c("alpha", "beta", "phi", "theta", "d"),
    all = c("alpha", "beta", "phi", "theta", "d", "u0")
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Returns information on the coefficients for updating purposes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - Added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.need.update <- function(model) {
  if (model != "BARC") model <- "arfima"
  need <- .names.coefs(type = model)
  scalar <- need %in% c("alpha", "d", "u0")
  key <- rep("dummy", length(need))
  key[need == "beta"] <- "nreg"
  key[need == "phi"] <- "p"
  key[need == "theta"] <- "q"
  return(list(need.update = need, scalar = scalar, npar.key = key))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function returns the current implemented models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.current.models <- function(type = "iid") {
  w <- .btsr.models[, "Distribution"] != ""
  models <- .btsr.models[w, "Long-Memory"]

  switch(
    EXPR = type,
    base = models,
    iid = .btsr.models[w, "i.i.d."], # iid samples
    reg = .btsr.models[w, "Regression"], # regression models
    arma = gsub("(FIMA)$", "", models), # short/long memory models
    nonu = c("MARFIMA", "ULARFIMA") # one parameter
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function returns the type of bounds for the implemented models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.current.models.bounds <- function(model) {
  bounds <- "both"
  gam <- unname(.btsr.models[.btsr.models[, "Distribution"] == "Gamma", -1])
  gam <- c(gam, paste0(gam, "V"))
  if (model %in% gam) {
    bounds <- "only.lower"
  }
  return(bounds)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Returns information on the valid links
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.valid.links <- function(type = "all") {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  0 = polynomial/linear: a real, b real
  #      g(x) = ax^b,
  #      from (l,u) to (min(g(l),g(u)), max(g(l),g(u)))
  #
  #  1 = logit:
  #      g(x) = log((x-l)/(u-x))
  #      from (l,u) to (-Inf, Inf)
  #
  #  2 = log:
  #      g(x) = log(x-l)
  #      from (l, u) to (-Inf, log(u-l))
  #
  #  3 = loglog:
  #      g(x) = log(-log((x-l)/(u-l)))
  #      from (l,u) to (-Inf, Inf)
  #
  #  4 = cloglog:
  #      g(x) = log(-log(1-(x-l)/(u-l)))
  #      from (l,u) to (-Inf, Inf)
  #
  #  5 = SIP = shifted inverse power
  #      g(x) = 1/(a+x)^b, a = 0,1 and b real
  #      from (l,u) to (g(u), g(l))
  #
  #      default link for beta regression:
  #      g(x) = 1/(1+x)
  #
  #      default link for gamma regression:
  #      g(x) = 1/x
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  names <- c(
    "linear", "polynomial",
    "logit", "log", "loglog", "cloglog",
    "SIP", "SIP1", "SIP0"
  )
  if (type == "names") {
    return(names)
  }
  # corresponding FORTRAN codes
  codes <- c(
    linear = 0, polynomial = 0,
    logit = 1, log = 2, loglog = 3, cloglog = 4,
    SIP = 5, SIP1 = 5, SIP0 = 5
  )
  if (type == "code") {
    return(data.frame(link = names, code = codes[names]))
  }
  if (type == "code.base") {
    basenames <- names[-c(1, 8, 9)]
    return(data.frame(link = basenames, code = codes[basenames]))
  }

  # can be used as gi1 and gi2? (iid models can use any link)
  allowed.unbounded.gi1 <- c(
    linear = TRUE, polynomial = TRUE,
    logit = FALSE, log = FALSE, loglog = FALSE, cloglog = FALSE,
    SIP = TRUE, SIP1 = TRUE, SIP0 = TRUE
  )
  allowed.bounded.gi1 <- c(
    linear = FALSE, polynomial = FALSE,
    logit = TRUE, log = FALSE, loglog = TRUE, cloglog = TRUE,
    SIP = FALSE, SIP1 = FALSE, SIP0 = FALSE
  )
  allowed.only.lower.gi1 <- c(
    linear = FALSE, polynomial = FALSE,
    logit = FALSE, log = TRUE, loglog = FALSE, cloglog = FALSE,
    SIP = FALSE, SIP1 = FALSE, SIP0 = FALSE
  )
  allowed.only.lower.gi2 <- c(
    linear = TRUE, polynomial = TRUE,
    logit = FALSE, log = TRUE, loglog = FALSE, cloglog = FALSE,
    SIP = TRUE, SIP1 = TRUE, SIP0 = TRUE
  )

  # can be used as g2(nu), with 0 < nu < Inf;
  # will g2(nu) have the same type of bounds as nu?
  allowed.only.lower.g2 <- c(
    linear = TRUE, polynomial = TRUE,
    logit = FALSE, log = TRUE, loglog = FALSE, cloglog = FALSE,
    SIP = TRUE, SIP1 = TRUE, SIP0 = TRUE
  )
  preserve.bounds <- c(
    linear = TRUE, polynomial = TRUE,
    logit = FALSE, log = FALSE, loglog = FALSE, cloglog = FALSE,
    SIP = TRUE, SIP1 = FALSE, SIP0 = TRUE
  )
  from.lower.to.bounds <- c(
    linear = FALSE, polynomial = FALSE,
    logit = FALSE, log = FALSE, loglog = FALSE, cloglog = FALSE,
    SIP = FALSE, SIP1 = TRUE, SIP0 = FALSE
  )

  return(data.frame(
    link = names,
    codes = codes[names],
    allowed.unbounded.gi1 = allowed.unbounded.gi1[names],
    allowed.bounded.gi1 = allowed.bounded.gi1[names],
    allowed.only.lower.gi1 = allowed.only.lower.gi1[names],
    allowed.only.lower.gi2 = allowed.only.lower.gi2[names],
    allowed.only.lower.g2 = allowed.only.lower.g2[names],
    preserve.bounds = preserve.bounds[names],
    from.lower.to.bounds = from.lower.to.bounds[names]
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function returns default lower and upper bounds for the
# support of y and the name of required parameters for some models
#  - KARFIMA models need update before sending to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.dist.defaults <- function(model) {
  # required parameters that will have default values in case
  # the user does not provide
  p.names <- switch(model,
    KARFIMA  = c("rho", "y.lower", "y.upper"),
    UWARFIMA = c("rho")
  )
  # default range and default parameter values
  pdist <- switch(model,
    BARC = c(y.lower = 0, y.upper = 1),
    BARFIMA = c(y.lower = 0, y.upper = 1),
    KARFIMA = c(rho = 0.5, y.lower = 0, y.upper = 1),
    GARFIMA = c(y.lower = 0, y.upper = Inf),
    MARFIMA = c(y.lower = 0, y.upper = 1),
    ULARFIMA = c(y.lower = 0, y.upper = 1),
    UWARFIMA = c(rho = 0.5, y.lower = 0, y.upper = 1)
  )

  return(list(p.names = p.names, pdist = pdist))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Returns the default links, limits and default constants
#  - g1 is a dummy argument to pass to FORTRAN
#  - for BARC models, limits for g2, g21 and g22 are ignored.
#  - limits for g2(nu) will be updated in FORTRAN
#  - for KARFIMA, limits for g11 and g12 need update before sending to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1.0.0:
#  - added this function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.links.default <- function(model, is.iid = TRUE, endsV = FALSE) {
  # default error scale for part 1 (predictive scale)
  error.scale <- 1

  # default links for all models
  linkg <- list(
    g1 = "polynomial", g11 = "polynomial", g12 = "polynomial", g13 = error.scale,
    g2 = "polynomial", g21 = "polynomial", g22 = "polynomial", g23 = "polynomial"
  )

  # default configs for all models
  lconfig <- rbind(
    g1 = c(0, 1, 1, 1), # dummy argument to pass to FORTRAN
    g11 = c(0, 1, 1, 1), # mu, y
    g12 = c(0, 1, 1, 1), # y
    g13 = c(-Inf, Inf, 1, 1), # error in part 1 (will be updated in FORTRAN)
    g2 = c(0, Inf, 1, 1), # nu
    g21 = c(0, Inf, 1, 1), # vt = g2(nu) (will be updated in FORTRAN)
    g22 = c(0, Inf, 1, 1), # vt = g2(nu) (will be updated in FORTRAN)
    g23 = c(0, Inf, 1, 2) # error in part 2 (will be updated in FORTAN)
  )
  colnames(lconfig) <- c("lower", "upper", "ctt", "power")

  # for iid models keep the default
  if (model != "BARC" && !is.iid) {
    # update g11 and g12
    linkg[c("g11", "g12")] <-
      switch(model,
        GARFIMA = "log",
        MARFIMA = "cloglog",
        "logit"
      )
  }

  # update g2, g12 and g22 only for time varying nu
  if (model != "BARC" && endsV) {
    linkg$g2 <- "default"
    linkg[c("g21", "g22")] <-
      switch(model,
        BARFIMA = "logit",
        "log"
      )
  }

  if (linkg$g2 == "default" && model == "GARFIMA") lconfig["g2", "ctt"] <- 0
  if (model == "GARFIMA") lconfig[c("g11", "g12"), "upper"] <- Inf # mu and y
  if (model == "BARC") {
    linkg$h <- "polynomial"
    lconfig <- rbind(lconfig, h = c(0, 1, 1, 1)) # h(T(u0))
  }

  return(
    list(linkg = linkg, error.scale = error.scale, lconfig = lconfig)
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function: Used to check if theta is compatible with the map
# selected by the user
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  1 = (kx)(mod 1).  k integer
#  2 = Rafael's map.  0 <= theta <= 1
#  3 = logistic map. 0 <= theta  <= 4
#  4 = Manneville-Pomeau. 0 < theta < 1
#  5 = Lasota-Mackey's map. No theta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.barc.starting.theta <- function(map) {
  theta <- c(3, 0.5, 3.5, 0.5, 0)
  theta[map]
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function provides starting values based on the link function used
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.link.starting.values <- function(link, l = 0, u = 1) {
  out <- switch(link,
    polynomial = 0,
    logit = (u - l) / 2,
    log = l + 1,
    loglog = l + (u - l) * exp(-1),
    cloglog = l + (u - l) * (1 - exp(-1)),
    0
  )
  return(out)
}
