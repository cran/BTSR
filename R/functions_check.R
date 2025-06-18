# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function validate parameters and assign defaults
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.dist.parameters <- function(model, p.names, pdist, param.list) {
  # if there are no parameters to update, return the default values
  if (is.null(p.names)) {
    return(pdist)
  }

  # check if any parameter value was provided by the user and updates
  # the vector of defaults
  nm <- names(pdist)
  msg <- NULL
  for (param in p.names) {
    if (is.null(param.list[[param]])) { # Assign default
      msg <- paste0(
        msg, paste0(
          if (!is.null(msg)) "\n ",
          if (is.null(msg)) " ",
          param, " is missing. Assuming ", param, " = ", pdist[param]
        )
      )
    } else {
      pdist[param] <- param.list[[param]]
    }
  }
  if (!is.null(msg)) .warning.with.message(msg)
  return(pdist)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Checks if the data has any NA and if it is the the correct range.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.data <- function(yt, lower, upper, openIC) {
  out <- list()

  # checking for NA's
  if (sum(is.na(yt)) > 0) {
    out$conv <- 123
    out$message <- "NA in the data"
    .stop.with.message(" NA's in the data are not allowed")
  }

  # checking if y is in the correct range
  if (openIC[1]) {
    a <- !(min(yt) > lower)
  } else {
    a <- !(min(yt) >= lower)
  }
  if (openIC[2]) {
    b <- !(max(yt) < upper)
  } else {
    b <- !(max(yt) <= upper)
  }
  if (a || b) {
    out$conv <- 123
    out$conv.message <- "Out of range"
    .stop.with.message(
      paste0(" OUT OF RANGE. yt must be bewteen ", lower, " and ", upper)
    )
  }

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Check if the coefficient is missing and set the value to zero
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @importFrom stats setNames
.check.xreg.format <- function(xreg) {
  # set default value
  if (is.null(xreg)) {
    return(setNames(vector("list", length(.parts)), .parts))
  }

  # for compatibility with the old structure:
  #   - checking the structure and converting to the new format
  if (!("list" %in% class(xreg))) {
    return(list(part1 = xreg, part2 = NULL))
  }

  # checking compatibility with the new structure
  if (!any(.parts %in% names(xreg))) {
    .stop.with.message(" Unknow structure for xreg")
  }

  return(xreg)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Check the structure and provides default initial values for xreg.start
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.xreg.start <- function(xreg.start, p, xregar, xreg, nreg) {
  # Check if the required structure is provided
  xreg.start <- .check.xreg.format(xreg.start)

  # checking compatibility of starting values (if provided)
  # Checking starting values
  #   - If the model dos not include xreg in the AR regression,
  #     there is no use for initial values of xreg.
  for (i in 1:length(.parts)) {
    key <- .parts[i]
    if (!is.null(xreg.start[[key]]) && xregar[i]) {
      if (length(xreg.start[[key]]) != nreg[i]) {
        .stop.with.message(
          paste0(" Part ", i, ":, \n length(xreg.start) is not the same as nreg")
        )
      }
    }

    if (!xregar[i]) { # set dummy values
      xreg.start[[key]] <- rep(0, max(1, nreg[i]))
    } else {
      if (is.null(xreg.start[[key]])) { # set the default values
        xreg[[key]] <- as.matrix(xreg[[key]])
        xreg.start[[key]] <- apply(xreg[[key]], 2, function(x) mean(x[1:p[i]]))
      }
    }
  }
  return(xreg.start)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function checks the compatibility of beta and xreg
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.size.beta <- function(beta, nreg) {
  for (k in 1:length(nreg)) {
    key <- paste0("part", k)
    if (nreg[k] > 0) {
      if (is.null(beta[[key]])) {
        .stop.with.message(
          paste0(
            " Part ", k, ":",
            "\n xreg is provided and beta is missing with no default"
          )
        )
      } else {
        if (length(beta[[key]]) != nreg[k]) {
          .stop.with.message(
            paste0(
              " Part ", k, ":",
              "\n beta and nreg are not compatible",
              "\n length of beta = ", length(beta[[key]]),
              "\n number of regressors = ", nreg[k]
            )
          )
        }
      }
    }
  }
  return(TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function
#   - checks and converts the coefficients to the new format
#   - sets default values based on the model
# type indicates if the object is a list of coefficients, lower/upper bounds,
# lags, fixed.values or fixed.lags
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.coefs.defaults <- function(model, object, type = "generic",
                                  is.coef = FALSE, is.sim = FALSE) {
  # check if any value is provided
  if (is.null(object) || length(object) == 0) {
    # For simulation purposes
    #  - the coefs are required. Other objects are ignored
    # For extracting or fitting
    #  - all objects can be NULL
    if (is.coef && is.sim) {
      .stop.with.message(
        paste0(" `coefs` (list of coefficients) is missing with no default")
      )
    } else {
      return(list(part1 = list(), part2 = list()))
    }
  }

  # if the object is provided, check the structure
  #  - the object must be a list
  if (!"list" %in% class(object)) {
    .stop.with.message(
      paste0(" ", type, " must be a named list")
    )
  }

  # check the object names and the expected names
  ob.names <- names(object)
  cf.names <- .names.coefs(type = model)

  # list of models related to one parametric distributions
  has.nu <- !(.get.base.model(model) %in% .current.models("nonu"))

  # check the structure
  new <- any(.parts %in% ob.names)
  old <- any(c(cf.names, "nu") %in% ob.names)
  if (!old && !new) {
    .stop.with.message(
      paste0(" coefficients structure unknown in object", type)
    )
  }

  # in case the old structure is provided, convert it to the new one
  # if the new structure is provides, update it keeping only the expected names
  if (old) {
    #  - for simulation, if the model has the coefficient nu,
    #    then nu must be provided
    if (is.sim && has.nu && !("nu" %in% ob.names)) {
      .stop.with.message(" 'nu' is missing with no default")
    }
    # create the object in the new format
    object <- list(
      part1 = .build.list(mc = object, args = cf.names),
      part2 = list(alpha = object$nu)
    )
  } else {
    # return only the expected names
    for (part in .parts) {
      if (!is.null(object[[part]])) {
        object[[part]] <- .build.list(mc = object[[part]], args = cf.names)
      }
    }
  }
  if (is.sim && model == "BARC" && is.null(object$part1$u0)) {
    .stop.with.message(" 'u0' is missing with no default")
  }
  if (!has.nu) object$part2 <- NULL

  # for simulation, coefs cannot be an empty list
  if (is.sim && is.coef) {
    #  - part 2 only needs to be checked if the model has the parameter nu
    #  - alpha only need to be checked if the model is iid
    check.all <- c(TRUE, has.nu)
    par <- c("mu", "nu")
    for (k in 1:length(.parts)) {
      part <- .parts[k]
      if (check.all[k]) {
        if (is.null(unlist(object[[part]]))) {
          # coefs cannot be null
          .stop.with.message(
            paste0(
              " coefficients of part ", k, " are missing with no default",
              "\n Valid arguments are:",
              "\n  - ", paste0(cf.names, collapse = ", ")
            )
          )
        }
        if (model %in% .current.models("iid")) {
          # alpha must be positive
          if (object[[part]]["alpha"] <= 0) {
            .stop.with.message(
              paste0(" The coefficients provided imply that ", par[k], " < 0")
            )
          }
        }
      }
    }
  }

  return(object)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Not supposed to be called by the user.
# Perform a sequence of testes and set some default values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.coefs <- function(model, coefs.configs, order, type, ignore.coef = FALSE) {
  is.fit <- type == "fit"
  is.sim <- type == "sim"

  # check the coefficients format and set default values (as NULL)
  # - for simulation and extraction at least one must be provided
  # - to fit a model, both can be null since the check is done before
  #   computing starting values
  if (!is.fit &&
    is.null(coefs.configs$coefs) &&
    is.null(coefs.configs$fixed.values)) {
    msg1 <- ifelse(is.sim,
      " Please provide the following:",
      " Please provide at least one of the following:"
    )
    .stop.with.message(
      paste0(
        msg1,
        "\n  - `coefs`: a named list containing the coefficients",
        if (!is.sim) "\n  - `fixed.values`: a named list containing the
     fixed values for the coefficients"
      )
    )
  }

  # for compatibility with the old structure:
  #   - checking the structure and converting to the new format
  #   - set default values for specific models (e.g. iid)
  #   - only checks if nu is provided for simulation subroutines
  #     (assuming the user knows how to use the extraction function)

  if (!ignore.coef) {
    # when fitting a model, skip this step if a list of starting values is
    # not provided by the user.
    coefs <- .check.coefs.defaults(
      model = model, object = coefs.configs$coefs, type = "coefs",
      is.coef = TRUE, is.sim = is.sim
    )
  } else {
    coefs <- NULL
  }

  if (is.sim) {
    # when simulating a model, lags, fixed.values and fixed.lags are not required
    lags <- NULL
    fixed.values <- NULL
    fixed.lags <- NULL
  } else {
    # when extracting or fitting a model, lags, fixed.values and fixed.lags
    # must be checked
    lags <- .check.coefs.defaults(
      model = model, object = coefs.configs$lags, type = "lags"
    )
    fixed.values <- .check.coefs.defaults(
      model = model, object = coefs.configs$fixed.values, type = "fixed.values"
    )
    fixed.lags <- .check.coefs.defaults(
      model = model, object = coefs.configs$fixed.lags, type = "fixed.values"
    )
  }
  if (type == "extract" && model == "BARC") {
    if (is.null(coefs$part1$u0) && is.null(fixed.values$part1$u0)) {
      .stop.with.message(" 'u0' is missing with no default")
    }
  }

  # if not ARMA-type, set p and q to zero (ignore the user input)
  if (.get.model.type(model) != "ARMA") {
    order$p <- c(0, 0)
    order$q <- c(0, 0)
  }

  if (is.fit) {
    # if starting values are not provided, p or q must be provided
    if (ignore.coef) {
      if (is.null(order$p)) {
        .stop.with.message(
          " Please provide 'p' or a list of starting values for phi"
        )
      }

      if (is.null(order$q)) {
        .stop.with.message(
          " Please provide 'q' or a list of starting values for theta"
        )
      }
    }

    temp <- grep("ima", model, ignore.case = TRUE)
    if (length(temp) == 0) {
      # if not *ARFIMA-type, ignore input and set d to FALSE
      order$d <- c(FALSE, FALSE)
    } else {
      # check if d must be included in the model.
      #  - if not: set d as fixed value and equal to zero
      #            remove d from coefs list
      #  - if yes: d cannot be fixed and estimated at the same time.
      part <- .parts
      for (i in 1:length(part)) {
        key <- part[i]
        if (order$d[i] == FALSE) {
          fixed.values[[key]]$d <- 0
          if (!ignore.coef) coefs[[key]]$d <- NULL
        } else {
          check <- !is.null(fixed.values[[key]]$d)
          check <- check && !is.null(coefs[[key]]$d)
          if (check) {
            .stop.with.message(
              paste0(
                " Part ", i, " - ",
                "\n An initial value for 'd' was provided: ", coefs[[key]]$d,
                "\n but 'd' was also fixed as ", fixed.values[[key]]$d, ".",
                "\n  - If you wish to fix d, remove d from starting values or set d = FALSE. ",
                "\n  - If you wish to fit d, remove d from the list of fixed values"
              )
            )
          }
        }
      }
    }
  }

  # set default values for order and update based on the input
  or <- cbind(
    nreg = c(0L, 0L),
    p = c(0L, 0L),
    q = c(0L, 0L),
    inf = c(0L, 0L)
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # helper function: update the order of the model based on the input
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  .update.order <- function(kvalue, parname, cf, fl, fv) {
    # current parts in the model
    keys <- .parts

    # update the order "k"
    if (is.null(kvalue)) {
      is.null <- is.null(fl) && is.null(fv)
      # loop in parts
      kvalue <- NULL
      for (i in 1:length(keys)) {
        key <- keys[i]
        if (is.null) {
          ptemp <- 0
        } else {
          ptemp <- max(length(fv[[key]][[parname]]), length(fl[[key]][[parname]]))
        }
        kvalue <- c(kvalue, length(cf[[key]][[parname]]) + ptemp)
      }
    }
    if (length(kvalue) == 1) kvalue <- c(kvalue, kvalue)
    return(as.integer(kvalue))
  }

  # update inf and nreg
  keys <- c("nreg", "inf")
  vals <- c(0, 1000) # default values (inf might be updated in FORTRAN)
  for (i in 1:length(keys)) {
    key <- keys[i]
    if (is.null(order[[key]])) order[[key]] <- vals[i]
    or[, key] <- .update.order(
      kvalue = order[[key]], parname = key, cf = NULL, fl = NULL, fv = NULL
    )
  }

  # update p and q
  knames <- c("p", "q")
  parname <- c("phi", "theta")
  for (i in 1:length(knames)) {
    key <- knames[i]
    or[, key] <- .update.order(
      kvalue = order[[key]], parname = parname[i],
      cf = coefs, fl = fixed.lags, fv = fixed.values
    )
  }

  # if the model is not time varying the orders in part 2 must be zero.
  if (!endsWith(model, "V")) or[2, ] <- 0L
  if (model == "BARC") or[, "inf"] <- 0L

  # convert xregar
  xregar <- coefs.configs$xregar
  if (is.null(xregar)) xregar <- TRUE
  if (length(xregar) == 1) {
    xregar <- rep(xregar, length.out = 2)
  }
  xregar[or[, "nreg"] == 0 | or[, "p"] == 0] <- FALSE

  return(list(
    order = or, xregar = xregar, coefs = coefs, lags = lags,
    fixed.values = fixed.values, fixed.lags = fixed.lags, fitd = order$d
  ))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function checks if the link object has the correct format and
# all the elements needed:
#  - version 0.0.x: link is a vector
#  - version 1.0.x: link must be a list and have the elements
#                   g11, g12, g2, g21 and g22
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.link.format <- function(link, is.config = FALSE) {
  if (is.null(link)) {
    return(list())
  }

  #  For compatibility with the previous version of the package:
  #   - if link is a character or a two character vector
  #   - if configs is a vector of size 1 or 2.
  if (!("list" %in% class(link))) {
    if (length(link) == 1) link <- list(link, link)
    if (length(link) == 2) {
      link <- as.list(link)
      names(link) <- c("g11", "g12")
      other <- .current.link.names[!.current.link.names %in% c("g11", "g12")]
      if (is.config) {
        link[other] <- 1
      } else {
        link[other] <- "linear"
      }
      return(link)
    }
  }

  if (is.config) {
    return(link)
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # check the arguments
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (!all(names(link) %in% .current.link.names)) {
    nm <- names(link)
    nmnot <- nm[!nm %in% .current.link.names]
    .stop.with.message(
      paste0(
        " Unused arguments in linkg: \n ",
        paste0("'", nmnot, "'", collapse = ", ")
      )
    )
  }

  return(link)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Checks if the linkg object has all elements
# and sets some defaults based on the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.link.elements <- function(model, link, configs, lk.default, cf.default, escale) {
  # update escale.
  if (!is.null(escale)) lk.default$g13 <- escale

  # for iid samples all links are the identity
  if (model %in% .current.models("iid")) {
    return(list(link = lk.default, configs = cf.default))
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # helper functions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # copies the information of non-missing link to missing link
  .copy.link <- function(lk, from, to, cf) {
    # copy link, power, and the constant
    lk[[to]] <- lk[[from]]
    cf$power[[to]] <- cf$power[[from]]
    cf$ctt[[to]] <- cf$ctt[[from]]
    return(list(link = lk, configs = cf))
  }
  # check if any link is missing,
  # If one link is missing, makes a copy of the other one.
  # if both are present, just return the object
  .copy.check <- function(g1, g2, link, cf.link) {
    # checks if any link is missing
    l1 <- is.null(link[[g1]])
    l2 <- is.null(link[[g2]])
    if (l1) {
      temp <- .copy.link(lk = link, from = g2, to = g1, cf = cf.link)
    } else if (l2) {
      temp <- .copy.link(lk = link, from = g1, to = g2, cf = cf.link)
    } else {
      temp <- list(link = link, configs = cf.link)
    }
    return(list(link = temp$link, configs = temp$configs))
  }
  # updates the values in the default configurations list using the values
  # in the input link. In case any value is missing, keeps the default values
  .update.cf.default <- function(cf.df, cf, keys) {
    for (key in keys) {
      if (!is.null(cf$power[[key]])) cf.df$power[[key]] <- cf$power[[key]]
      if (!is.null(cf$ctt[[key]])) cf.df$ctt[[key]] <- cf$ctt[[key]]
    }
    return(cf.df)
  }

  # g1 is dummy
  # check if g11 and g12 are provided
  # if they are missing, use the default values based on the model
  g11 <- !is.null(link$g11)
  g12 <- !is.null(link$g12)
  if (g11 || g12) {
    # if one is missing, make a copy
    temp <- .copy.check(g1 = "g11", g2 = "g12", link = link, cf.link = configs)
    lk.default[c("g11", "g12")] <- temp$link[c("g11", "g12")]
    cf.default <- .update.cf.default(cf.default, temp$configs, keys = c("g11", "g12"))
  }

  # only BARC has link h
  if (model == "BARC") {
    if (!is.null(link[["h"]])) lk.default[["h"]] <- link[["h"]]
    cf.default <- .update.cf.default(cf.default, configs, keys = "h")
  }

  # if nu is fixed then g2, g21, g22 and g23 are the identity function
  if (any(!endsWith(model, "V"))) {
    # model can be a vector with aliases
    return(list(link = lk.default, configs = cf.default))
  }

  # check if g2 is provided
  if (!is.null(link$g2)) lk.default$g2 <- link$g2
  cf.default <- .update.cf.default(cf.default, configs, keys = "g2")

  # check if g21 and g22 are provided
  # if they are missing, use the default values based on the model
  g21 <- !is.null(link$g21)
  g22 <- !is.null(link$g22)
  if (g21 || g22) {
    temp <- .copy.check(g1 = "g21", g2 = "g22", link = link, cf.link = configs)
    lk.default[c("g21", "g22")] <- temp$link[c("g21", "g22")]
    cf.default <- .update.cf.default(cf.default, temp$configs, keys = c("g21", "g22"))
  }

  # check if g23 is provided
  if (!is.null(link$g23)) lk.default$g23 <- link$g23
  cf.default <- .update.cf.default(cf.default, configs, keys = c("g23"))

  return(list(link = lk.default, configs = cf.default))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Checks if the link can be evaluated
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.link.eval <- function(model, link, allow.any) {
  # Here model is the particular model

  # check the type of bounds for the current model
  ybounds <- .current.models.bounds(model)

  # Check if g11 matches the model and g12 can be evaluated
  # bounds = "only.lower"  (GAMMA-based models)
  # - g11 and g12 cannot depend on upper bounds,  0 < y, mu < Inf
  #
  # bounds = "both"
  # - g11inv must be bounded,  a < y, mu < b
  #   However, for iid samples we can use unbounded functions and
  #   control the lower and upper bounds for the parameters.
  lk <- .get.link.info(c(link$g11, link$g12))
  if (!allow.any[1]) {
    if (ybounds == "only.lower") {
      g11test <- !lk[lk$link == link$g11, "allowed.only.lower.gi1"]
      g12test <- !lk[lk$link == link$g12, "allowed.only.lower.gi2"]
      if (g11test || g12test) {
        .stop.with.message(
          paste0(
            " The selected links are not compatible.",
            "\n g11 = ", link$g11, ", g12 = ", link$g12,
            "\n For the selected model, 0 < mu, y < Inf, so that",
            if (g11test) "\n  - g11 cannot depend on upper bounds",
            if (g12test) "\n  - g12 cannot depend on upper bounds",
            "\n Please select another link"
          )
        )
      }
    } else {
      g11test <- !lk[lk$link == link$g11, "allowed.bounded.gi1"]
      if (g11test) {
        .stop.with.message(
          paste0(
            " The selected links are not compatible.",
            "\n g11 = ", link$g11,
            "\n For the selected model, mu and y are bounded",
            if (g11test) "\n  - g11inv must convert the data back to the original scale",
            "\n Please select another link"
          )
        )
      }
    }
  }

  # models that do not have nu or nu is fixed
  nonu <- model %in% .current.models("nonu")
  if (nonu || !endsWith(model, "V")) {
    return(TRUE)
  }

  lk <- .get.link.info(c(link$g2, link$g21, link$g22, link$g23))

  #  For all implemented models, 0 < nu < inf.
  #   - g2 cannot depend on upper bounds
  if (!lk[lk$link == link$g2, "allowed.only.lower.g2"]) {
    .stop.with.message(
      paste0(
        " g2 = ", link$g2, " is not compatible",
        "\n For the selected model 0 < nu < Inf.",
        "\n The selected link depends on upper bounds.",
        "\n Please select another link:",
        "\n  - g2 cannot depend on upper bounds",
        "\n  - ginv must map the values of g2(nu) back to (0, Inf)."
      )
    )
  }

  # check if g21 and g22 are compatible with g2
  #  g2 : (0,Inf) -> (c,Inf)
  #     g21 must take (c,Inf) to (-Inf, Inf)
  #     g22 can be any function accepting only lower bounds.
  if (lk[lk$link == link$g2, "from.lower.to.bounds"]) {
    if (!lk[lk$link == link$g21, "allowed.bounded.gi1"]) {
      # default link for beta takes (0,Inf) to (0,1)
      # so g21 must be able to deal with bounded intervals
      .stop.with.message(
        paste0(
          " The selected links are not compatible.",
          "\n g2 = ", link$g2, ", g21 = ", link$g21,
          "\n For the selected model, 0 < g2(nu) < 1, so that",
          "\n  - g21inv must convert the data back to the original scale",
          "\n Please select another link"
        )
      )
    }
  } else if (lk[lk$link == link$g2, "preserve.bounds"]) {
    g21test <- !lk[lk$link == link$g21, "allowed.only.lower.gi1"]
    pb <- lk[lk$link == link$g21, "preserve.bounds"]
    c1 <- (g21test || pb) && !allow.any[2]
    g22test <- !lk[lk$link == link$g22, "allowed.only.lower.gi2"]
    if (c1 || g22test) {
      .stop.with.message(
        paste0(
          " The selected links are not compatible.",
          "\n g2 = ", link$g2, ", g21 = ", link$g21, ", g22 = ", link$g22,
          "\n For the selected model, 0 < g2(nu) < Inf, so that",
          if (g21test) "\n  - g21 cannot depend on upper bounds",
          if (!g21test && pb) "\n  - g21 must map (0,Inf) to (-Inf, Inf)",
          if (g22test) "\n  - g22 cannot depend on upper bounds",
          "\n Please select another link"
        )
      )
    }
  } else {
    g21test <- !lk[lk$link == link$g21, "allowed.unbounded.gi1"]
    g22test <- !lk[lk$link == link$g22, "allowed.unbounded.gi1"]
    if (g21test || g22test) {
      .stop.with.message(
        paste0(
          " The selected links are not compatible.",
          "\n g2 = ", link$g2, ", g21 =", link$g21, ", g22 = ", link$g22,
          "\n For the selected model, -Inf < g2(nu) < Inf, so that",
          if (g21test) "\n  - g21 cannot depend on lower and/or upper bounds",
          if (g22test) "\n  - g22 cannot depend on lower and/or upper bounds",
          "\n Please select another link"
        )
      )
    }
  }

  # check if g23 is compatible with escale
  if (link$g13 > 0 && !lk[lk$link == link$g23, "allowed.unbounded.gi1"]) {
    .stop.with.message(
      paste0(
        " The selected link g23 is not compatible with the error.scale:",
        "\n error.scale = ", link$g13, " and  g23 =", link$g23,
        "\n  - For the selected scale, -Inf < e1(t) < Inf",
        "\n  - g23 cannot depend on lower and/or upper bounds",
        "\n Please select another link."
      )
    )
  }

  return(TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Checks compatibility
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.link <- function(model, link, lk.default, allow.any) {
  # check and convert the link to the new list format
  linkg <- .check.link.format(link$linkg)
  link$configs.linkg$ctt <- .check.link.format(
    link = link$configs.linkg$ctt, is.config = TRUE
  )
  link$configs.linkg$power <- .check.link.format(
    link = link$configs.linkg$power, is.config = TRUE
  )
  if (model == "BARC") {
    linkg[["h"]] <- link$linkh
    a <- link$configs.linkh[["ctt"]]
    if (!is.null(a)) link$configs.linkg$ctt[["h"]] <- a
    b <- link$configs.linkh[["power"]]
    if (!is.null(b)) link$configs.linkg$power[["h"]] <- b
  }

  # check if all elements are present. Set missing elements to default values
  link <- .check.link.elements(
    model = model, link = linkg, configs = link$configs.linkg,
    lk.default = lk.default$linkg, cf.default = lk.default$configs.linkg,
    escale = link$error.scale
  )
  configs <- link$configs
  link <- link$link

  # update name, if necessary (based on BASE model)
  w <- which(link == "default")
  if (length(w) > 0) {
    base <- .get.base.model(model)
    nome <- switch(base,
      BARFIMA = "SIP1",
      GARFIMA = "SIP0",
      "polynomial"
    )
    ctt <- switch(base,
      BARFIMA = 1,
      GARFIMA = 0,
      1
    )
    link[w] <- nome
    configs$ctt[[w]] <- ctt
    configs$power[[w]] <- 1
  }

  w <- which(link == "SIP" & configs$ctt == 0)
  if (length(w) > 0) link[w] <- "SIP0"
  w <- which(link == "SIP" & configs$ctt == 1)
  if (length(w) > 0) link[w] <- "SIP1"

  # check if the link can be evaluated for the current model
  ok <- .check.link.eval(model = model, link = link, allow.any = allow.any)

  return(list(link = link, configs = configs))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function: Used to check if theta is compatible with the map
# selected by the user
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  1 = (kx)(mod 1).  k integer
#  2 = Rafael's map.  0 <= theta <= 1
#  3 = logistic map. 0 <= theta  <= 4
#  4 = Manneville-Pomeau. 0 < theta < 1
#  5 = Lasota-Mackey's map. No theta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.map <- function(map, theta) {
  if ((map != 5) && is.null(theta)) {
    .stop.with.message(" 'theta' is missing with no default")
  }
  r <- .barc.maps[map, "r"]

  if (length(theta) > r && r == 0) {
    .warning.with.message(" Parameter 'theta' will be ignored.")
  }

  if (r == 0) {
    return(invisible(list(map = as.integer(map), theta = 0, r = 0L)))
  }

  # checking the length of theta
  if (length(theta) != r) {
    .stop.with.message(paste0(
      "Length of the provided 'theta' = ", length(theta),
      "\n  - The selected map requires length ", r, ".",
      "\n  - Please provide a new theta with length ", r, "."
    ))
  }

  # checking the range of theta
  # !!! if theta is a vector this part of the code will need revision !!!!!
  if (theta < .barc.maps[map, "lower"] || theta > .barc.maps[map, "upper"]) {
    .stop.with.message(" Theta is out of range")
  }

  invisible(list(map = as.integer(map), theta = theta[1:r], r = as.integer(r)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Checks the lower and upper limits for theta in BARC models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.lu.theta <- function(map, lower, upper) {
  maps <- .barc.maps

  # checking if lower and upper limits are provided
  fix.lower <- fix.upper <- FALSE
  if (is.null(lower)) fix.lower <- TRUE
  if (is.null(upper)) fix.upper <- TRUE

  # checking if lower and upper are in the correct range
  if (!is.null(lower)) {
    if (lower < maps[map, "lower"] || lower > maps[map, "upper"]) {
      fix.lower <- TRUE
    }
  }
  if (!is.null(upper)) {
    if (upper < maps[map, "lower"] || upper > maps[map, "upper"]) {
      fix.upper <- TRUE
    }
  }

  # if needed, fix the wrong values
  if (fix.lower) lower <- maps[map, "lower"]
  if (fix.upper) upper <- maps[map, "upper"]

  if (lower > upper) {
    .stop.with.message(" Please check the lower and upper limits for theta")
  }

  if (lower == upper && map != 5) {
    msg <- "lower and upper limits for theta are the same.\n "
    msg <- paste0(msg, "The range for the selected map is,\n")
    .warning.with.message(paste0(
      msg, "lower = ", maps[map, "lower"], "\n",
      "upper = ", maps[map, "upper"]
    ))
  }
  invisible(list(lower = lower, upper = upper))
}
