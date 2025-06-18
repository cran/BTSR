#' @title
#' Retrieve Default Arguments for BTSR Package Functions
#'
#' @description
#' Extracts and displays the default argument values for any function in the
#' BTSR package, including both exported and non-exported functions.
#'
#' @param fun Character string specifying the function name to examine
#'   (case-sensitive).
#'
#' @return Invisibly returns a data frame with two columns:
#' \itemize{
#'   \item \code{Argument}: Name of the function parameter
#'   \item \code{Default}: Default value or `"(no default)"` string if no
#'   default exists
#' }
#'   The function primarily prints a formatted table of the results to the
#'   console.
#'
#' @examples
#' \dontrun{
#' # View defaults for BTSR.fit function
#' get.defaults("BARFIMA.fit")
#'
#' # Capture the results for later use
#' defaults <- get.defaults("BARFIMA.fit")
#' }
#'
#' @export
#' @importFrom utils getFromNamespace
get.defaults <- function(fun) {
  ns <- getNamespace("BTSR")
  f <- try(get(fun, envir = ns), silent = TRUE)

  if (inherits(f, "try-error")) {
    stop(fun, " is not a function in the BTSR package")
  }
  if (!is.function(f)) {
    stop(fun, " is not a function")
  }

  args <- formals(f)

  check <- c(
    endsWith(fun, ".sim"),
    endsWith(fun, ".extract"),
    endsWith(fun, ".fit")
  )
  if (any(check) & !startsWith(fun, "BARC") & !startsWith(fun, "btsr")) {
    check2 <- c(
      endsWith(fun, "V.sim"),
      endsWith(fun, "V.extract"),
      endsWith(fun, "V.fit")
    )

    if (!startsWith(fun, "BARFIMA") || any(check2)) {
      w <- which(check)
      f2 <- paste0("BARFIMA", c(".sim", ".extract", ".fit")[w])
      f2 <- get(f2, envir = ns)

      args2 <- formals(f2)
      args2[names(args)] <- args
      args <- args2
      if (w %in% c(1, 2)) args <- args[names(args) != "..."]
    }
  }

  format_default <- function(x) {
    if (is.name(x) && as.character(x) == "") {
      return("(no default)")
    }

    if (is.list(x)) {
      if (length(x) == 0) {
        return("list()")
      }

      elements <- character(length(x))
      for (i in seq_along(x)) {
        nm <- if (is.null(names(x))) "" else paste0(names(x)[i], " = ")
        val <- if (is.list(x[[i]])) "(list)" else deparse(x[[i]])
        elements[i] <- paste0(nm, val)
      }
      return(paste0("list(", paste(elements, collapse = ", "), ")"))
    }

    txt <- paste(deparse(x), collapse = " ")
    gsub("\\s+", " ", txt)
  }

  arg_table <- data.frame(
    Argument = names(args),
    Default = sapply(args, format_default),
    stringsAsFactors = FALSE
  )

  # Calculate column widths
  arg_width <- max(nchar("Argument"), nchar(arg_table$Argument))
  def_width <- nchar("Default")

  # Print header with optimized underline length
  cat("\n")
  cat(format("Argument", width = arg_width), " ", "Default\n")
  cat(rep("-", arg_width), "  ", rep("-", def_width + 2), "\n", sep = "")

  # Print each argument
  for (i in seq_len(nrow(arg_table))) {
    arg <- format(arg_table$Argument[i], width = arg_width)
    def <- arg_table$Default[i]

    if (def == "list()") {
      cat(arg, " : list()\n", sep = "")
      next
    }

    if (grepl("^list\\(", def)) {
      elements <- strsplit(sub("^list\\((.*)\\)$", "\\1", def), ", ")[[1]]

      if (length(elements) == 1) {
        cat(arg, " : list(", elements[1], ")\n", sep = "")
      } else {
        cat(arg, " : list(", elements[1], ",\n", sep = "")
        for (j in 2:(length(elements) - 1)) {
          cat(rep(" ", arg_width), "        ", elements[j], ",\n", sep = "")
        }
        cat(rep(" ", arg_width), "        ", elements[length(elements)], ")\n", sep = "")
      }
    } else {
      def_lines <- strwrap(def,
        width = getOption("width") - arg_width - 3,
        indent = 0, exdent = 0
      )
      cat(arg, " : ", def_lines[1], "\n", sep = "")
      if (length(def_lines) > 1) {
        for (line in def_lines[-1]) {
          cat(rep(" ", arg_width), " : ", line, "\n", sep = "")
        }
      }
    }
  }

  invisible(arg_table)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Performs several checks to make sure that
# the correct type of variables will be passed to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.model.output <- function(model, obj, type, configs, debug) {
  is.sim <- type == "sim"
  is.fit <- type == "fit"
  is.barc <- model == "BARC"

  # update configs values
  nm <- c("lconfig")
  if (!is.barc) nm <- c("order", nm)
  configs[nm] <- obj[nm]

  # check if there was any problem in FORTRAN
  if ((!is.fit && obj$conv > 0) || (is.fit && obj$conv > 2)) {
    obj$model <- model
    obj$message <- .code.to.message(
      code = obj$conv, method = configs$control$method
    )
    # if type = "fit" and code = 1, do not break the code
    class(obj) <- c("try-error", class(obj))
    return(invisible(obj))
  }

  # define start and end for each time series based on the type of object
  ini <- ifelse(is.sim, configs$burn + 1, 1)
  end <- ifelse(is.sim, configs$burn + configs$n, configs$n)

  yt <- obj$ts[ini:end, "yt"]
  if (is.sim && !configs$complete) {
    # if complete = TRUE returns the full model.
    #   - for the general model, there is nothing to do
    #   - if nu is fixed, part2 is constant. In this case, fix the output for
    #     compatibility with previous versions
    # otherwise only yt is returned.

    class(yt) <- c(class(yt), tolower(.get.base.model(model)))
    return(invisible(yt))
  }

  out <- list(model = model)

  debug <- ifelse(is.null(debug), FALSE, debug)
  if (debug) {
    on.exit(
      {
        out$out.Fortran <- obj
        return(out)
      },
      add = TRUE
    )
  }

  class(out) <- c(class(out), tolower(.get.base.model(model)))
  if (is.fit) {
    class(out) <- c(class(out), "btsr")
    if (is.barc) {
      d <- FALSE
    } else {
      d <- !(configs$d$nfix == 1 & configs$d$fvalues == 0)
    }
    out$call <- .print.call(
      model = model, p = configs$order[, "p"], q = configs$order[, "q"],
      d = d, nreg = configs$order[, "nreg"]
    )
    out$n <- configs$n
  }

  if (is.sim || is.fit) out$yt <- yt

  if (is.sim && is.barc) {
    out$u0 <- ifelse(configs$burn == 0, obj$u0, obj$ts[configs$burn + 1, "Tt"])
  }

  nm <- .get.output.names(model = model)
  if (is.fit) {
    # g11(y) and g12(y)
    l1 <- c(obj$link["g11"], obj$lconfig["g11", ])
    l2 <- c(obj$link["g12"], obj$lconfig["g12", ])
    if (identical(unname(l1), unname(l2))) {
      nms <- "g11yt"
    } else {
      nms <- c("g11yt", "g12yt")
    }
    out$gyt <- obj$ts[, nms]

    # regressors
    if (sum(obj$order[, "nreg"] > 0)) {
      if (length(nm$xreg) == 1) {
        out$xreg <- obj$xreg1
      } else {
        w <- obj$order[, "nreg"] > 0
        out$xreg[.parts[w]] <- as.list(obj[nm$xreg[w]])
      }
    } else {
      configs$xreg <- NULL
    }

    out$control <- configs$control

    #  Convergence
    out$convergence <- obj$conv
    out$message <- .code.to.message(code = obj$conv, method = configs$control$method)
    out$counts <- obj$neval

    if (sum(configs$npar) > 0) {
      #  Coefficients: starting values
      out$start <- configs$coefs

      #  Coefficients: final values
      out$coefficients <- obj$coefs
    }
  }

  if (!is.sim) {
    # g21(g(nu)) and g22(g(nu))
    l1 <- c(obj$link["g21"], obj$lconfig["g21", ])
    l2 <- c(obj$link["g22"], obj$lconfig["g22", ])
    if (!identical(unname(l1), unname(l2))) {
      out$g22varthetat <- obj$ts[, "g22varthetat"]
    }
  }

  if (is.fit) {
    # returns a matrix with the time series
    out$fitted.values <- obj$ts[, nm$fitted]
  } else {
    # returns an entry for each time series
    out[nm$fitted] <- as.list(as.data.frame(obj$ts[ini:end, nm$fitted]))
  }
  out$etat <- obj$ts[ini:end, nm$eta]

  # error term
  out$error <- obj$ts[ini:end, nm$error]

  if (is.sim) {
    return(out)
  }

  if (is.fit) {
    # residuals
    if (obj$link["g13"] == 1) {
      out$residuals <- obj$ts[ini:end, "yt"] - obj$ts[ini:end, "mut"]
    } else {
      out$residuals <- obj$ts[ini:end, "error1"]
    }
  }

  # forecast
  if (obj$length[2] > 0) {
    out$forecast <- obj$forecast[, nm$forecast]
    if (sum(obj$order[, "nreg"] > 0)) {
      if (length(nm$xnew) == 1) {
        out$xnew <- obj$xnew1
      } else {
        w <- obj$order[, "nreg"] > 0
        out$xnew[.parts[w]] <- as.list(obj[nm$xnew[w]])
      }
    }
  }

  #  likelihood, information matrix and extra arguments
  out$sll <- NULL
  out$score <- NULL
  out$info.Matrix <- NULL
  if (obj$extras["llk"] == 1) out$sll <- obj$sll
  if (obj$extras["sco"] == 1) {
    out$score <- obj$U
    names(out$score) <- configs$coefsnames
  }
  if (obj$extras["info"] == 1) {
    out$info.Matrix <- as.matrix(obj$K)
    colnames(out$info.Matrix) <- configs$coefsnames
    rownames(out$info.Matrix) <- configs$coefsnames
  }
  if (obj$extras["extra"] == 1) {
    out[c("D", "T", "E", "h")] <- obj[c("D", "T", "E", "h")]
  }

  if (is.fit) {
    # information for printing summary
    out$link <- cbind(
      link = obj$link,
      ctt = obj$lconfig[, "ctt"],
      power = obj$lconfig[, "power"]
    )

    #  Extra information for prediction
    nms <- names(out)
    nmsc <- names(configs)
    nmse <- !(nmsc %in% nms)
    out$configs[nmsc[nmse]] <- configs[nmse]
    names(out)[names(out) == "yt"] <- "series"
  }

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Extract the necessary arguments to pass to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.predict.configs <- function(obj, newdata, nnew) {
  # Update the object to get the variables to pass to FORTRAN
  nms <- names(obj)
  nms <- nms[nms != "configs"]
  temp <- obj$configs
  obj <- obj[nms]
  obj[names(temp)] <- temp

  is.barc <- obj$model == "BARC"

  # update gy
  if (!is.matrix(obj$gyt)) {
    obj$gyt <- cbind(g11yt = obj$gyt, g12yt = obj$gyt)
  }

  if (!is.barc) {
    # update eta
    if (!is.matrix(obj$etat)) {
      obj$etat <- cbind(
        eta1t = obj$etat,
        eta2t = numeric(obj$n)
      )
    }

    # update fitted.values
    if (!is.matrix(obj$fitted.values)) {
      obj$fitted.values <- cbind(
        mut = obj$fitted.values,
        nut = numeric(obj$n),
        varthetat = numeric(obj$n)
      )
    }

    # update g(gnu)
    if (is.null(obj$g22varthetat)) {
      obj$g22varthetat <- obj$etat[, "eta2t"]
    }

    # error
    if (!is.matrix(obj$error)) {
      obj$error <- cbind(error1 = obj$error, error2 = numeric(obj$n))
    }
  }

  # convert nnew, if necessary
  temp <- .convert.xreg(
    model = obj$model, xreg = obj$xreg, xnew = newdata,
    n = obj$n, nnew = nnew, skip.forecast = FALSE
  )
  obj[c("nnew", "xreg", "xnew")] <- temp[c("nnew", "xreg", "xnew")]
  obj[c("llk", "sco", "info", "extra")] <- 0L
  names(obj)[names(obj) == "coefficients"] <- "coefs"

  return(obj)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Performs several checks to make sure that
# the correct type of variables will be passed to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @importFrom stats setNames
.get.model.configs <- function(model, mc, type) {
  is.sim <- type == "sim"
  is.fit <- type == "fit"
  is.barc <- model == "BARC"

  # Get default values and codes to pass to FORTRAN
  #  - model code,
  #  - distribution parameters and limits (need update)
  #  - limits and constants for the links (need updae)
  out <- .get.model.settings(model = model)

  # extracts the arguments from the object constructed from "..."
  args <- .get.arguments(type = type)
  pdist <- .build.list(mc = mc, args = args$pdist)
  series <- .build.list(mc = mc, args = args$series)
  order.configs <- .build.list(mc = mc, args = args$order.configs)
  data.start <- .build.list(mc = mc, args = args$data.start)
  link.configs <- .build.list(mc = mc, args = args$link.configs)
  coefs.configs <- .build.list(mc = mc, args = args$coefs.configs)
  if (!is.sim) {
    extra.configs <- .build.list(mc = mc, args = args$extra.configs)
    control <- mc$control
  }

  # Updating the distribution related parameters
  out$pdist <- .check.dist.parameters(
    model = out$model, p.names = out$p.names,
    pdist = out$pdist, param.list = pdist
  )

  # update default map value for BARC models
  if (is.barc) {
    if (!is.null(mc$map)) out$map <- mc$map
    order.configs$q <- .barc.maps$r[.barc.maps$map == out$map]
  }

  # Setting yt related configurations
  if (is.sim) {
    # - final sample size and burn-in (must be integer)
    out$n <- ifelse(is.null(series$n), 1L, as.integer(series$n))
    out$burn <- ifelse(is.null(series$burn), 0L, as.integer(series$burn))
  } else {
    # - check if yt is provided
    if (is.null(series$yt)) .stop.with.message(" yt is missing! ")
    out$n <- length(series$yt)
    out$burn <- 0L

    # - checking if the data has NA's or any out-of-range value
    data.ok <- .check.data(
      yt = series$yt, openIC = c(TRUE, TRUE),
      lower = out$pdist[["y.lower"]], upper = out$pdist[["y.upper"]]
    )
  }
  N <- out$n + out$burn

  # Regressors.
  #  for simulation:
  #    - forecast is not required
  #    - xnew is not required
  #  for extract and fit:
  #    - xnew is needed by the FORTRAN subroutine
  #    - skip.forecast must be set as FALSE
  skip.forecast <- is.sim
  temp <- .convert.xreg(
    model = model, n = N, nnew = series$nnew,
    xreg = series$xreg, xnew = series$xnew,
    skip.forecast = skip.forecast
  )
  out[c("nnew", "xreg", "xnew")] <- temp[c("nnew", "xreg", "xnew")]
  order.configs$nreg <- temp$nreg

  # Parameters configuration
  ignore.coef <- FALSE

  if (is.fit) {
    # - check if starting values must be ignored
    if (is.null(coefs.configs$ignore.start)) {
      coefs.configs$ignore.start <- FALSE
    }

    # - in case starting values are provided and must be ignored, set
    #   the coefs list to NULL. Save starting values in the coefs list.
    if (coefs.configs$ignore.start) coefs.configs$start <- NULL
    coefs.configs$coefs <- coefs.configs$start
    ignore.coef <- is.null(coefs.configs$start)

    # - check if d = TRUE/FALSE is provided
    if (is.null(order.configs$d)) order.configs$d <- out$fitd
    if (length(order.configs$d) == 1) {
      order.configs$d <- c(order.configs$d, order.configs$d)
    }
  }

  # - set default values for some models (e.g. iid)
  # - the default is to check the coefs list, but when fitting a
  #   model, starting values are computed only after checking the
  #   format for other coefficient related lists
  # - identify the order of the model
  temp <- .check.coefs(
    model = model, coefs.configs = coefs.configs, order = order.configs,
    type = type, ignore.coef = ignore.coef
  )
  coefs.configs[names(temp)] <- temp

  # update d = TRUE/FALSE after checking user provided starting values
  # (in case they exist) and fixed values
  out$fitd <- coefs.configs$fitd
  out$order <- coefs.configs$order
  xregar <- coefs.configs$xregar
  out$xregar <- as.integer(xregar)

  # check compatibility (only for simulation)
  if (is.sim) {
    betaok <- .check.size.beta(
      beta = lapply(
        coefs.configs$coefs[.parts], `[[`, "beta"
      ),
      nreg = out$order[, "nreg"]
    )
  }

  # Link configurations
  #  - check link and update defaults (if necessary)
  #  - save the list with link names
  #  - convert to FORTRAN format
  if (model == "BARC") {
    allow.any <- TRUE
  } else {
    allow.any <- apply(
      out$order[, c("nreg", "p", "q")], 1, function(x) sum(x) == 0
    )
  }
  temp <- .check.link(
    model = model,
    link = link.configs,
    lk.default = list(
      linkg = out$linkg,
      configs.linkg = list(
        ctt = out$lconfig[, "ctt"], power = out$lconfig[, "power"]
      )
    ),
    allow.any = allow.any
  )
  linkg <- temp$link
  configs.linkg <- temp$configs
  out$linkg <- c(
    g13 = .convert.escale(escale = linkg$g13),
    .convert.link(linkg[names(linkg) != "g13"])
  )
  out$lconfig[names(temp$configs$ctt), "ctt"] <- temp$configs$ctt
  out$lconfig[names(temp$configs$power), "power"] <- temp$configs$power

  # check if the bounds for mu and y need update
  if (!all(out$pdist[1:2] == out$lconfig["g11", 1:2])) {
    out$lconfig[c("g11", "g12"), 1] <- out$pdist[["y.lower"]]
    out$lconfig[c("g11", "g12"), 2] <- out$pdist[["y.upper"]]
  }

  # Setting and converting initial values
  data.start$xreg.start <-
    .check.xreg.start(
      xreg.start = data.start$xreg.start, p = out$order[, "p"],
      xregar = xregar, xreg = out$xreg, nreg = out$order[, "nreg"]
    )
  data.start <- .convert.data.start(
    linkg = linkg, data.start = data.start, nreg = out$order[, "nreg"]
  )
  out[names(data.start)] <- data.start

  if (is.fit && ignore.coef) {
    # checking if parameter initialization is required.
    # in case no starting values were not provided, uses the default values.
    # in case ignore.start = TRUE, starting values are ignored and recalculated.
    # partial starting values are not allowed.
    series$xreg <- setNames(vector("list", length(.parts)), .parts)
    for (i in 1:length(.parts)) {
      if (out$order[i, "nreg"] > 0) {
        series$xreg[[.parts[i]]] <- out$xreg[[.parts[i]]]
      }
    }
    coefs.configs$coefs <- coefs.start(
      model = model, xreg = series$xreg, yt = series$yt, y.start = out$y.start,
      y.lower = out$pdist[["y.lower"]], y.upper = out$pdist[["y.upper"]],
      p = out$order[, "p"], d = out$fitd, map = out$map, q = out$order[, "q"],
      linkg = linkg[.current.link.names], configs.linkg = configs.linkg,
      lags = coefs.configs$lags, fixed.values = coefs.configs$fixed.values,
      fixed.lags = coefs.configs$fixed.lags
    )
    # for compatibility with old version, the output of coefs.start depends
    #  on the model, so conversion might be necessary
    coefs.configs$coefs <- .check.coefs.defaults(
      model = model, object = coefs.configs$coefs,
      type = "coefs", is.coef = TRUE
    )
  }

  if (is.barc) {
    # get the value of theta
    if (is.sim) {
      theta <- coefs.configs$coefs$part1$theta
    } else {
      # check if theta is fixed or not
      if (out$order[1, "q"] > 0) {
        theta <- c(
          coefs.configs$fixed.values$part1$theta,
          coefs.configs$coefs$part1$theta
        )
      } else {
        theta <- NULL
      }
    }

    # check if theta and the map are compatible.
    # Converts the map to FORTRAN format
    cf <- .check.map(map = out$map, theta = theta)
    out$map <- cf$map

    # if r = 0 use dummy values to pass to FORTRAN
    if (cf$r == 0 && out$order[, "q"] > 0) {
      out$order[, "q"] == 0L
      if (is.sim) {
        out$theta <- c(0, 0)
      } else {
        for (part in .parts) {
          for (li in c("coefs", "lags", "fixed.lags", "fixed.values")) {
            coefs.configs[[li]][[part]]["theta"] <- NULL
          }
        }
      }
    }

    # if u0 was not provided use the default value
    if (is.null(c(
      coefs.configs$coefs$part1$u0, coefs.configs$fixed.values$part1$u0
    ))) {
      coefs.configs$coefs$part1$u0 <- pi / 4
    }
  }

  # organizing the values to be passed to FORTRAN
  temp <- .convert.coefs(
    model = model, order = out$order, is.sim = is.sim,
    coefs = coefs.configs$coefs, lags = coefs.configs$lags,
    fixed.values = coefs.configs$fixed.values,
    fixed.lags = coefs.configs$fixed.lags
  )
  out[names(temp)] <- temp

  if (is.fit) {
    # Convert the bounds bound
    temp <- .convert.bounds(
      model = model, nfix = out$nfix, order = out$order,
      lower = coefs.configs$lower, upper = coefs.configs$upper,
      map = out$map
    )
    out[names(temp)] <- temp
  }

  # Other configurations
  if (!is.null(out$d)) {
    for (i in 1:2) {
      if (is.sim) {
        check <- abs(out$d[i]) > 0
      } else {
        check <- (out$d$nfix[i] == 0) || (abs(out$d$fvalues[i, 1]) > 0)
      }
      if (check && out$order[i, "inf"] < 100) {
        .warning.with.message(
          paste0(
            "\n Part ", i, ":",
            "\n d != 0 and inf = ", out$order[i, "inf"], ".",
            "\n Be carefull, this value of inf may be too small"
          )
        )
      }
    }
  }

  if (!is.sim) {
    # set defaults
    if (sum(out$npar) == 0) {
      extra.configs[c("sco", "info", "extra")] <- 0L
      out$coefs <- 0
      if (is.fit) {
        # dummy variables
        out$nbd <- 0L
        out$lower <- 0
        out$upper <- 0
      }
    }
    extra.configs <- .convert.extra.configs(extra.configs)
    out[names(extra.configs)] <- extra.configs
    if (is.fit) {
      out$control <- .update.control(control)
      if (out$control$method.code == 0L) out$sco <- 1L
    }
  }

  # remove unused from main list
  remove <- c("p.names", "error.scale")
  if (is.sim) remove <- c(remove, "xnew")
  out <- out[setdiff(names(out), remove)]
  invisible(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# This function checks the input name and finds the base model
#  - long-memory models:
#    *ARFIMAV: General model with time varying parameters.
#    *ARFIMA: Classical long memory models with nu fixed.
#
#  - short-memory models:
#    *ARMAV: Short memory models (d = 0) with time varying parameters.
#    *ARMA: Classical short memory models (d = 0) with nu fixed.
#
#  - regression models:
#    *REGV: Classical regression model with time varying parameters.
#    *REG: Classical regression model with nu fixed.
#
#  - iid samples:
#    BETA, GAMMA, etc.: iid samples.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.base.model <- function(model) {
  if (model == "BARC") {
    return(model)
  }

  msl <- .current.models("arma")
  mreg <- .current.models("reg")
  miid <- paste0(.current.models("iid"), "IID")

  model <- paste0(model[1], "IID")
  MODEL <- startsWith(
    model,
    c(
      msl, # long and short memory models
      mreg, # regression models
      miid # iid samples
    )
  )

  # update the name of the base model
  MODEL <- rep(.current.models("base"), 3)[which(MODEL == TRUE)]

  return(MODEL)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Finds the type of the model for printing purposes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.model.type <- function(model) {
  if (model == "BARC") {
    return("ARMA")
  }

  msl <- .current.models("arma")
  mreg <- .current.models("reg")
  miid <- paste0(.current.models("iid"), "IID")

  model <- paste0(model[1], "IID")
  MODEL <- startsWith(
    model,
    c(
      msl, # long and short memory models
      mreg, # regression models
      miid # iid samples
    )
  )

  # update the name of the base model
  MODEL <- c(
    rep("ARMA", length(msl)),
    rep("REG", length(mreg)),
    rep("iid", length(miid))
  )[which(MODEL == TRUE)]

  return(MODEL)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal function.
# Finds the type of the model for printing purposes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.dist.name <- function(model) {
  if (model == "BARC") {
    return("BETA")
  }

  MODEL <- startsWith(
    model,
    c(
      .current.models("arma"), # long and short memory models
      .current.models("reg")
    )
  )
  # update the name of the base model
  MODEL <- rep(.current.models("iid"), 2)[which(MODEL == TRUE)]

  return(MODEL)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Get model configurations to pass to FORTRAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.model.settings <- function(model) {
  # check if is iid before conversion
  is.iid <- model %in% .current.models("iid")
  is.barc <- model == "BARC"
  endsV <- endsWith(model, "V")

  # BARC models do not have d
  fitd <- !is.barc

  # find the base model
  if (fitd) {
    model <- .get.base.model(model)
    fitd <- c(fitd, fitd)
  }

  out <- list(
    model = model,
    fitd = fitd,
    code = .name.to.code(model) # FORTRAN code
  )
  if (is.barc) out$map <- .default.map.barc

  # Distribution related configurations
  out <- c(out, .dist.defaults(model))

  # Link related configurations
  out <- c(out, .links.default(model = model, is.iid = is.iid, endsV = endsV))

  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal function.
# Not supposed to be called by the user.
# Returns information on a given link (or list of links)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.link.info <- function(link, only.code = FALSE) {
  # in case a list of names is given
  link <- unlist(link)

  # link information
  info <- .valid.links(ifelse(only.code, "code", "all"))

  # if only codes are required (e.g. to pass to FORTRAN)
  if (only.code) {
    codes <- link %in% info$link
    for (i in 1:length(link)) {
      codes[i] <- ifelse(codes[i], info$code[which(link[i] == info$link)], NA)
    }
    return(list(codes = codes))
  }

  # if complete information is required (e.g, to perform checks)
  # here we don't need duplicates
  invisible(info[info$link %in% unique(link), ])
}
