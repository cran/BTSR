##----------------------------------------------------------------------------
## internal function.
## Not supposed to be called by the user.
## This function checks if the required model is implemented
##----------------------------------------------------------------------------
.check.model <- function(model, type){
  if(!model %in% c("BARFIMA", "GARFIMA", "KARFIMA", "UWARFIMA"))
    stop("The selected model is not implemented yet")

  if(!type %in% c("sim", "fit", "extract"))
    stop("wrong type. Must be one of: sim, extract, fit")

  fun <-  switch(EXPR = model,
                 BARFIMA = "barfimaR",
                 GARFIMA = "garfimaR",
                 KARFIMA = "karfimaR",
                 UWARFIMA = "uwarfimaR")

  tp <- NULL
  if(type == "sim") tp <- "sim"

  fun <- paste0(tp,fun)
  return(fun)

}

##----------------------------------------------------------------------------
## internal function.
## Not supposed to be called by the user.
## This function is used to avoid problems if the user defines one of
## these variable as NULL.
##----------------------------------------------------------------------------
.fix.null.configs <- function(coefs, lags, fixed.values, fixed.lags, lower, upper){

  # default settings
  out <- list(coefs = NULL, lags = list(), fixed.lags = list(),
              fixed.values = list(), lower = list(), upper = list())

  # updating values passed by the user
  if(length(coefs) > 0) out$coefs[names(coefs)] <- coefs
  if(length(lags) > 0) out$lags[names(lags)] <- lags
  if(length(fixed.values) > 0) out$fixed.values[names(fixed.values)] <- fixed.values
  if(length(fixed.lags) > 0) out$fixed.lags[names(fixed.lags)] <- fixed.lags
  if(length(lower) > 0) out$lower[names(lower)] <- lower
  if(length(upper) > 0) out$upper[names(upper)] <- upper
  return(out)
}



##-------------------------------------------------------------------------
## internal function.
## Not supposed to be called by the user.
## Initializes xreg and xnew to pass to FORTRAN.
## This function creates dummy matrices of size 1 x 1 when
## the model does not have regressors or forecast is not required.
##-------------------------------------------------------------------------
.xreg.convert <- function(xreg, xnew, n, nnew, skip.forecast){

  out <- c()
  out$nnew <- as.integer(nnew)

  ##------------------------------
  ##  nreg = 0
  ##------------------------------
  if(is.null(xreg)){
    ##----------------------------
    ## no regressors in the model
    ## xnew will be ignored
    ##----------------------------
    out$xreg <- matrix(0, ncol = 1, nrow = max(1,n))
    out$nreg <- 0L
    ##----------------------------------
    ## checking if forecast is required
    ##----------------------------------
    if(skip.forecast){
      # for compatibility with other functions
      out$xnew <- matrix(0, ncol = 1, nrow = 1)
      out$nnew <- 0L
    }else{
      out$xnew <- matrix(0, ncol = 1, nrow = max(1,nnew))}
    return(invisible(out))
  }

  ##---------------------------
  ## nreg > 0
  ##---------------------------
  out$xreg <- as.matrix(xreg)
  out$nreg <- ncol(out$xreg)
  if(nrow(out$xreg) != n)
    stop("xreg and y do not have the same number of observations")

  ##----------------------------------
  ## checking if forecast is required
  ##----------------------------------
  if(skip.forecast){
    # for compatibility with other functions
    out$xnew <- matrix(0, ncol = 1, nrow = 1)
    out$nnew <- 0L
    return(invisible(out))
  }

  ##----------------------------------
  ## if forecast is required
  ##----------------------------------
  if(is.null(xnew)){
    out$nnew <- as.integer(nnew)
    out$xnew <- matrix(0, ncol = 1, nrow = max(1,nnew))
  }else{
    xnew <- as.matrix(xnew)
    if(ncol(xnew) != out$nreg){
      ##------------------------------------------------------
      ## if xnew is not compatible with xreg
      ##------------------------------------------------------
      stop("number of columns in xnew and xreg are not the same")
    }else{
      out$xnew <- xnew
      out$nnew <- nrow(xnew)
    }
  }
  invisible(out)
}



##---------------------------------------------------------------------------
## internal function.
## Initializes y.start and xreg.start to pass to FORTRAN.
##
## The default is to set the initial values as zero.
##
## For now, if initilization is required, the user must also provide the initial
## values for these variables.
##
## To do: to implement a function to initialize y.start and xreg.start based
## on the model selected.
##---------------------------------------------------------------------------
.data.start.convert <- function(y.start, xreg.start, nreg, xregar, y.default){

  # If the model dos not include xreg in the AR regression,
  # there is no use for initial values of xreg.
  if(!xregar)  xreg.start = rep(0, max(1,nreg))
  else xreg.start = NULL

  # checking compatibility
  if(!is.null(xreg.start) & xregar)
    if(length(xreg.start) != nreg)
      stop("length(xreg.start) is not the same as nreg")

  out <- c()

  # if y.start and/or xreg.start are provided, uses this values
  # otherwise, the initial values are set as the default value
  out$y.start <- ifelse(is.null(y.start), y.default, y.start)
  if(is.null(xreg.start)) out$xreg.start <- rep(0, max(1,nreg))
  else out$xreg.start <- xreg.start

  invisible(out)
}



##-------------------------------------------------------------------------
## internal function.
## Checks if the data has any NA and if it is the the correct range.
##-------------------------------------------------------------------------
.data.check <- function(yt, lower, upper, openIC){

  out <- list()

  ##----------------------
  ## checking for NA's
  ##----------------------
  if(sum(is.na(yt)) > 0){
    out$conv <- 1
    out$message <- "NA in the data"
    warning("NA's are not allowed", immediate. = TRUE)
    return(out)
  }

  ##--------------------------------------
  ## checking if y is in the correct range
  ##--------------------------------------
  if(openIC[1]) a <- (min(yt) <= lower + .Machine$double.eps)
  else a <- (min(yt) < lower)
  if(openIC[2]) b <- (max(yt) >= upper - .Machine$double.eps)
  else b <- (max(yt) > upper)
  if (a | b){
    out$conv <- 1
    out$conv_message <- "Out of range"
    warning(paste("OUT OF RANGE. yt must be bewteen ", lower,
                  " and ", upper, sep = ""), immediate. = TRUE)
    return(out)
  }

  return(out)
}



##-------------------------------------------------------------------------
## Internal function.
## Initializes the variables with fixed values
## and fixed lags to pass to FORTRAN.
##-------------------------------------------------------------------------
.coefs.convert <- function(parname, fvalues, flags, coefs, lags, npar){
  out <- c()

  ##----------------------
  ## npar = nfix + nfit
  ##----------------------
  if(is.null(npar)) stop("npar is missing")
  # if only npar is provided, the parameter is assumed to be
  # fixed and it will be set as zero.
  if(is.null(c(fvalues, flags, coefs, lags)) & npar > 0) fvalues = rep(0, npar)

  ##---------------------------------------
  ## non-fixed values
  ## lags = position of non-fixed values
  ##---------------------------------------
  lc = length(coefs)
  ll = length(lags)
  if(!is.null(coefs) & !is.null(lags)){
    if(lc != ll) stop("coefs and lags are not compatible.
                      Please, make sure that length(coefs) = length(lags)")}

  ##---------------------------------------
  ## fixed values
  ## flags = position of fixed values
  ##---------------------------------------
  lfv = length(fvalues)
  lfl = length(flags)
  if(!is.null(fvalues) & !is.null(flags)){
    if(lfv != lfl) stop("fixed.values and fixed.lags are not compatible.
                        Please, make sure that length(fvalues) = length(flags)")}


  if(!is.null(lags) & !is.null(flags))
    if(any(lags %in% flags)) stop("lags and flags have non-empty intersection")

  ##-------------------
  ## total
  ##-------------------
  lc = ll = max(lc, ll)
  lfv = lfl = max(lfv, lfl)

  if(npar != lc + lfv) stop("values provided are not compatible.
                      Please, check if fixed values/lags and
                      non-fixed values/lags were correctly informed.")

  if(lfv == 0){
    ##--------------------------------
    ## there are no fixed values:
    ##----------------------------------
    out$flags <- 0L    ## dummy passed to FORTRAN
    out$fvalues <- 0   ## dummy passed to FORTRAN
    out$nfix <- 0L
    lags = 1:npar
  }else{
    ##----------------------------------------------------
    ## if there are fixed values
    ##----------------------------------------------------
    if(is.null(lags) & is.null(flags) & npar > 1){
      stop("cannot decide which lags must be fixed/fitted.
            Please provide lags or flags")}

    if(npar == 1) flags = 1
    ##--------------
    ## lags
    ##--------------
    all = 1:npar
    if(is.null(lags) & lc > 0) lags = all[-flags] # flags was provided
    if(is.null(flags)) flags = all[-lags]         # lags was provided
    ##----------------
    ## fixed values
    ##----------------
    if(is.null(fvalues)) fvalues <- rep(0,lfv)
    out$flags <- as.integer(flags)
    out$fvalues <- fvalues
    out$nfix <- as.integer(lfv)
  }

  ##-----------------------------------------
  ## checking for non-fixed parameter values
  ##-----------------------------------------
  if(npar == lfv) out$coefs <- NULL
  else{
    ## if the non-fixed values were not provided, set as 0.
    if(is.null(coefs)) out$coefs <- rep(0, lc)
    else out$coefs <- coefs
  }

  if(lc > 0){
    if(length(lags) == 1) out$coefsnames <- parname
    else out$coefsnames <- paste(parname,"(",lags,")",sep = "")
  }else{out$coefsnames <- NULL}

  invisible(out)
}



##-------------------------------------------------------------------------
## Internal function.
## Initializes the variables with fixed values and fixed lags
## to pass to FORTRAN.
##-------------------------------------------------------------------------
.coefs.convert.all <- function(model, coefs,lags, fixed.values, fixed.lags, p, q, nreg){

  out <- c()
  par <- NULL
  nm = NULL

  ##-------------------------------------------------------
  ##  checking for and parsing fixed and non-fixed values.
  ##-------------------------------------------------------
  ## alpha - the intercept
  out$alpha <- .coefs.convert(parname = "alpha", fvalues = fixed.values$alpha,
                              flags = NULL, coefs = coefs$alpha,
                              lags = NULL, npar = 1)
  par <- c(par, alpha = out$alpha$coefs)
  nm <- c(nm, out$alpha$coefsname)

  ## beta - coefficients associated to Xreg
  out$beta <- .coefs.convert(parname = "beta", fvalues = fixed.values$beta,
                             flags = fixed.lags$beta, coefs = coefs$beta,
                             lags = lags$beta, npar = nreg)
  par <- c(par, beta = out$beta$coefs)
  nm <- c(nm, out$beta$coefsname)
  if(nreg - out$beta$nfix > 0){
    if(length(out$beta$coefs) < nreg - out$beta$nfix)
      stop("missing some values of beta in the parameter list")}

  ## phi - AR coefficients
  out$phi <- .coefs.convert(parname = "phi", fvalues = fixed.values$phi,
                            flags = fixed.lags$phi, coefs = coefs$phi,
                            lags = lags$phi, npar = p)
  par <- c(par, phi = out$phi$coefs)
  nm <- c(nm, out$phi$coefsname)
  if(p - out$phi$nfix > 0){
    if(length(out$phi$coefs) < p - out$phi$nfix)
      stop("missing some values of phi in the parameter list")}

  ## theta - MA coefficients or map parameter in BARC models
  out$theta <- .coefs.convert(parname = "theta", fvalues = fixed.values$theta,
                              flags = fixed.lags$theta, coefs = coefs$theta,
                              lags = lags$theta, npar = q)
  par <- c(par, theta = out$theta$coefs)
  nm <- c(nm, out$theta$coefsname)
  if(p - out$theta$nfix > 0){
    if(length(out$theta$coefs) < q - out$theta$nfix)
      stop("missing some values of theta in the parameter list")}

  if(!(model == "BARC")){
    ## d - long memory parameter
    out$d <- .coefs.convert(parname = "d", fvalues = fixed.values$d,
                            flags = NULL, coefs = coefs$d,
                            lags = NULL, npar = 1)
    par <- c(par, d = out$d$coefs)
    nm <- c(nm, out$d$coefsname)
  }

  ## nu - dispersion parameter
  out$nu <- .coefs.convert(parname = "nu", fvalues = fixed.values$nu,
                           flags = NULL, coefs = coefs$nu,
                           lags = NULL, npar = 1)
  par <- c(par, nu = out$nu$coefs)
  nm <- c(nm, out$nu$coefsname)

  if(!is.null(par)) names(par) <- nm
  out$coefs <- par
  out$coefsnames <- nm

  invisible(out)
}



##-------------------------------------------------------------------------
## Internal function.
## Convert the bounds to pass to FORTRAN.
##-------------------------------------------------------------------------
.bounds.convert <- function(npar, lower, upper){
  ##----------------------------------
  ##           bounds
  ##----------------------------------
  ## 0 = no bounds
  ## 1 = lower bound only
  ## 2 = lower and upper bounds
  ## 3 = upper bound only
  ##----------------------------------
  out <- c()

  out$nbd <- integer(npar)
  out$lower <- numeric(npar)
  out$upper <- numeric(npar)
  if(is.null(lower)) lower <- rep(-Inf, npar)
  if(is.null(upper)) upper <- rep(Inf, npar)

  w1 <- (lower > -Inf)&(upper == Inf)
  w2 <- (lower > -Inf)&(upper < Inf)
  w3 <- (lower == -Inf)&(upper < Inf)
  ##----------------------------
  ## lower bound only
  ##----------------------------
  if(sum(w1) > 0){
    out$nbd[w1] <- 1L
    out$upper[w1] <- 0    ## will be ignored by the function
  }
  ##----------------------------
  ## upper and lower bound
  ##----------------------------
  if(sum(w2) > 0){
    out$nbd[w2] <- 2L
    out$lower[w2] <- lower[w2]
    out$upper[w2] <- upper[w2]
  }
  ##----------------------------
  ## upper bound only
  ##----------------------------
  if(sum(w3) > 0){
    out$nbd[w3] <- 3L
    out$lower[w3] <- 0    ## will be ignored by the function
  }
  ##-------------------------------------------
  ## no bounds (FORTRAN does not accept Inf)
  ##-------------------------------------------
  w0 <- (out$nbd == 0L)
  if(sum(w0) > 0){
    out$lower[w0] <- 0
    out$upper[w0] <- 0
  }

  invisible(out)
}



