##----------------------------------------------------------
##
##    Simulation functions
##
##----------------------------------------------------------

##-------------------------------------------------------------------------
## internal function.
## Performs several checks to make sure that
## the correct type of variables will be passed to FORTRAN
##-------------------------------------------------------------------------
.sim.configs <- function(model, xreg, y.start, xreg.start, linkg, n, burn,
                         coefs, xregar, error.scale, seed, rngtype, y.default){

  if(missing(seed)) seed = NULL
  if(is.null(seed)) seed = stats::runif(1, min = 1000, max = 10000)

  ##-------------------
  ## initial values:
  ##-------------------
  if(length(coefs$beta) == 0) xregar = FALSE
  out <- .data.start.convert(y.start = y.start, xreg.start = xreg.start,
                             nreg = length(coefs$beta), xregar = xregar,
                             y.default = y.default)

  ##--------------------------------------------------------
  ## final sample size and burn-in  (must be integer)
  ##--------------------------------------------------------
  out$n <- as.integer(n)
  out$burn <- as.integer(burn)

  ##--------------------------------------------------------
  ## link function for mu and y  (must be integer)
  ##--------------------------------------------------------
	dummy <- .link.check(model = model, link = linkg)
  out$linkg <- .link.convert(linkg)
  if(any(is.na(out$linkg)))
    stop(paste("at least one link requested is not implemented yet", sep = ""))
  ## if only one link is provided, uses the same link for mu and y.
  if(length(linkg) == 1) out$linkg <- c(out$linkg, out$linkg)

  ##-----------
  ## coefs
  ##-----------
  pars <- c("alpha", "beta", "phi", "nu")
  out[pars] <- coefs[pars]
  if(is.null(coefs$nu)) stop("nu is missing with no default")
  ##---------------------------------------
  ##  if alpha is missing, set alpha = 0
  ##---------------------------------------
  if(is.null(out$alpha)) out$alpha <- 0
  ##--------------------------------------------------
  ##  xreg and beta - forecast is not required here
  ##  if xreg = NULL, set beta = 0 (dummy)
  ##--------------------------------------------------
  cr <- .xreg.convert(xreg = xreg, xnew = NULL, n = burn + n,
                      nnew = 0, skip.forecast = TRUE)
  out[names(cr)] <- cr
  if(!is.null(xreg)){
    if(nrow(out$xreg) < burn+n) stop("nrow(xreg) < burn + n")
    if(is.null(out$beta)) stop("beta is missing with no default")
    if(length(out$beta) < out$nreg) stop("beta and nreg are not compatible")
  }else{out$beta <- 0}
  ##-----------------------------------------
  ## if phi is missing, set phi = 0 (dummy)
  ##-----------------------------------------
  out$p <- length(coefs$phi)
  if(is.null(out$phi)) out$phi <- 0
  ##--------------------------------------------------------------
  ## xregar (if FALSE, xreg is not included in the AR recursion)
  ##--------------------------------------------------------------
  out$xregar <- as.integer(xregar)
  ##--------------------------------------------------------------
  ##  error.scale  = 1 => e = g(y) - g(mu)
  ##               = 0 => e = y - mu
  ##--------------------------------------------------------------
  out$error.scale <- as.integer(error.scale)

  ##----------------------------------------------------
  ## setting the seed and rngtype  (must be integer)
  ##----------------------------------------------------
  out$seed <- .seed.start(seed = seed, rngtype = rngtype)
  out$rngtype <- as.integer(rngtype)

  ## theta parameter for BARC model is checked in check.configs.barc
  ## this function will be called by the function specific to BARC models
  if(model == "BARC") return(invisible(out))

  ##-------------------------------
  ## theta and q - MA component
  ##-------------------------------
  out$theta <- coefs$theta
  out$q <- length(coefs$theta)
  if(is.null(out$theta)) out$theta <- 0
  ##------
  ## d
  ##------
  out$d <- ifelse(is.null(coefs$d), 0, coefs$d)

  invisible(out)
}



##---------------------------------------------------------------------------
## internal function:
## Interface between R and FORTRAN
## Also used to summarize the results of the simulation and return
## only the relevant variables
##---------------------------------------------------------------------------
.btsr.sim <- function(model, inf, configs, complete, debug){

  if(abs(configs$d) > 0 & inf < 100){
    warning(paste("non-zero d and inf = ", inf,
                  ". Be carefull, this value may be too small",
                  sep = ""), immediate. = TRUE)}

  fun <- .check.model(model[1],"sim")

  out <- .Fortran(fun,
                  n = configs$n,
                  burn = configs$burn,
                  pdist = configs$nu,
                  alpha = configs$alpha,
                  nreg = configs$nreg,
                  beta = configs$beta,
                  p = configs$p,
                  phi = configs$phi,
                  q = configs$q,
                  theta = configs$theta,
                  d = configs$d,
                  linkg = configs$linkg,
                  xreg  = configs$xreg,
                  xregar = configs$xregar,
                  yt = numeric(configs$n+configs$burn),
                  ystart = configs$y.start,
                  xstart = configs$xreg.start,
                  mut = numeric(configs$n+configs$burn),
                  etat = numeric(configs$n+configs$burn),
                  error = numeric(configs$n+configs$burn),
                  escale = configs$error.scale,
                  ns = length(configs$seed),
                  seed = configs$seed,
                  rngtype = configs$rngtype,
                  inf = as.integer(inf),
                  rev = 1L)

  if(out$rev == 1){
    warning("Revision Required. Try changing the link functions\n", immediate. = TRUE)
    return(invisible(out))
  }

  ##-----------------------------------------------
  ## if complete = TRUE returns the full model.
  ## otherwise only yt is returned
  ##-----------------------------------------------
  ini <- configs$burn + 1
  end <- configs$burn + configs$n
  if(complete){
    final <- list(model = model,
                  yt = out$yt[ini:end],
                  mut = out$mu[ini:end],
                  etat = out$eta[ini:end],
                  error = out$error[ini:end],
                  xreg = out$xreg[ini:end,])
    if(out$nreg == 0) final$xreg <- NULL
    if(debug) final$out.Fortran <- out
  }
  else final <- out$yt[ini:end]

  invisible(final)

}
