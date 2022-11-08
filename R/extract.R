##-------------------------------------------------------------------------
## internal function.
## Performs several checks to make sure that
## the correct type of variables will be passed to FORTRAN
## Returns a list with the arguments to be passed to simulation
## and/or fitting subroutines
##-------------------------------------------------------------------------
.extract.configs <- function(model, yt, y.start, y.lower, y.upper, openIC,
                             xreg, xnew, nnew, xreg.start, linkg,
                             p, q, inf, m, xregar, error.scale,
                             coefs, lags, fixed.values, fixed.lags,
                             llk, sco, info, extra){

  ##----------------------------------------------------------
  ## checking if the data has NA's or any value outside (0,1)
  ##----------------------------------------------------------
  out <- .data.check(yt = yt, lower = y.lower, upper = y.upper, openIC = openIC)
  if(!is.null(out$conv)) return(invisible(out))
  out$n <- as.integer(length(yt))
  out$y.lower = y.lower
  out$y.upper = y.upper

  ##--------------------------------------------------
  ## The code allows for different links for y and mu
  ##--------------------------------------------------
	dummy <- .link.check(model = model, link = linkg)
  if(length(linkg) == 1) linkg <- c(linkg, linkg)
  out$linkg <- .link.convert(link = linkg)

  ##--------------------------------------------------
  ## Regressors.
  ## xnew is needed by the FORTRAN subroutine so
  ## skip.forecast must be set as FALSE
  ##--------------------------------------------------
  temp <- .xreg.convert(xreg = xreg, xnew = xnew, n = out$n,
                        nnew = nnew, skip.forecast = FALSE)
  out[names(temp)] <- temp

  ##---------------------------------------------------------
  ## initial values: using y.default = y.lower -1
  ## assures that the Fortran subroutine will set g(y) = 0
  ##---------------------------------------------------------
  if(out$nreg == 0) xregar = FALSE
  temp <- .data.start.convert(y.start = y.start, xreg.start = xreg.start,
                              nreg = out$nreg, xregar = xregar,
                              y.default = y.lower - 1)
  out[names(temp)] <- temp

  ##----------------------------------------------------------------------------
  ## parameters initialization and fixed values identification
  ##----------------------------------------------------------------------------
  if(is.null(coefs$nu)){
    if(is.null(fixed.values$nu)) stop("nu is missing with no default")}

  temp <- .coefs.convert.all(model = model, p = p, q = q, nreg = out$nreg,
                             coefs = coefs, lags = lags,
                             fixed.values = fixed.values, fixed.lags = fixed.lags)
  out[names(temp)] <- temp
  out$p <- as.integer(p)

  ##-------------------------------
  ## Other configurations
  ##-------------------------------
  out$inf <- as.integer(inf)
  if(!is.null(out$d)){
    if((out$d$nfix == 0 | out$d$fvalues != 0) & out$inf < 100)
      warning(paste("non-zero d and inf = ", inf,
                    ". Be carefull, this value may be too small",
                    sep = ""), immediate. = TRUE)}

  out$m <- as.integer(m)
  out$error.scale <- as.integer(error.scale)
  out$xregar <- as.integer(xregar)

  out$llk <- as.integer(llk)
  out$sco <- as.integer(sco)
  out$info <- as.integer(info)
  out$extra <- as.integer(extra)

  if(!(model == "BARC")) out$q <- as.integer(q)
  out$npar <- length(out$coefs)
  if(out$npar == 0) out$coefs = 0

  invisible(out)
}



##---------------------------------------------------------------------------
## internal function:
## Interface between R and FORTRAN
## Also used to summarize the results of the extraction and return
## only the relevant variables
##---------------------------------------------------------------------------
.btsr.extract <- function(model, yt, configs, debug){

  fun <- .check.model(model[1], "extract")

  if(configs$npar == 0){
    configs$sco <- 0L
    configs$info <- 0L
    configs$extra <- 0L
  }
  temp <-  .Fortran(fun,
                    n = configs$n,
                    yt = yt,
                    gyt = numeric(configs$n),
                    ystart = configs$y.start,
                    nreg = as.integer(configs$nreg),
                    xreg = configs$xreg,
                    xstart = configs$xreg.start,
                    mut = numeric(configs$n),
                    etat = numeric(configs$n),
                    error = numeric(configs$n),
                    escale = configs$error.scale,
                    nnew = configs$nnew,
                    xnew = configs$xnew,
                    ynew = numeric(max(1,configs$nnew)),
                    linkg = configs$linkg,
                    npar = max(1L, configs$npar),
                    coefs = configs$coefs,
                    fixa = configs$alpha$nfix,
                    alpha = configs$alpha$fvalues,
                    fixb = configs$beta$nfix,
                    flagsb = configs$beta$flags,
                    beta = configs$beta$fvalues,
                    p = configs$p,
                    fixphi = configs$phi$nfix,
                    flagsphi = configs$phi$flags,
                    phi = configs$phi$fvalues,
                    xregar = configs$xregar,
                    q = configs$q,
                    fixtheta = configs$theta$nfix,
                    flagstheta = configs$theta$flags,
                    theta = configs$theta$fvalues,
                    fixd = configs$d$nfix,
                    d = configs$d$fvalues,
                    fixnu = configs$nu$nfix,
                    pdist = configs$nu$fvalues,
                    inf = configs$inf,
                    m = configs$m,
                    llk = configs$llk,
                    sll = 0,
                    sco = configs$sco,
                    U = numeric(max(1,configs$npar*configs$sco)),
                    info = configs$info,
                    K = diag(max(1,configs$npar*configs$info)),
                    extra = configs$extra,
                    Drho = matrix(0, max(1,configs$n*configs$extra),
                                  max(1,(configs$npar-1+configs$nu$nfix)*configs$extra)),
                    T = numeric(max(1, configs$n*configs$extra)),
                    E = matrix(0, max(1,configs$n*configs$extra),
                               1+2*(1-configs$nu$nfix)*configs$extra),
                    h = numeric(max(1, configs$n*configs$extra)))

  out <- list(model = model)
  vars <- c("coefs","yt", "xreg", "gyt", "mut", "etat", "error")
  out[vars] <- temp[vars]
  if(configs$nreg == 0) out$xreg = NULL

  if(configs$llk == 1) out$sll <- temp$sll
  if(configs$sco == 1){
    out$score <- temp$U
    names(out$score) <- names(configs$coefs)
  }
  if(configs$info == 1){
    out$info.Matrix <- as.matrix(temp$K)
    colnames(out$info.Matrix) <- names(configs$coefs)
    rownames(out$info.Matrix) <- names(configs$coefs)
  }
  if(configs$extra == 1){
    out[c("Drho", "T", "E", "h")] <- temp[c("Drho", "T", "E", "h")]
  }
  if(configs$nnew > 0) out$yt.new <- temp$ynew
  if(debug) out$out.Fortran <- temp
  invisible(out)
}
