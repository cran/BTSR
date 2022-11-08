##-------------------------------------------------------------------------
##  Similar to make.link()
##-------------------------------------------------------------------------
##'  Given the name of a link, this function returns a link function,
##'  an inverse link function, the derivative   \eqn{d\eta / d\mu}{deta/dmu}
##'  and the derivative \eqn{d\mu / d\eta}{dmu/deta}.
##'
##' @title Create a Link for BTSR models
##'
##' @param link  character; one of  \code{"linear"}, \code{"logit"},
##' \code{"log"}, \code{"loglog"}, \code{"cloglog"}. See \sQuote{Details}.
##'
##' @return An object of class \code{"link-btsr"}, a list with components
##'
##'  \item{linkfun}{Link function \code{function(mu)}}
##'  \item{linkinv}{Inverse link function \code{function(eta)}}
##'  \item{linkdif}{Derivative \code{function(mu)} \eqn{d\eta / d\mu}{deta/dmu}}
##'  \item{mu.eta}{Derivative \code{function(eta)} \eqn{d\mu / d\eta}{dmu/deta}}
##'  \item{name}{a name to be used for the link}
##'
##'@details The available links are:
##'
##'  linear: \eqn{f(x) = ax}, for \eqn{a} real.  The parameter is set using the
##'  argument \code{ctt.ll}, when invoking the functions created by \code{link.btsr}
##'
##'  logit:  \eqn{f(x) = log(x/(1-x))}
##'
##'  log:    \eqn{f(x) = log(x)}
##'
##'  loglog: \eqn{f(x) = log(-log(x))}
##'
##'  cloglog: \eqn{f(x) = log(-log(1-x))}
##'
##' @examples
##' mylink <- BTSR::link.btsr("linear")
##' y = 0.8
##' a = 3.4
##' gy = a*y
##'
##' mylink$linkfun(mu = y, ctt.ll = a); gy
##' mylink$linkinv(eta = gy, ctt.ll = a); y
##' mylink$diflink(mu = y, ctt.ll = a); a
##' mylink$mu.eta(eta = gy, ctt.ll = a); 1/a
##'
##'
##' @export
link.btsr <- function(link){

  ##--------------------------------------------------
  ## linkfun:  Link function function(mu)
  ## linkinv:  Inverse link function function(eta)
  ## mu.eta:   Derivative function(eta) dmu/deta
  ## diflink:  Derivative function(mu) deta/dmu
  ##--------------------------------------------------

  # convert character to number
  linktemp <- .link.convert(link)

  # defines g(mu)
  linkfun <- function(mu,...){
    args <- list(...)
    ctt.ll <- 1
    yl <- 0
    yu <- 1
    if(link == "linear" & !is.null(args[['ctt.ll']])) ctt.ll <- args$ctt.ll
    if(!is.null(args[['y.lower']])) yl <- max(args$y.lower, .Machine$double.xmin)
    if(!is.null(args[['y.upper']])) yu <- min(args$y.upper, .Machine$double.xmax)
    n <- length(mu)
    .Fortran("linkR", link = linktemp, a = ctt.ll, ylim = c(yl, yu),
             n = n, ilk = 0L, y = mu, lk = 1L,
             gy = numeric(n), dl = 0L, dlink = 1)$gy
  }

  # defines g^{-1}(eta)
  linkinv <- function(eta,...){
    args <- list(...)
    ctt.ll <- 1
    yl <- 0
    yu <- 1
    if(link == "linear" & !is.null(args[['ctt.ll']])) ctt.ll <- args$ctt.ll
    if(!is.null(args[['y.lower']])) yl <- max(args$y.lower, .Machine$double.xmin)
    if(!is.null(args[['y.upper']])) yu <- min(args$y.upper, .Machine$double.xmax)
    n <- length(eta)
    .Fortran("linkR", link = linktemp, a = ctt.ll,  ylim = c(yl, yu),
             n = n, ilk = 1L, y =  numeric(n), lk = 0L,
             gy = eta, dl = 0L, dlink = 1)$y
  }

  # defines dg/dmu
  diflink <- function(mu,...){
    args <- list(...)
    ctt.ll <- 1
    yl <- 0
    yu <- 1
    if(link == "linear" & !is.null(args[['ctt.ll']])) ctt.ll <- args$ctt.ll
    if(!is.null(args[['y.lower']])) yl <- max(args$y.lower, .Machine$double.xmin)
    if(!is.null(args[['y.upper']])) yu <- min(args$y.upper, .Machine$double.xmax)
    n <- length(mu)
    .Fortran("linkR", link = linktemp, a = ctt.ll,  ylim = c(yl, yu),
             n = n, ilk = 0L, y = mu, lk = 0L, gy = 1,
             dl = 1L, dlink =  numeric(n))$dlink
  }

  # defines dmu/deta = 1/g'(ginv(eta))
  mu.eta <- function(eta,...){
    1/diflink(mu = linkinv(eta = eta,...),...)
  }

  #environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(diflink) <- asNamespace("BTSR")
  structure(list(linkfun = linkfun,
                 linkinv = linkinv,
                 diflink = diflink,
                 mu.eta = mu.eta,
                 name = link), class = "link-btsr")



}


##-------------------------------------------------------------------------
## internal function. Converts the link to the corresponding integer
## to be passed to FORTRAN
##-------------------------------------------------------------------------
.link.convert <- function(link){

  ##------------------------------------------
  ##          Links
  ##------------------------------------------
  ##  0 = linear: f(x) = ax, a real
  ##  1 = logit:  f(x) = log(x/(1-x))
  ##  2 = log:    f(x) = log(x)
  ##  3 = loglog: f(x) = log(-log(x))
  ##  4 = cloglog: f(x) = log(-log(1-x))
  ##------------------------------------------
  links <- matrix(0:4, ncol = 1)
  nm <- c("linear", "logit",  "log", "loglog",  "cloglog")
  rownames(links) <- nm

  lk <- numeric(length(link))
  for(i in 1:length(lk)){
    lk[i] <- ifelse(link[i] %in% rownames(links), links[link[i],1], NA)
    if(is.na(lk[i])){
      mes <- paste(link[i], "link not available, available links are: ")
      for(j in 1:nrow(links)) mes <- paste(mes, "'", nm[j], "', ", sep = "")
      stop(mes)
    }
  }
  as.integer(lk)
}



##-------------------------------------------------------------------------
## internal function. Checks compatibility
##-------------------------------------------------------------------------
.link.check <- function(model, link){

  lk = unique(link)
  ok = TRUE
  for(i in 1:length(lk)){
    if(model == "GARFIMA"){
      # g1 and g2 cannot be from (a,b) -> (-Inf, Inf) because
      # the data does not have an upper bound
      if(lk[i] %in% c("logit","loglog", "cloglog")) ok = FALSE
    }else{
      # g2 can be anything, g1 must be from (a,b) -> (-Inf, Inf)
      if(i == 1 & lk[i] %in% c("log")) ok = FALSE
    }
  }

  if(!ok) stop("The selected model and link are not compatible")

  invisible(ok)

}
