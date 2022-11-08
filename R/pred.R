#' @title Predict method for BTSR
#'
#' @description Predicted values based on btsr object.
#'
#' @param object Object of class inheriting from \code{"btsr"}
#' @param newdata A matrix with new values for the regressors.   If omitted
#' and \code{"xreg"} is present in the model, the fitted values are returned.
#' If the model does not include regressors, the functions will use
#' the value of \code{nnew}.
#' @param nnew number of out-of-sample forecasts required. If \code{newdata} is
#' provided, \code{nnew} is ignored.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#'   \code{predict.btsr} produces predicted values, obtained by evaluating
#'   the regression function in the frame \code{newdata}.
#'
#'   If \code{newdata} is omitted the predictions are based on the data
#'   used for the fit.
#'
#'   For now, prediction intervals are not provided.
#'
#' @return A list with the following arguments
#'
#'  \item{series}{The original time series yt.}
#'
#'  \item{xreg}{The original regressors (if any).}
#'
#'  \item{fitted.values}{The in-sample forecast given by \eqn{\mu_t}.}
#'
#'  \item{etat}{In-sample values of \eqn{g(\mu[t])}.}
#'
#'  \item{error}{The error term (depends on the argument \code{error.scale})}
#'
#'  \item{residuals}{The (in-sample) residuals, that is, the observed minus the predicted values.
#'  Same as error when \code{error.scale} = 0}
#'
#'  \item{forecast}{The predicted values for yt.}
#'
#'  \item{TS}{only for \code{"BARC"} models. The iterated map.}
#'
#'  \item{Ts.forecast}{only for \code{"BARC"} models. The predicted values
#'  of the iterated map.}
#'
#' @examples
##'  #------------------------------------------------------------
##'  # Generating a Beta model were mut does not vary with time
##'  # yt ~ Beta(a,b), a = mu*nu, b = (1-mu)*nu
##'  #------------------------------------------------------------
##'
##' y <- btsr.sim(model= "BARFIMA", linkg = "linear",
##'                n = 100, seed = 2021,
##'                coefs = list(alpha = 0.2, nu = 20))
##'
##' # fitting the model
##' f <- btsr.fit(model = "BARFIMA", yt = y, report = TRUE,
##'              start = list(alpha = 0.5, nu = 10),
##'              linkg = "linear", d = FALSE)
##'
##' pred = predict(f, nnew = 5)
##' pred$forecast
##'
##' @export
##'
predict.btsr <-
  function(object, newdata, nnew = 0,...){

    out <- list()
    nms.out <- c("series", "xreg", "fitted.values", "etat", "error", "residuals", "forecast")
    if(object$model == "BARC") nms.out = c(nms.out, "Ts", "Ts.forecast")

    if(missing(newdata)) newdata = NULL

    if(is.null(newdata) & nnew <= 0){
      #------------------------------------------------------
      # New data was not provided.
      # Extracting existing components and returning
      #------------------------------------------------------
      out[nms.out] <- object[nms.out]
    }else{

      if(is.null(newdata) & object$configs$nreg > 0)
        stop("Please, provide the new values for the regressors")

      #------------------------------------------------------
      # New data was provided.
      # Making the necessary calculations
      #------------------------------------------------------
      xnew = NULL
      if(!is.null(newdata)){
        xnew = as.matrix(newdata)
        nnew = nrow(xnew)
      }

      nms = names(object)
      obj <- object[nms[nms != "configs"]]
      obj[names(object$configs)] <- object$configs
      temp <- .xreg.convert(xreg = object$xreg, xnew = xnew,
                            n = obj$n, nnew = nnew, skip.forecast = FALSE)
      obj[names(temp)] <- temp
      obj[c("llk", "sco", "info", "extra")] <- 0L
      obj$coefs <- obj$coefficients

      if(obj$model == "BARC")
        temp <- .barc.predict(obj, TRUE)
      else
        temp <- .btsr.predict(obj, TRUE)


      out[nms.out] <- obj[nms.out]
      out$forecast <- temp$yt.new
      out$xnew <- NULL
      if(obj$nreg > 0) out$xnew <- xnew
      else out$xreg <- NULL
      if(obj$model == "BARC")  out$Ts.forecast <- temp$Ts.new
    }

    out
  }



##---------------------------------------------------------------------------
## internal function:
## Interface between R and FORTRAN
##---------------------------------------------------------------------------
.btsr.predict <- function(object, debug){

  if(! object$model %in% c("BARFIMA", "GARFIMA", "KARFIMA"))
    stop("The selected model is not implemented yet")

  temp <-  .Fortran("btsrpredictR",
                    n = object$n,
                    series = object$series,
                    ylower = object$y.lower,
                    yupper = object$y.upper,
                    gy = object$gyt,
                    nreg = object$nreg,
                    xreg = object$xreg,
                    escale = object$error.scale,
                    error = object$error,
                    nnew = object$nnew,
                    xnew = object$xnew,
                    ynew = numeric(max(1,object$nnew)),
                    linkg = object$linkg,
                    npar = max(1L,object$npar),
                    coefs = object$coefs,
                    fixa = object$alpha$nfix,
                    alpha = object$alpha$fvalues,
                    fixb = object$beta$nfix,
                    flagsb = object$beta$flags,
                    beta = object$beta$fvalues,
                    p = object$p,
                    fixphi = object$phi$nfix,
                    flagsphi = object$phi$flags,
                    phi = object$phi$fvalues,
                    xregar = object$xregar,
                    q = object$q,
                    fixtheta = object$theta$nfix,
                    flagstheta = object$theta$flags,
                    theta = object$theta$fvalues,
                    fixd = object$d$nfix,
                    d = object$d$fvalues,
                    fixnu = object$nu$nfix,
                    nu = object$nu$fvalues[1],
                    inf = object$inf)

  out <- list(model = object$model,
              yt.new = temp$ynew)
  if(debug) out$out.Fortran <- temp
  invisible(out)
}
