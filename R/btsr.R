##' @title
##' BTSR: Bounded Time Series Regression.
##'
##' @description
##' The BTSR package provides functions to simulate, estimate and forecast a
##' wide range of regression based dynamic models for bounded time series. The
##' package covers the most commonly applied models in the literature.
##' The package's main calculations are done in FORTRAN, which translates into
##' very fast algorithms.
##'
##' @author Taiane Schaedler Prass \email{taianeprass@@gmail.com}
##' @docType package
##' @name BTSR.Package
##' @aliases BTSR
##' @useDynLib BTSR, .registration=TRUE
##'
##' @section The BTSR structure:
##'
##' The general structure of the deterministic part of a BTSR model is
##'
##'  \deqn{g_1(\mu_t) = \alpha + X_t\beta +
##'  \sum_{j=1}^p \phi_j[g_2(y_{t-j}) - I_{xregar}X_{t-j}\beta] + h_t}
##'
##' where
##'  \itemize{
##'    \item \eqn{I_{xregar}} is 0, if \code{xreg} is not included in the AR part of the model and 1,
##' otherwise
##'
##'   \item the term \eqn{h_t} depends on the argument \code{model}:
##'    \itemize{
##'     \item for BARC models: \eqn{h_t =  h(T^{t-1}(u_0))}
##'     \item otherwise: \eqn{h_t =  \sum_{k = 1}^\infty c_k r_{t-k}}
##'    }
##'
##'   \item \eqn{g_1} and \eqn{g_2} are the links defined in \code{linkg}.
##'   Notice that \eqn{g_2} is only used in the AR part of the model and, typically,
##'   \eqn{g_1 = g_2}.
##'
##'   \item \eqn{r_t} depends on the \code{error.scale} adopted:
##'   \itemize{
##'     \item  if \code{error.scale = 0}: \eqn{r_t = y_t - \mu_t} (data scale)
##'     \item if \code{error.scale = 1}:  \eqn{r_t = g_1(y_t) - g_1(\mu_t)}
##'      (predictive scale)
##'   }
##'
##'    \item \eqn{c_k} are the coefficients of \eqn{(1-L)^d\theta(L)}.
##'    In particular, if \eqn{d = 0}, then \eqn{c_k = \theta_k}, for
##'    \eqn{k = 1, \dots, q}, and 0 otherwise.
##' }
##'
NULL
##> NULL
