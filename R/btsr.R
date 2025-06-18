#' @title
#' \packageTitle{BTSR}
#'
#' @name btsr-package
#'
#' @aliases BTSR BTSR-Package
#'
#' @description
#' The BTSR package provides a unified framework for simulating, fitting, and
#' forecasting bounded time series regression models. It supports a wide
#' range of models, including i.i.d., regression, ARMA-like, and ARFIMA-like
#' models, with a focus on bounded time series data.
#'
#' Key features of the BTSR package include
#' \itemize{
#'  \item Simulation of bounded time series data using various models.
#'  \item Estimation of model parameters using efficient algorithms.
#'  \item Forecasting future values based on fitted models.
#'  \item Support for both short-memory and long-memory models.
#'  \item Flexible link functions and error scales.
#' }

#' @template section_mathematical_notation
#' @template section_btsr_structure
#'
#' @examples
#' #----------------------------
#' # Quickstart examples.
#' #----------------------------
#'
#' # Example 1: Simulate i.i.d. samples
#' set.seed(1234)
#' y1 <- btsr.sim(model = "BETA", n = 1000, coefs = list(alpha = 0.2, nu = 20))
#' hist(y1)
#'
#' # Example 2: Simulate ARMA-like model with fixed nu
#' y2 <- btsr.sim(
#'   model = "BARMA", n = 100, link = "logit",
#'   coefs = list(alpha = 0.2, phi = 0.5, theta = 0.3, nu = 20)
#' )
#' plot(y2, type = "l")
#'
#' @seealso
#' For detailed examples and usage instructions, see the documentation for
#' individual functions
#' \itemize{
#'  \item [btsr.sim]: Simulate bounded time series data.
#'
#'  \item [btsr.extract]: Extract components of a BTSR model, for a given set of
#'   parameters
#'
#'  \item [btsr.fit]: Fit a BTSR model to data.
#'
#'  \item [predict]: Forecast future values using a fitted model.
#'
#'  \item [arguments][arguments]: Shared documentation for arguments
#' }
#'
#' @author
#' Taiane Schaedler Prass \email{taianeprass@@gmail.com},
#' Guilherme Pumi \email{guipumi@@gmail.com}
#'
#' @importFrom Rdpack reprompt
#
#' @references
#'  \insertRef{bayer2017}{BTSR}
#'
#'  \insertRef{pumi2019}{BTSR}
#'
#'  \insertRef{pumi2021}{BTSR}
#'
#'  \insertRef{pumi2024uw}{BTSR}
#'
#'  \insertRef{pumi2024uwc}{BTSR}
#'
#'  \insertRef{pumi2025ul}{BTSR}
#'
#'  \insertRef{pumi2025matsu}{BTSR}
#'
#'  \insertRef{prass2025}{BTSR}
#'
#' @useDynLib BTSR, .registration=TRUE
#'
#' @keywords package distribution regression
"_PACKAGE"
