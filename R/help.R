#' @title Table of available model
#'
#' @description
#' This function returns the table of available models.
#'
#' @param do.plot logical; if `TRUE` returns a plot with the output, otherwise
#'   prints the results in the console.
#'
#' @return
#' `NULL` (invisibly). The function is called for its side effects
#' (printing/plotting).
#'
#' @export
BTSR.models <- function(do.plot = interactive()) {
  .BTSR.model.table(do.plot)
}


#' @title
#' Available models in BTSR package
#'
#' @name arguments.model
#'
#' @description
#' The BTSR package supports a variety of models, including
#' \itemize{
#'   \item i.i.d structure,
#'   \item regression models,
#'   \item short- and long-memory time series models
#'   \item chaotic processes.
#' }
#'
#' This documentation describes
#' \itemize{
#'  \item the `model` argument and available model strings,
#'  \item default configurations for specific models,
#'  \item how to reproduce models from literature.
#' }
#'
#' @template param_models
#'
#' @seealso [BTSR.models], [BTSR.model.defaults], [get.defaults]
#'
NULL
#> NULL

#' @title
#' Shared documentation for the time series
#'
#' @name arguments.series
#'
#' @description
#' This is the common documentation for arguments related to the
#' observed/simulated time series and its conditional distribution.
#'
#' @template param_series
NULL
#> NULL

#' @title
#' Shared documentation for regressors
#'
#' @name arguments.regressors
#'
#' @description
#' This is the common documentation for arguments related to the regressors.
#'
#' @template param_regressors
NULL
#> NULL

#' @title
#' Shared documentation for models order
#'
#' @name arguments.order
#'
#' @description
#' This is the common documentation for arguments related to order of
#' polynomials and truncation points for infinite sums, presented in BTSR models.
#'
#' @template param_order
NULL
#> NULL
#>

#' @title
#' Shared documentation for coefficients
#'
#' @name arguments.coefs
#'
#' @description
#' This is the common documentation for arguments related to the coefficients in
#' BTSR models.
#'
#' @template param_coefficients
NULL
#> NULL

#' @title
#' Available map functions in BTSR package
#'
#' @name arguments.map
#'
#' @description
#' This documentation describes the `map` argument in [BARC][BARC.functions]
#' models and the map functions implemented in the BTSR package.
#'
#' @template param_map
NULL
#> NULL

#' @title
#' Shared documentation for link functions
#'
#' @name arguments.link
#'
#' @description
#' This is the common documentation for arguments related link functions in BTSR
#' models.
#'
#' @template param_link
NULL
#> NULL

#' @title
#' Shared documentation for log-likelihood
#'
#' @name arguments.loglik
#'
#' @description
#' This is the common documentation for arguments related the log-likelihood
#' functions, score vector and information matrix for BTSR models.
#'
#' @template param_loglik
NULL
#> NULL

#' @title
#' Shared documentation for configuration related parameteres
#'
#' @name arguments.configs
#'
#' @description
#' This is the common documentation for arguments related the configurations
#' for fitting models and printing reports.
#'
#' @template param_configs
NULL
#> NULL
