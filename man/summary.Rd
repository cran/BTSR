% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_misc.R
\name{summary}
\alias{summary}
\alias{summary.btsr}
\alias{print.summary.btsr}
\title{Summary Method of class BTSR}
\usage{
\method{summary}{btsr}(object, robust = FALSE, outer = FALSE,
  full.report = TRUE, ...)

\method{print}{summary.btsr}(x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
}
\arguments{
\item{object}{object of class \code{"btsr"}.}

\item{robust}{logical. If \code{TRUE} the robust covariance matrix is
computed}

\item{outer}{logical. If \code{TRUE} uses the outer product of the gradient
to compute the covariance matrix. If \code{robust = TRUE}, \code{outer} is
used as a second option (in case of error computing the robust version)}

\item{full.report}{logical. If \code{TRUE} prints a more detailed report.}

\item{...}{further arguments passed to or from other methods.}

\item{x}{an object of class \code{"summary.btsr"}, usually, a result of a call to
\code{summary.btsr}.}

\item{digits}{minimal number of significant digits, see \link{print.default}.}

\item{signif.stars}{logical. If \code{TRUE}, \sQuote{significance stars} are
printed for each coefficient.}
}
\value{
The function \code{summary.btsr} computes and returns a list of summary statistics
of the fitted model given in \code{object}. Returns a list of class
\code{summary.btsr}, which contains the following components

\itemize{
\item \code{model}: the corresponding model.

\item \code{call}: the matched call.

\item \code{residuals}: The (in-sample) residuals, that is, the observed values
\eqn{Y_t} minus the fitted values \eqn{\mu_t}. The same as the \code{error} term
if \code{error.scale = 0}.

\item \code{coefficients}: a \eqn{k \times 4}{k x 4} matrix with columns for the
estimated coefficient, its standard error, z-statistic and corresponding
(two-sided) p-value.

\item \code{sigma.res}: the square root of the estimated variance of the error
term in \code{residuals}
\deqn{
 \hat\sigma^2 = \displaystyle\frac{1}{n-k}\sum_{i=1}^{n-k}{e_i^2},
 }
where \eqn{e_i} is the \eqn{i}-th residual.

\item \code{df}: degrees of freedom, a 2-vector \eqn{(k, n-k)}, the first being
the number of estimated coefficients.

\item \code{vcov}: a \eqn{k \times k}{k \times k} matrix of (unscaled)
covariances. The inverse ov the information matrix.

\item \code{loglik}: the sum of the log-likelihood values (\eqn{L})

\item \code{aic}: the AIC value. \eqn{AIC = -2L + 2k}.

\item \code{bic}: the BIC value. \eqn{BIC = -2L + k\log(n)}.

\item \code{hqc}: the HQC value. \eqn{HQC = -2L + k\log(\log(n))}.
}
}
\description{
\code{summary} method for class \code{"btsr"}.
}
\details{
\code{print.summary.btsr} tries to be smart about formatting the coefficients,
standard errors, etc. and additionally provides \sQuote{significance stars}.
}
