% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_misc.R
\name{fit.control}
\alias{fit.control}
\title{Default control list}
\usage{
fit.control(control = list())
}
\arguments{
\item{control}{a list with configurations to be passed to the optimization
subroutines. Missing arguments will receive default values. See
\sQuote{Details}.}
}
\value{
Returns a list with all arguments in \sQuote{Details}.
}
\description{
Sets default values for constants used by the optimization functions in
FORTRAN.
}
\details{
The \code{control} argument is a list that can supply any of the following
components
\describe{
\item{\code{method}}{The optimization method. Current available options are
\code{"L-BFGS-B"} and \code{"Nelder-Mead"}. Default is \code{"L-BFGS-B"}.}

\item{\code{maxit}}{The maximum number of iterations. Defaults is \code{1000}.}

\item{\code{iprint}}{The frequency of reports if \code{control$trace} is positive.
Defaults is -1 (no report).
\itemize{
\item For \code{"L-BFGS-B"} method:
\itemize{
\item \code{iprint < 0}  no output is generated;
\item \code{iprint = 0}  print only one line at the last iteration;
\item \verb{0 < iprint < 99} print also f and |proj g| every iprint iterations;
\item \code{iprint = 99} print details of every iteration except n-vectors;
\item \code{iprint = 100}  print also the changes of active set and final x;
\item \code{iprint > 100}  print details of every iteration including x and g;
}

\item For \code{"Nelder-Mead"} method:
\itemize{
\item \code{iprint < 0} No printing
\item \code{iprint = 0} Printing of parameter values and the function value
after initial evidence of convergence.
\item \code{iprint > 0} As for \code{iprint = 0} plus progress reports after every
\code{iprint} evaluations, plus printing for the initial simplex.
}
}}

\item{\code{factr}}{controls the convergence of the \code{"L-BFGS-B"}  method.
Convergence occurs when the reduction in the objective is within this
factor of the machine tolerance. The iteration will stop when
\deqn{\dfrac{(f^k - f^{k+1})}{\max\{|f^k|,|f^{k+1}|,1\}} \le
  factr \times epsmch}
where \eqn{epsmch} is the machine precision, which is automatically
generated by the code. Typical values for \code{factr}: 1.e+12 for low accuracy;
1.e+7 for moderate accuracy; 1.e+1 for extremely high accuracy. Default is
\code{1e7}, that is a  tolerance of about \code{1e-8}.}

\item{\code{pgtol}}{helps control the convergence of the \code{"L-BFGS-B"} method.
It is a tolerance on the projected gradient in the current search
direction. the iteration will stop when
\deqn{\max\{|\text{proj }g_i |, i = 1, ..., n\} \le pgtol}
where \eqn{\text{proj }g_i} is the ith component of the projected gradient.
Default is \code{1e-12}.}

\item{\code{stopcr}}{The criterion applied to the standard deviation of the
values of objective function at the points of the simplex, for
\code{"Nelder-Mead"} method.}
}
}
\examples{
BTSR::fit.control()

}
