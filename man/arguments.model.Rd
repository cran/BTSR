% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help.R
\name{arguments.model}
\alias{arguments.model}
\title{Available models in BTSR package}
\arguments{
\item{model}{character string (case-insensitive) indicating the model to be
fitted to the data. Must be one of the options listed in the Section
\link[=arguments.model]{Supported Models}.}
}
\description{
The BTSR package supports a variety of models, including
\itemize{
\item i.i.d structure,
\item regression models,
\item short- and long-memory time series models
\item chaotic processes.
}

This documentation describes
\itemize{
\item the \code{model} argument and available model strings,
\item default configurations for specific models,
\item how to reproduce models from literature.
}
}
\section{Supported Models}{
Internally, all models are handled by the same function and all models can be
obtained from the more general case \code{"*ARFIMAV"}. When a particular model
(e.g. \code{"BREG"} or \code{"BARMA"}) is invoked some default values are assumed.

The following table summarizes the available distributions and the
corresponding string to generate each model type. The character \code{V} at
the end of the string indicates that \eqn{\nu} is time-varying.

\if{html}{\out{<div class="sourceCode">}}\preformatted{+--------------+--------+------------+---------+-----------+---------+
| Distribution | i.i.d. | Regression | Short   | Long      | Chaotic |
|              | sample |            | Memory  | Memory    |         |
+--------------+--------+------------+---------+-----------+---------+
| Beta         | BETA   | BREG       | BARMA   | BARFIMA   | BARC    |
|              |        | BREGV      | BARMAV  | BARFIMAV  |         |
+--------------+--------+------------+---------+-----------+---------+
| Gamma        | GAMMA  | GREG       | GARMA   | GARFIMA   |         |
|              |        | GREGV      | GARMAV  | GARFIMAV  |         |
+--------------+--------+------------+---------+-----------+---------+
| Kumaraswamy  | KUMA   | KREG       | KARMA   | KARFIMA   |         |
|              |        | KREGV      | KARMAV  | KARFIMAV  |         |
+--------------+--------+------------+---------+-----------+---------+
| Matsuoka     | MATSU  | MREG       | MARMA   | MARFIMA   |         |
+--------------+--------+------------+---------+-----------+---------+
| Unit-Lindley | UL     | ULREG      | ULARMA  | ULARFIMA  |         |
+--------------+--------+------------+---------+-----------+---------+
| Unit-Weibull | UW     | UWREG      | UWARMA  | UWARFIMA  |         |
|              |        | UWREGV     | UWARMAV | UWARFIMAV |         |
+--------------+--------+------------+---------+-----------+---------+
}\if{html}{\out{</div>}}
}

\section{Default values}{
All models are special cases of the general \code{"*ARFIMAV"} structure. When a
specific model is selected via \code{model = "NAME"}, the package automatically
applies these default configurations (any parameter that does not appear in
the equations below is ignored)

\strong{i.i.d samples} (e.g., BETA, GAMMA,...)
\deqn{
 \eta_{1t} = \alpha_1 = \mu, \quad \eta_{2t} = \alpha_2  = \nu.
}
Fixed

\if{html}{\out{<div class="sourceCode">}}\preformatted{p <- q <- d <- 0
xreg <- NULL
linkg <- list(g11 = "linear", g2 = "linear",
              g21 = "linear", g23 = "linear")
}\if{html}{\out{</div>}}

\strong{Regression models with \eqn{\nu_t} constant over time} (e.g., BREG,
GREG,...)
\deqn{
 \eta_{1t} =
 g_{11}(\mu_t) = \alpha_1 + \boldsymbol{X}_{1t}'\boldsymbol{\beta}_1,
 \quad
 \eta_{2t} = \alpha_2 = \nu.
 }
Fixed

\if{html}{\out{<div class="sourceCode">}}\preformatted{p <- q <- d <- 0
xreg <- list(part1 = "user's regressors", part2 = NULL)
linkg <- list(g11 = "user's choice", g12 = "linear",
              g2 = "linear", g21 = "linear", g23 = "linear")
}\if{html}{\out{</div>}}

\strong{Regression models with \eqn{\nu_t} varying on time} (e.g. BREGV,
GREGV)
\deqn{
  \eta_{1t}  =
  g_{11}(\mu_t) = \alpha_1 + \boldsymbol{X}_{1t}'\boldsymbol{\beta}_1,
  \quad
  \eta_{2t}  =
  g_{21}(g_2(\nu_t)) = \alpha_2 + \boldsymbol{X}_{2t}'\boldsymbol{\beta}_2.
 }
Fixed

\if{html}{\out{<div class="sourceCode">}}\preformatted{p <- q <- d <- 0
linkg <- list(g11 = "user's choice", g12 = "linear",
              g2 = "user's choice", g21 = "user's choice",
              g22 = "linear", g23 = "linear")
}\if{html}{\out{</div>}}

\strong{Short-memory models with \eqn{\nu} constant over time (ARMA-like)} (e.g. BARMA,
GARMA,...)
\deqn{
  \begin{aligned}
  \eta_{1t} & =
  g_{11}(\mu_t) = \alpha_1 + \boldsymbol{X}_{1t}'\boldsymbol{\beta}_1 +
  \sum_{i=1}^{p_1} \phi_{1i}\bigl(g_{12}(Y_{t-i})-
  I_{X_1}\boldsymbol{X}_{1(t-i)}'\boldsymbol{\beta}_1\bigr) +
  \sum_{k=1}^{q_1} \theta_{1k} e_{1,t-k}, \\
  \eta_{2t} &  = \alpha_2 = \nu.
  \end{aligned}
  }
Fixed

\if{html}{\out{<div class="sourceCode">}}\preformatted{d <- 0
xreg <- list(part1 = "user's regressors", part2 = NULL)
linkg <- list(g11 = "user's choice", g12 = "user's choice",
              g2 = "linear", g21 = "linear", g23 = "linear")
}\if{html}{\out{</div>}}

\strong{Short-memory models with \eqn{\nu_t} varying on time} (e.g. BARMAV,
GARMAV,...)
\deqn{
  \begin{aligned}
  \eta_{1t} & =
  g_{11}(\mu_t) =\alpha_1 + \boldsymbol{X}_{1t}'\boldsymbol{\beta}_1 +
  \sum_{i=1}^{p_1} \phi_{1i}\big(g_{12}(Y_{t-i})-
  I_{X_1}\boldsymbol{X}_{1(t-i)}'\boldsymbol{\beta}_1\big) +
  \sum_{k=1}^{q_1} \theta_{1k} r_{t-k},\\
  \vartheta_t & = g2(\nu_t)\\
  \eta_{2t} &
  = g_{21}(\vartheta_t) =\alpha_2 + \boldsymbol{X}_{2t}' \boldsymbol{\beta}_2
  + \sum_{i=1}^{p_2} \phi_{2i}\big(g_{22}(\vartheta_{t-i})-
  I_{X_2}\boldsymbol{X}_{2(t-i)}'\boldsymbol{\beta}_2\big) +
  \sum_{k=1}^{q_2} \theta_{2k}  g_{23}(e_{1,t-k}).
  \end{aligned}
 }
Fixed

\if{html}{\out{<div class="sourceCode">}}\preformatted{d <- 0
}\if{html}{\out{</div>}}

\strong{Long-memory models with \eqn{\nu} constant over time (ARFIMA-like
models)} (e.g. BARFIMA, GARFIMA,...)
\deqn{
 \begin{aligned}
 \eta_{1t} & =
 g_{11}(\mu_t) =\alpha_1 + \boldsymbol{X}_{1t}'\boldsymbol{\beta}_1 +
 \sum_{i=1}^{p_1} \phi_{1i}\big(g_{12}(Y_{t-i})-
 I_{X_1}\boldsymbol{X}_{1(t-i)}'\boldsymbol{\beta}_1\big) +
 \sum_{k=1}^\infty c_{1k} r_{t-k},\\
 \eta_{2t} &  =\alpha_2 =  \nu.
 \end{aligned}
}
Fixed

\if{html}{\out{<div class="sourceCode">}}\preformatted{p <- c("user's p", 0)
q <- c("user's q", 0)
d <- c("user's d", 0)
xreg <- list(part1 = "user's regressors", part2 = NULL)
linkg <- list(g11 = "user's choice", g12 = "user's choice",
              g2 = "linear", g21 = "linear", g23 = "linear")
}\if{html}{\out{</div>}}
}

\section{Reproducing Models from the Literature}{
This section summarizes how to replicate well-known time series models from
the literature using the BTSR package. For each model type, we provide the
necessary parameter settings and references to the original publications.
These configurations act as templates, helping users correctly apply the
package to reproduce results or extend established models.

Key arguments (e.g., \code{error.scale}, \code{xregar}, \code{y.lower},
\code{y.upper}, \code{rho}) should be set to match the specifications in the
referenced articles. While we focus on the \verb{btsr.*} functions (see
\link{BTSR.functions}), all models can also be implemented using the corresponding
parent model functions (for details, see \link{BTSR.parent.models}).

\strong{i.i.d. samples}: The arguments \code{error.scale} and \code{xregar} are
ignored.
\itemize{
\item Beta distribution with parameters  \code{shape1} and \code{shape2} compatible
with the one from \link{rbeta}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "BETA"
alpha = shape1/(shape1 + shape2)
nu = shape1 + shape2
}\if{html}{\out{</div>}}

\item Gamma distribution with parameters \code{shape} and \code{scale} compatible with
the one from \link{rgamma}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "GAMMA"
alpha = shape*scale
nu = shape
}\if{html}{\out{</div>}}

\item Kumaraswamy distribution with shape parameters \code{shape1} and \code{shape2}
(respectively denoted by \eqn{a} and \eqn{b} in
\insertCite{kumaraswamy1980;nobrackets}{BTSR}):

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "KUMA"
alpha = (y.lower - y.upper)*(1 - (1-rho)^1/shape2)*1/shape1 + y.lower
nu = shape1
}\if{html}{\out{</div>}}

\strong{Warning:} Choose \eqn{\mu},  \eqn{\nu} and \eqn{\rho} carefully
since \eqn{|\log(1-\rho)| >> |\log(1 - \mu^\nu)|} may cause numerical
instability.

\item Matsuoka distribution with shape parameter \code{shape}
\insertCite{matsuoka2024}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "MATSU"
alpha = (shape/(shape+1))^(3/2)
}\if{html}{\out{</div>}}

\item Unit-Lindley distribution with parameter \code{theta}
\insertCite{mazucheli2018}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "UL"
alpha = 1/(1 + theta)
}\if{html}{\out{</div>}}

\item Unit-Weibull distribution with parameter \code{mu}, \code{beta} and \code{tau} from
\insertCite{mazucheli2019}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "UW"
alpha = mu
nu = beta
rho = tau
}\if{html}{\out{</div>}}

}

\strong{Regression models}: the argument \code{error.scale} and all entries but
\code{g11} in \code{linkg} are ignored
\itemize{
\item Beta regression \insertCite{ferrari2004}{BTSR}: \code{model = "BREG"}

\item Kumaraswamy regression \insertCite{mitnikbaek2013}{BTSR}:
\code{model = "KREG"}.

\item Unit-Lindley regression \insertCite{mazucheli2018}{BTSR}:
\code{model = "ULREG"}.

\item Unit-Weibull regression \insertCite{mazucheli2019}{BTSR}:
\code{model = "UWREG"}.
}

\strong{ARMA-like models}
\itemize{
\item BARMA model \insertCite{rocha2009,rocha2017}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "BARMA"
error.scale = 1
xregar = TRUE
}\if{html}{\out{</div>}}

\item KARMA model \insertCite{bayer2017}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "KARMA"
error.scale = 1
xregar = TRUE
y.lower = 0
y.upper = 1
rho = 0.5
}\if{html}{\out{</div>}}

\item GARMA model \insertCite{prass2025}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "GARMA"
error.scale = 0
}\if{html}{\out{</div>}}

\item MARMA model \insertCite{pumi2025matsu}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "MARMA"
error.scale = 1
xregar = TRUE
}\if{html}{\out{</div>}}

\item ULARMA model \insertCite{pumi2025ul}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "ULARMA"
error.scale = 1
xregar = TRUE
}\if{html}{\out{</div>}}

}

\strong{ARFIMA-like models}
\itemize{
\item BARFIMA model \insertCite{pumi2019}{BTSR}:

\if{html}{\out{<div class="sourceCode">}}\preformatted{model = "BARFIMA"
error.scale = 1
xregar = TRUE
d = TRUE (for fitting)
}\if{html}{\out{</div>}}

}

\strong{Chaotic models}
\itemize{
\item BARC model \insertCite{pumi2021}{BTSR}: set \code{model = "BARC"} and
\code{error.scale = 1}.
}
}

\references{
\insertAllCited{}
}
\seealso{
\link{BTSR.models}, \link{BTSR.model.defaults}, \link{get.defaults}
}
