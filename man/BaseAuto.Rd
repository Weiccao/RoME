% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_fore.R
\name{BaseAuto}
\alias{BaseAuto}
\title{Forecast HTS with ARIMA or ETS method}
\usage{
BaseAuto(x, method.forecast = "ARIMA", term = 1)
}
\arguments{
\item{x}{a \code{hts} object}

\item{method.forecast}{Forecast model. For more details, see \code{\link{auto.arima}} and \code{\link{ets}}.}

\item{term}{Number of ahead forecast steps.}
}
\value{
a list object contains: \item{base}{h-step forecasts results}
\item{residual}{in sample forecast error}
\item{x}{the original time series}
}
\description{
Return base forecast result of \code{arima} or \code{ets} method and forecast residuals.
}
