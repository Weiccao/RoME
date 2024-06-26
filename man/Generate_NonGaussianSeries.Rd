% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_hts.R
\name{Generate_NonGaussianSeries}
\alias{Generate_NonGaussianSeries}
\title{Function to generate \code{hts} with non-Gaussian residuals.}
\usage{
Generate_NonGaussianSeries(
  size.n = 9,
  size.c = 3,
  size.T = 192,
  alpha = 0.8,
  rho = 0.4,
  type.epsilon
)
}
\arguments{
\item{size.n}{The total number of \code{ts}.}

\item{size.c}{The number of \code{ts} with non-Gaussian residuals.}

\item{size.T}{The number of observations of each \code{ts}.}

\item{alpha}{Value of the autoregressive coefficient.}

\item{rho}{The correlation coefficient parameter between any two standard normality \code{ts}.}

\item{type.epsilon}{The distribution of non-Gaussian error term.}
}
\value{
\code{hts} with non-Gaussian residuals
}
\description{
Function to generate \code{hts} with non-Gaussian residuals.
}
