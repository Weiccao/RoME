% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_hts.R
\name{Generation_VariationOfNormality}
\alias{Generation_VariationOfNormality}
\title{Function to generate \code{hts} with different variation parameter}
\usage{
Generation_VariationOfNormality(
  size.n = 6,
  size.T = 192,
  alpha = 0.6,
  rho = 0.4,
  sigma
)
}
\arguments{
\item{size.n}{The total number of \code{ts}.}

\item{size.T}{The number of observations of each \code{ts}.}

\item{alpha}{Value of the autoregressive coefficient.}

\item{rho}{The correlation coefficient parameter between any two standard normality \code{ts}.}

\item{sigma}{The value of variation parameter.}
}
\value{
\code{hts} with different variation parameter.
}
\description{
Function to generate \code{hts} with different variation parameter
}
