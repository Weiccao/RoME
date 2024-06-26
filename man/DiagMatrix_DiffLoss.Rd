% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RoME_fore.R
\name{DiagMatrix_DiffLoss}
\alias{DiagMatrix_DiffLoss}
\title{Calculate the diagonal matrix for LQA algorithm}
\usage{
DiagMatrix_DiffLoss(residual, loss.diff, varsigma = 1e-08, tau = 1e-16)
}
\arguments{
\item{residual}{The forecast residuals.}

\item{loss.diff}{The type of loss function.}

\item{varsigma}{A positive perturbation, the default is \code{varsigma = 0}.}

\item{tau}{A positive constant, the default is \code{tau = 1e-16}.}
}
\description{
Calculate the diagonal matrix for LQA algorithm
}
