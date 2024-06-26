#' Calculate Root Mean Square Error (RMSE)
#'
#' Function to calculate RMSE error
#'
#' @param x Residuals of forecast.
#' @param range.term Number of forecast step.
#' @param range.m Number of \code{ts} in each hierarchy.
#' @param weight.m Weights matrix for different series with different ahead-forecasting steps.
#'
#' @return RMSE error of RoME forecast.
#' @export
#'

RMSE <- function(x, range.term, range.m, weight.m = NULL) {

  x <- rbind(x); size.m <- ncol(x); length.term <- length(range.term); length.m <- length(range.m)

  if (is.null(weight.m))
    weight.m <- rep(1, size.m)

  x.squared <- x^2
  x.weighted <- t(t(x.squared) * weight.m)

  x.mse <- matrix(nrow = length.term, ncol = length.m)
  for (i in 1:length.term)
    for (j in 1:length.m)
      x.mse[i, j] <- mean(x.squared[range.term[[i]], range.m[[j]]])

  x.rmse <- sqrt(x.mse)

  return(x.rmse)
}
