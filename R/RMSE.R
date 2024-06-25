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
#' @examples
#' Generate_NonGaussianSeries(type.epsilon = 'Mixture',times.repeat = 1)
#' data <- read.table('NonGaussianSeries_epsilon_Mixture_y.txt')
#' hierarchy <- rbind(rep(1, 9), c(rep(1, 3), rep(0, 6)), c(rep(0, 3), rep(1, 3), rep(0, 3)), c(rep(0, 6), rep(1, 3)), diag(1, 9))
#' data.train <- as.matrix(data[1:180, ]) %*% t(hierarchy)
#' data.test <- as.matrix(data[-c(1:180), ]) %*% t(hierarchy)
#' sol.base <- BaseAuto(data.all, term = 12)
#' W <- t(sol.base$residual) %*% sol.base$residual / 180
#' lambda <- LambdaShrink(sol.base$residual)
#'
#' # forecast with RoME
#' sol.RoME <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'))
#' res <- sol.RoME$RoME - data.test
#' rmse <- RMSE(res, range.term = list(1:12),range.m = list(1,2:4,5:13))
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
