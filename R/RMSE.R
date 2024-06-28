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
#' data("Tourism")
#' data <- Tourism$data
#' size.n = ncol(data)
#' H1 <- Tourism$H1
#' diag1<-diag(1,size.n)
#' S <- as.matrix(rbind(H1,diag1))
#' weight = 1 / rowSums(S)
#' data.all = as.matrix(data) %*% t(S)
#' size.T = nrow(data)
#' size.m = ncol(data)

#' order.train = 1:96
#' order.test = (1:12)+96

#' data.train = data.all[order.train, ]
#' data.test = data.all[order.test, ]

#' sol.base = BaseAuto(x = data.train, method.forecast = 'ETS', term = 12)
#' sol.Wt = CovarianceW(data.train, method.forecast = 'ETS', is.rolling = FALSE)
#' sol.LAD <- RoME_fore(base = sol.base$base, hierarchy = S, loss.diff = DiffLoss(type = 'LAD'),
#' design.W = 'OLS', matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
#' residuals = sol.LAD$RoME - data.test
#' rmse <- RMSE(residuals, range.term = list(1:12),
#' range.m = list(1,2:8,9:35,36:111), weight.m = weight)
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
