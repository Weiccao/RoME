#' Forecast with ARIMA or ETS method
#'
#' Return forecast result of \code{ARIMA} or \code{ETS} model
#'
#' @param x a univariate time series.
#' @param method.forecast Forecast model. For more details, see \code{\link{auto.arima}} and \code{\link{ets}}.
#' @param term Number of ahead forecast steps.
#'
#' @return Original time series and the h-step forecasts.
#' @export
#'
#' @examples
#' # Arima
#' ForecastAuto(WWWusage, method.forecast = 'ARIMA', term = 6)
#'
#' # ETS
#' ForecastAuto(WWWusage, method.forecast = 'ETS', term = 6)
#'
ForecastAuto <- function(x, method.forecast = 'ARIMA', term = 0) {

  if (method.forecast == 'ARIMA')
    sol.forecast <- auto.arima(y = x)

  if (method.forecast == 'ETS')
    sol.forecast <- ets(y = x)

  x.forecast <- as.numeric(fitted(sol.forecast))

  if (term) {
    x.forward <- forecast(sol.forecast, h = term)
    x.forecast <- c(x.forecast, as.numeric(x.forward$mean))
  }

  return(x.forecast)
}

#' Forecast HTS with ARIMA or ETS method
#'
#' Return base forecast result of \code{arima} or \code{ets} method and forecast residuals.
#'
#' @param x a \code{hts} object
#' @param method.forecast Forecast model. For more details, see \code{\link{auto.arima}} and \code{\link{ets}}.
#' @param term Number of ahead forecast steps.
#'
#' @return a list object contains: \item{base}{h-step forecasts results}
#' \item{residual}{in sample forecast error}
#' \item{x}{the original time series}
#'
#' @export
#'
#' @examples
#' Generate_NonGaussianSeries(type.epsilon = 'Mixture',times.repeat = 1)
#' data <- read.table('NonGaussianSeries_epsilon_Mixture_y.txt')
#' hierarchy <- rbind(rep(1, 9), c(rep(1, 3), rep(0, 6)), c(rep(0, 3), rep(1, 3), rep(0, 3)), c(rep(0, 6), rep(1, 3)), diag(1, 9))
#' data.train <- as.matrix(data[1:180, ]) %*% t(hierarchy)
#' data.test <- as.matrix(data[-c(1:180), ]) %*% t(hierarchy)
#' sol.base <- BaseAuto(data.all, term = 12)
#'

BaseAuto <- function(x, method.forecast = 'ARIMA', term = 1) {

  size.T <- nrow(x); size.m <- ncol(x)

  x.all <- matrix(nrow = size.T + term, ncol = size.m)
  for (i in 1:size.m)
    x.all[, i] <- ForecastAuto(x = x[, i], method.forecast = method.forecast, term = term)

  x.forecast <- NULL
  if (term)
    x.forecast <- rbind(x.all[-(1:size.T), ])

  x.residual <- x.all[1:size.T, ] - x

  list.return <- list(base = x.forecast,
                      residual = x.residual,
                      x = x)
  return(list.return)
}

#' Bottom-Up HTS forecast model
#'
#' Return forecast result of bottom-up method.
#'
#' @param base base forecast of \code{hts}
#' @param hierarchy summing matrix
#'
#' @return Forecast results of bottom-up method.
#' @export
#' @examples
#'
#' # generate hts
#' Generate_NonGaussianSeries(type.epsilon = 'Mixture')
#' data <- read.table('NonGaussianSeries_epsilon_Mixture_y.txt')
#' hierarchy <- rbind(rep(1,9), c(rep(1,3), rep(0,6)), c(rep(0,3), rep(1,3), rep(0,3)), c(rep(0,6), rep(1,3)), diag(1,9))
#' data.all <- as.matrix(data[1:180, ]) %*% t(hierarchy)
#'
#' # forecast
#' sol.base <- BaseAuto(data.all, term = 12)
#' sol.bu <- BottomUp(sol.base, hierarchy)


BottomUp <- function(base, hierarchy) {

  size.mstar <- nrow(hierarchy) - ncol(hierarchy)

  sol.bu <- base[, -(1:size.mstar)] %*% t(hierarchy)

  return(sol.bu)
}
