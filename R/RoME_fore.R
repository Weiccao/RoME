#' Calculate the covariance matrix forecast error
#'
#' Return covariance matrix of one step forecast error and shrink parameter.
#'
#' @import MASS
#' @import hts
#' @import forecast
#' @import zoo
#' @import stats
#'
#' @param x a \code{hts} object.
#' @param method.forecast Forecast model. For more details, see \code{\link{auto.arima}} and \code{\link{ets}}
#' @param is.rolling Whether rolling or not.
#' @param rolling.window The window size of rolling.
#'
#' @return a list object contains:
#' \item{W}{covariance matrix}
#' \item{lambda}{shrink parameter}
#' \item{sd}{standard error of residuals}
#'

CovarianceW <- function(x, method.forecast = 'ARIMA', is.rolling = TRUE, rolling.window = NULL) {

  size.T <- nrow(x); size.m <- ncol(x)

  if (is.rolling) {

    if (is.null(rolling.window))
      rolling.window <- size.T - size.m

    times.rolling <- size.T - rolling.window; residuals.x <- matrix(nrow = times.rolling, ncol = size.m)
    for (i in 1:times.rolling) {

      rows.i <- 1:rolling.window + i - 1; x.i <- x[rows.i, ]

      base.i <- BaseAuto(x = x.i, method.forecast = method.forecast, term = 1)
      residuals.x[i, ] <- base.i$base - x[rolling.window + i, ]
    }

    W <- crossprod(residuals.x) / times.rolling
  }

  else {

    sol.base <- BaseAuto(x = x, method.forecast = method.forecast, term = 0)
    residuals.x <- sol.base$residual

    W <- crossprod(as.matrix(residuals.x)) / size.T
  }

  residuals.x <- as.matrix(residuals.x)
  list.return <- list(W = W,
                      lambda = LambdaShrink(residual = residuals.x),
                      sd = sd(residuals.x))
  return(list.return)
}


#' Calculate the sqrt value of covariance matrix of forecast error
#'
#' @param x  The covariance matrix.
#' @param tau A positive constant, the default is \code{tau = 1e-16}.

SqrtMatrix <- function(x, tau = 1e-16) {

  sol.eigen <- eigen(x)

  eigen.min <- min(sol.eigen$values)
  if (eigen.min <= 0)
    sol.eigen$values <- sol.eigen$values - eigen.min + tau

  eigen.sqrt <- sqrt(sol.eigen$values)
  x.sqrt <- sol.eigen$vectors %*% diag(eigen.sqrt) %*% t(sol.eigen$vectors)

  return(x.sqrt)
}

#' Function to calculate shrinkage parameter
#'
#' @param residual The forecast residuls.
#'
LambdaShrink <- function(residual) {
  n <- nrow(residual)
  tar <- diag(apply(residual, 2, crossprod) / n)
  covm <- crossprod(residual) / n; corm <- cov2cor(covm)
  xs <- scale(residual, center = FALSE, scale = sqrt(diag(covm)))
  v <- (crossprod(xs^2) - 1 / n * (crossprod(xs))^2) / (n * (n - 1)); diag(v) <- 0
  corapn <- cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v) / sum(d); lambda <- max(min(lambda, 1), 0)
  return(lambda)
}


#' Define the derivation form of loss function
#'
#' Define the first order derivation form of loss function.
#'
#' @param type The type of loss function.
#' If \code{type = 'LS'}, return the first-derivative of least square loss function.
#' If \code{type = 'LAD'}, return the first-derivative of least absolute deviation loss function.
#' If \code{type = 'Huber'}, return the first-derivative of Huber loss function.
#' @param k.Huber Constant for Huber loss function. The default is \code{k.Huber = 1.345}.
#'
#' @return The first derivative form of loss function, a \code{function} object.
#'

DiffLoss <- function(type = 'LS', k.Huber = 1.345) {

  if (type == 'LS')
    loss.diff <- function(t) t

  if (type == 'LAD')
    loss.diff <- function(t) sign(t)

  if (type == 'Huber') {
    loss.diff <- function(t) {
      ft.LS <- t * (abs(t) <= k.Huber); ft.LAD <- k.Huber * sign(t) * (abs(t) > k.Huber)
      ft.Huber <- ft.LS + ft.LAD
    }
  }

  return(loss.diff)
}

#' Calculate the diagonal matrix for LQA algorithm
#'
#' @param residual The forecast residuals.
#' @param loss.diff The type of loss function.
#' @param varsigma A positive perturbation, the default is \code{varsigma = 0}.
#' @param tau A positive constant, the default is \code{tau = 1e-16}.

DiagMatrix_DiffLoss <- function(residual, loss.diff, varsigma = 1e-8, tau = 1e-16) {

  residuals <- abs(residual)

  vector.diff <- (residuals + varsigma) / (loss.diff(residuals) + tau)
  diag.diff <- diag(as.numeric(vector.diff))

  return(diag.diff)
}

#' Robust optimal reconciliation for hierarchical time series forecasting
#'
#' The RoME model for HTS forecasting. Return the forecast result and computation time.
#'
#' @param base The value of base forecasts.
#' @param hierarchy The hierarchy structure of \code{hts}.
#' @param loss.diff Type of loss function.
#' @param design.W The alternative designs for the covariance matrix of base forecast errors.
#' @param matrix.W The covariance matrix of first-step-ahead base forecasts of error term.
#' @param lambda.shrink Shrink parameter for the shrink covariance design.
#' @param threshold.convergence The convergence threshold, the default set is 1e-4.
#' @param limit.iteration The limit iteration times, the default set is 1000.
#' @param varsigma A positive perturbation, the default is \code{varsigma = 0}.
#' @param tau A positive constant, the default is \code{tau = 1e-16}.
#' @return A list object, which is basically a list consisting of:
#' \item{RoME}{the forecast result of RoME method}
#' \item{time}{the computation time of RoME}
#' \item{call}{a \code{list} object contains arguments for RoME algorithm}
#'
#' @export
#'


RoME_fore <- function(base, hierarchy, loss.diff, design.W = 'OLS', matrix.W = NULL, lambda.shrink = NULL,
                 varsigma = 0, tau = 1e-16, threshold.convergence = 1e-4, limit.iteration = 1e+3) {

  term <- nrow(base); size.m <- nrow(hierarchy); size.n <- ncol(hierarchy); size.mstar <- size.m - size.n

  if (design.W == 'OLS')
    sqrt.W <- diag(1, nrow = size.m, ncol = size.m)

  if (design.W == 'WLSv')
    sqrt.W <- diag(sqrt(diag(matrix.W)))

  if (design.W == 'WLSs')
    sqrt.W <- diag(sqrt(rowSums(hierarchy)))

  if (design.W == 'Sample')
    sqrt.W <- SqrtMatrix(matrix.W)

  if (design.W == 'Shrink') {
    sqrt.W <- (1 - lambda.shrink) * matrix.W; diag(sqrt.W) <- diag(matrix.W)
    sqrt.W <- SqrtMatrix(sqrt.W)
  }

  forecasts.BU <- base[, -(1:size.mstar)] %*% t(hierarchy)

  time.start <- Sys.time()

  J <- cbind(matrix(0, nrow = size.n, ncol = size.mstar), diag(1, nrow = size.n, ncol = size.n))
  tU <- cbind(diag(1, nrow = size.mstar, ncol = size.mstar), rbind(-hierarchy[1:size.mstar, ]))

  forecasts.RoME <- matrix(nrow = term, ncol = size.m)
  for (j in 1:term) {

    x.begin <- rep(0, size.m); x.end <- forecasts.BU[j, ]
    for (k in 1:limit.iteration) {

      cost.begin <- x.begin - base[j, ]

      D <- DiagMatrix_DiffLoss(residual = cost.begin, loss.diff = loss.diff, varsigma = varsigma)
      WDW <- sqrt.W %*% D %*% sqrt.W
      tUWDWU <- tU %*% WDW %*% t(tU); diag(tUWDWU) <- diag(tUWDWU) + tau

      inv.tUWDWU <- try(solve(tUWDWU), silent = TRUE)
      if (!is.numeric(inv.tUWDWU))
        break

      P <- J - J %*% WDW %*% t(tU) %*% inv.tUWDWU %*% tU
      x.end <- hierarchy %*% P %*% base[j, ]

      test.convergence <- max(abs(x.end - x.begin))
      if (test.convergence < threshold.convergence)
        break

      x.begin <- x.end

      # if (k == limit.iteration)
      #   cat('Out of the iteration limit! Method: ', design.W, '.\n', sep = '')
    }

    forecasts.RoME[j, ] <- x.end
  }

  time.end <- Sys.time(); time.RoME <- as.numeric(time.end - time.start) / term

  list.return <- list(RoME = forecasts.RoME,
                      time = time.RoME,
                      call = list(loss = loss.diff,
                                  covariance = matrix.W,
                                  design = design.W,
                                  lambda = lambda.shrink))
  return(list.return)
}
