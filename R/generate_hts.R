#' Function to generate error terms with different distribution.
#'
#' @param type The distribution of error term.
#' @param size The number of observations.
#'
#' @return error term with specific distribution.

GenerateError <- function(type = 'Normal', size = 1) {

  if (type == 'Normal')
    errors <- rnorm(n = size, mean = 0, sd = 1)

  if (type == 'Mixture') {

    indicators.outlier <- as.logical(rbinom(n = size, size = 1, prob = .1)); number.outlier <- sum(indicators.outlier)

    errors <- rnorm(n = size, mean = 0, sd = 1)
    errors[indicators.outlier] <- rnorm(n = number.outlier, mean = 0, sd = 3)
  }

  if (type == 'Tstudent')
    errors <- rt(n = size, df = 3, ncp = 0)

  if (type == 'Cauchy')
    errors <- rcauchy(n = size, location = 0, scale = 1)

  return(errors)
}



#' Function to generate \code{hts} with non-Gaussian residuals.
#'
#' @param size.n The total number of \code{ts}.
#' @param size.c The number of \code{ts} with non-Gaussian residuals.
#' @param size.T The number of observations of each \code{ts}.
#' @param alpha Value of the autoregressive coefficient.
#' @param rho  The correlation coefficient parameter between any two standard normality \code{ts}.
#' @param type.epsilon The distribution of non-Gaussian error term.
#' @param times.repeat Replication time.
#'
#' @return \code{hts} with non-Gaussian residuals

Generate_NonGaussianSeries <- function(size.n = 9, size.c = 3, size.T = 192, alpha = .8, rho = .4, type.epsilon, times.repeat = 1e+3) {

  size.s <- size.n - size.c; mu.error <- rep(0, size.s); Sigma.error <- matrix(nrow = size.s, ncol = size.s)
  for (i in 1:size.s)
    for (j in 1:size.s)
      Sigma.error[i, j] <- rho^abs(i - j)

  for (t in 1:times.repeat) { cat('Non-Gaussian series. Replication: ', t, '.\n', sep = '')
    for (i in type.epsilon) {

      errors.c <- matrix(nrow = size.T, ncol = size.c)
      for (j in 1:size.c)
        errors.c[, j] <- GenerateError(type = i, size = size.T)

      errors.s <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error)

      y <- cbind(errors.c, errors.s); y <- .5 * y

      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]

      write.table(y, paste('NonGaussianSeries_epsilon', i, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

#' Function to generate \code{hts} with different variation parameter
#'
#' @param size.n The total number of \code{ts}.
#' @param size.T The number of observations of each \code{ts}.
#' @param alpha Value of the autoregressive coefficient.
#' @param rho  The correlation coefficient parameter between any two standard normality \code{ts}.
#' @param sigma The value of variation parameter.
#' @param times.repeat Replication time.
#'
#' @return \code{hts} with different variation parameter.

Generation_VariationOfNormality <- function(size.n = 6, size.T = 192, alpha = .6, rho = .4, sigma, times.repeat = 1e+3) {

  mu.error <- rep(0, size.n); Sigma.error <- matrix(nrow = size.n, ncol = size.n)
  for (i in 1:size.n)
    for (j in 1:size.n)
      Sigma.error[i, j] <- rho^abs(i - j)

  for (t in 1:times.repeat) { cat('Variation Of normality. Replication: ', t, '.\n', sep = '')
    for (i in sigma) {

      y <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error); y <- i * y

      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]

      write.table(y, paste('VariationOfNormality_sigma', i * 10, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

#' Function to generate \code{hts} with different hierarchy scale.
#'
#' @param size.n The total number of \code{ts}.
#' @param size.T The number of observations of each \code{ts}.
#' @param alpha Value of the autoregressive coefficient.
#' @param rho  The correlation coefficient parameter between any two standard normality \code{ts}.
#' @param times.repeat Replication time
#'
#' @return \code{hts} with different hierarchy scale.

Generation_ScaleOfHierarchy <- function(size.n = seq(10,50,5), size.T = 192, alpha = .8, rho = .4, times.repeat = 1e+3) {

  for (i in size.n) {

    mu.error <- rep(0, i); Sigma.error <- matrix(nrow = i, ncol = i)
    for (j in 1:i)
      for (k in 1:i)
        Sigma.error[j, k] <- rho^abs(j - k)

    for (t in 1:times.repeat) { cat('Scale of hierarchy. n: ', i, '. Replication: ', t, '.\n', sep = '')

      y <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error); y = .5 * y

      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]

      write.table(y, paste('ScaleOfHierarchy_n', i, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

#' Function to generate \code{hts} with different bad forecasts proportion.
#'
#' @param size.n The total number of \code{ts}.
#' @param size.c The number of \code{ts} with non-Gaussian residuals.
#' @param size.T The number of observations of each \code{ts}.
#' @param type.epsilon The distribution of non-Gaussian error term.
#' @param range.alpha Lower and upper limits of the autoregressive coefficient which draw from \code{unif} distribution.
#' @param times.repeat Replication time
#'
#' @return \code{hts} with different bad forecasts proportion.

Generation_ProportionOfBadForecasts <- function(size.n = 30, size.c, size.T = 192, type.epsilon = c('Mixture', 'Tstudent', 'Cauchy'),
                                                range.alpha = c(.6, .8), times.repeat = 1e+3) {

  for (i in size.c) {

    size.s <- size.n - i; mu.error <- rep(0, size.s); Sigma.error <- diag(1, size.s)

    for (t in 1:times.repeat) { cat('Proportion of bad forecasts. size.c: ', i, '. Replication: ', t, '.\n', sep = '')

      epsilons.c <- sample(type.epsilon, size = i, replace = TRUE); errors.c <- matrix(nrow = size.T, ncol = i)
      for (j in 1:i)
        errors.c[, j] <- GenerateError(type = epsilons.c[j], size = size.T)

      errors.s <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error)

      y <- cbind(errors.c, errors.s); y <- .5 * y

      alpha <- runif(n = size.n, min = range.alpha[1], max = range.alpha[2])
      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]

      write.table(y, paste('ProportionOfBadForecasts_changeable', i, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

#' Function to generate \code{hts} with different correlation coefficient.
#'
#' @param size.n The total number of \code{ts}.
#' @param size.T The number of observations of each \code{ts}.
#' @param alpha value of autoregressive coeffcient.
#' @param rho the correlation coeffcient parameter between any two bottom-level ts
#' @param times.repeat the replication time
#'
#' @return \code{hts} with different correlation coefficient.

Generation_Correlation <- function(size.n = 9, size.T = 192, alpha = .6, rho, times.repeat = 1e+3) {

  mu.error <- rep(0, size.n)

  for (i in rho) {

    Sigma.error <- matrix(i, nrow = size.n, ncol = size.n); diag(Sigma.error) <- 1

    for (t in 1:times.repeat) { cat('Correlation. rho: ', i, '. Replication: ', t, '.\n', sep = '')

      y <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error); y = .5 * y

      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]

      write.table(y, paste('Correlation_rho', i * 10, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

#' Function to generate complexity \code{hts} (change in scale, bad forecast proportion, etc.)
#'
#' @param size.n The total number of \code{ts}.
#' @param size.T The number of observations of each \code{ts}.
#' @param ratio.c Ratio of non Gaussian series.
#' @param range.alpha Lower and upper limits of the autoregressive coefficient which draw from \code{unif} distribution.
#' @param times.repeat Replication time
#'
#' @return \code{hts} with complexity structure.

Generation_LargeScale <- function(size.n, ratio.c = .4, size.T = 192, type.epsilon = c('Mixture', 'Tstudent', 'Cauchy'),
                                  range.alpha = c(.6, .8), times.repeat = 1e+3) {

  for (i in size.n) {

    size.c <- floor(i * ratio.c); size.s <- i - size.c; mu.error <- rep(0, size.s); Sigma.error <- diag(1, size.s)

    for (t in 1:times.repeat) { cat('Large scale. n: ', i, '. Replication: ', t, '.\n', sep = '')

      indicators.c <- sample(1:i, size = size.c); y <- matrix(nrow = size.T, ncol = i)
      for (j in indicators.c)
        y[, j] <- GenerateError(type = sample(type.epsilon, size = 1), size = size.T)

      y[, -indicators.c] <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error); y <- .5 * y

      alpha <- runif(n = i, min = range.alpha[1], max = range.alpha[2])
      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]

      write.table(y, paste('LargeScale_n', i, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}


