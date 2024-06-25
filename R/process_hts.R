#' generate non-Gaussain hts and forecast with RoME
#'
#' @param type.epsilon distribution of error term.
#' @param times.repeat the replication time.
#' @param size.T number of observations for each time series.
#' @param term forecast step.
#' @param method.forecast forecast method for base forecasts.
#' @param design.W type of design matrix.
#'
#' @return the RoME forecast result

Process_NonGaussianSeries <- function(type.epsilon = c('Mixture', 'Tstudent', 'Cauchy'), times.repeat = 1e+3,
                                      size.T = 192, term = 12, method.forecast, design.W) {

  hierarchy <- rbind(rep(1, 9), c(rep(1, 3), rep(0, 6)), c(rep(0, 3), rep(1, 3), rep(0, 3)), c(rep(0, 6), rep(1, 3)), diag(1, 9))

  for (i in type.epsilon) {

    y <- as.matrix(read.table(paste('NonGaussianSeries_epsilon', i, 'y.txt', sep = '_')))

    for (t in 1:times.repeat) { cat('Non-Gaussian series. epsilon: ', i, '. Replication: ', t, '.\n', sep = '')

      rows.t <- 1:size.T + (t - 1) * size.T; indicators.train <- 1:(size.T - term)
      yt <- y[rows.t, ] %*% t(hierarchy); yt.train <- yt[indicators.train, ]; yt.test <- yt[-indicators.train, ]

      sol.Wt <- CovarianceW(x = yt.train, method.forecast = method.forecast, is.rolling = FALSE)

      sol.base <- BaseAuto(x = yt.train, method.forecast = method.forecast, term = term)
      write.table(sol.base$base - yt.test,
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'base_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      sol.BU <- BottomUp(base = sol.base$base, hierarchy = hierarchy)
      write.table(sol.BU - yt.test,
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'BU_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      time.LS <- time.LAD.1 <- time.LAD.2 <- time.Huber <- NULL
      for (j in design.W) {

        sol.LS <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'),
                       design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LS$RoME - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'LS', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD.1 <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LAD'),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda, varsigma = 1e-8)
        write.table(sol.LAD.1$RoME - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'LAD1', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD.2 <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1e-4),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LAD.2$RoME - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'LAD2', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.Huber <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1.345),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.Huber$RoME - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'Huber', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        write.table((sol.LS$RoME + sol.LAD.2$RoME) / 2 - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'Averaged', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table((sol.LS$RoME + sol.LAD.2$RoME * 1:term) / (1 + 1:term) - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'Single', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table((sol.LS$RoME * term:1 + sol.LAD.2$RoME * 1:term) / (term + 1) - yt.test,
                    paste('NonGaussianSeries_epsilon', i, method.forecast, 'Double', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        time.LS <- c(time.LS, sol.LS$time)
        time.LAD.1 <- c(time.LAD.1, sol.LAD.1$time); time.LAD.2 <- c(time.LAD.2, sol.LAD.2$time)
        time.Huber <- c(time.Huber, sol.Huber$time)
      }

      write.table(t(time.LS),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'LS_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LAD.1),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'LAD1_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LAD.2),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'LAD2_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.Huber),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'Huber_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LS + time.LAD.1),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'Averaged_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LS + time.LAD.2),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'Single_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LS + time.LAD.1),
                  paste('NonGaussianSeries_epsilon', i, method.forecast, 'Double_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}


#' generate HTS with varies variation parameter and forecast with RoME
#'
#' @param sigma the sd value of time series
#' @param times.repeat the replication time.
#' @param size.T number of observations for each time series.
#' @param term forecast step.
#' @param method.forecast forecast method for base forecasts.
#' @param design.W type of design matrix.
#'
#' @return the RoME forecast result
#'
Process_VariationOfNormality <- function(sigma = c(.5, 1, 1.5, 2, 3), times.repeat = 1e+3, size.T = 192, term = 12, method.forecast, design.W) {

  hierarchy <- rbind(rep(1, 6), c(rep(1, 3), rep(0, 3)), c(rep(0, 3), rep(1, 3)), diag(1, 6))

  for (i in sigma) {

    y <- as.matrix(read.table(paste('VariationOfNormality_sigma', i * 10, 'y.txt', sep = '_')))

    for (t in 1:times.repeat) { cat('Variation of normality. sigma: ', i, '. Replication: ', t, '.\n', sep = '')

      rows.t <- 1:size.T + (t - 1) * size.T; indicators.train <- 1:(size.T - term)
      yt <- y[rows.t, ] %*% t(hierarchy); yt.train <- yt[indicators.train, ]; yt.test <- yt[-indicators.train, ]

      sol.Wt <- CovarianceW(x = yt.train, method.forecast = method.forecast, is.rolling = FALSE)

      sol.base <- BaseAuto(x = yt.train, method.forecast = method.forecast, term = term)
      write.table(sol.base$base - yt.test,
                  paste('VariationOfNormality_sigma', i * 10, method.forecast, 'base_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      sol.BU <- BottomUp(base = sol.base$base, hierarchy = hierarchy)
      write.table(sol.BU - yt.test,
                  paste('VariationOfNormality_sigma', i * 10, method.forecast, 'BU_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      for (j in design.W) {

        sol.LS <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'),
                       design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LS$RoME - yt.test,
                    paste('VariationOfNormality_sigma', i * 10, method.forecast, 'LS', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1e-4),
                        design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LAD$RoME - yt.test,
                    paste('VariationOfNormality_sigma', i * 10, method.forecast, 'LAD', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.Huber <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1.345),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.Huber$RoME - yt.test,
                    paste('VariationOfNormality_sigma', i * 10, method.forecast, 'Huber', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

      }
    }
  }
}

#' generate HTS with different hierarchy scale and forecast with RoME
#'
#' @param size.n the total number of time series.
#' @param size.each number of series in each node.
#' @param size.T number of observations for each time series.
#' @param times.repeat the replication time
#' @param term forecast step.
#' @param method.forecast forecast method for base forecasts.
#' @param design.W type of design matrix.
#'
#' @return the RoME forecast result
#'

Process_ScaleOfHierarchy <- function(size.n = 1:5 * 10, size.each = 5, times.repeat = 1e+3, size.T = 192, term = 12, method.forecast, design.W) {

  for (i in size.n) {

    y <- as.matrix(read.table(paste('ScaleOfHierarchy_n', i, 'y.txt', sep = '_')))

    number.middle <- i / size.each; nodes <- list(number.middle, rep(size.each, number.middle))
    hierarchy <- Nodes2Hierarchy(nodes = nodes)

    for (t in 1:times.repeat) { cat('Scale of hierarchy. n: ', i, '. Replication: ', t, '.\n', sep = '')

      rows.t <- 1:size.T + (t - 1) * size.T; indicators.train <- 1:(size.T - term)
      yt <- y[rows.t, ] %*% t(hierarchy); yt.train <- yt[indicators.train, ]; yt.test <- yt[-indicators.train, ]

      sol.Wt <- CovarianceW(x = yt.train, method.forecast = method.forecast, is.rolling = FALSE)

      sol.base <- BaseAuto(x = yt.train, method.forecast = method.forecast, term = term)
      write.table(sol.base$base - yt.test,
                  paste('ScaleOfHierarchy_n', i, method.forecast, 'base_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      sol.BU <- BottomUp(base = sol.base$base, hierarchy = hierarchy)
      write.table(sol.BU - yt.test,
                  paste('ScaleOfHierarchy_n', i, method.forecast, 'BU_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      time.LS <- time.LAD <- time.Huber <- NULL
      for (j in design.W) {

        sol.LS <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'),
                       design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LS$RoME - yt.test,
                    paste('ScaleOfHierarchy_n', i, method.forecast, 'LS', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1e-4),
                        design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LAD$RoME - yt.test,
                    paste('ScaleOfHierarchy_n', i, method.forecast, 'LAD', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.Huber <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1.345),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.Huber$RoME - yt.test,
                    paste('ScaleOfHierarchy_n', i, method.forecast, 'Huber', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        time.LS <- c(time.LS, sol.LS$time); time.LAD <- c(time.LAD, sol.LAD$time); time.Huber <- c(time.Huber, sol.Huber$time)
      }

      write.table(t(time.LS),
                  paste('ScaleOfHierarchy_n', i, method.forecast, 'LS_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LAD),
                  paste('ScaleOfHierarchy_n', i, method.forecast, 'LAD_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.Huber),
                  paste('ScaleOfHierarchy_n', i, method.forecast, 'Huber_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}


#' generate HTS with different proportions of bad forecasts and forecast with RoME
#'
#' @param size.c number of non Gaussian series
#' @param size.T number of observations for each time series.
#' @param times.repeat the replication time
#' @param term forecast step.
#' @param method.forecast forecast method for base forecasts.
#' @param design.W type of design matrix.
#'
#' @return the RoME forecast result
#'

Process_ProportionOfBadForecasts <- function(size.c = 1:9 * 3, times.repeat = 1e+3, size.T = 192, term = 12, method.forecast, design.W) {

  hierarchy <- matrix(0, nrow = 5, ncol = 30)
  for (i in 1:5)
    hierarchy[i, 1:6 + (i - 1) * 6] <- 1
  hierarchy <- rbind(rep(1, 30), hierarchy, diag(1, 30))

  for (i in size.c) {

    y <- as.matrix(read.table(paste('ProportionOfBadForecasts_changeable', i, 'y.txt', sep = '_')))

    for (t in 1:times.repeat) { cat('Proportion of bad forecasts. sice.c: ', i, '. Replication: ', t, '.\n', sep = '')

      rows.t <- 1:size.T + (t - 1) * size.T; indicators.train <- 1:(size.T - term)
      yt <- y[rows.t, ] %*% t(hierarchy); yt.train <- yt[indicators.train, ]; yt.test <- yt[-indicators.train, ]

      sol.Wt <- CovarianceW(x = yt.train, method.forecast = method.forecast, is.rolling = FALSE)

      sol.base <- BaseAuto(x = yt.train, method.forecast = method.forecast, term = term)
      write.table(sol.base$base - yt.test,
                  paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'base_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      sol.BU <- BottomUp(base = sol.base$base, hierarchy = hierarchy)
      write.table(sol.BU - yt.test,
                  paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'BU_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      for (j in design.W) {

        sol.LS <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'),
                       design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LS$RoME - yt.test,
                    paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'LS', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1e-4),
                        design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LAD$RoME - yt.test,
                    paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'LAD', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.Huber <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1.345),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.Huber$RoME - yt.test,
                    paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'Huber', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}


#' generate HTS with complicate covariance structure and forecast with RoME
#'
#' @param rho the correlation coeffcient parameter between any two bottom-level ts
#' @param size.T number of observations for each time series.
#' @param times.repeat the replication time
#' @param term forecast step.
#' @param method.forecast forecast method for base forecasts.
#' @param design.W type of design matrix.
#'
#' @return the RoME forecast result
#'

Process_Correlation <- function(rho = 0:9 * .1, times.repeat = 1e+3, size.T = 192, term = 12, method.forecast, design.W) {

  hierarchy <- rbind(rep(1, 9), c(rep(1, 3), rep(0, 6)), c(rep(0, 3), rep(1, 3), rep(0, 3)), c(rep(0, 6), rep(1, 3)), diag(1, 9))

  for (i in rho) {

    y <- as.matrix(read.table(paste('Correlation_rho', i * 10, 'y.txt', sep = '_')))

    for (t in 1:times.repeat) { cat('Correlation. rho: ', i, '. Replication: ', t, '.\n', sep = '')

      rows.t <- 1:size.T + (t - 1) * size.T; indicators.train <- 1:(size.T - term)
      yt <- y[rows.t, ] %*% t(hierarchy); yt.train <- yt[indicators.train, ]; yt.test <- yt[-indicators.train, ]

      sol.Wt <- CovarianceW(x = yt.train, method.forecast = method.forecast, is.rolling = FALSE)

      sol.base <- BaseAuto(x = yt.train, method.forecast = method.forecast, term = term)
      write.table(sol.base$base - yt.test,
                  paste('Correlation_rho', i * 10, method.forecast, 'base_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      sol.BU <- BottomUp(base = sol.base$base, hierarchy = hierarchy)
      write.table(sol.BU - yt.test,
                  paste('Correlation_rho', i * 10, method.forecast, 'BU_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      for (j in design.W) {

        sol.LS <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'),
                       design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LS$RoME - yt.test,
                    paste('Correlation_rho', i * 10, method.forecast, 'LS', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1e-4),
                        design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LAD$RoME - yt.test,
                    paste('Correlation_rho', i * 10, method.forecast, 'LAD', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.Huber <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1.345),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.Huber$RoME - yt.test,
                    paste('Correlation_rho', i * 10, method.forecast, 'Huber', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}


#' generate generate complexity hts (change in scale, bad forecast proportion, etc.) and forecast with RoME
#'
#' @param size.n the total number of time series.
#' @param size.each number of series in each node
#' @param size.T number of observations for each time series.
#' @param times.repeat the replication time
#' @param term forecast step.
#' @param method.forecast forecast method for base forecasts.
#' @param design.W type of design matrix.
#'
#' @return the RoME forecast result
#'

Process_LargeScale <- function(size.n = 1:6 * 20, size.each = 10, times.repeat = 1e+3, size.T = 192, term = 12, method.forecast, design.W) {

  for (i in size.n) {

    y <- as.matrix(read.table(paste('LargeScale_n', i, 'y.txt', sep = '_')))

    number.middle <- i / size.each; nodes <- list(number.middle, rep(size.each, number.middle))
    hierarchy <- Nodes2Hierarchy(nodes = nodes)

    for (t in 1:times.repeat) { cat('Large scale. n: ', i, '. Replication: ', t, '.\n', sep = '')

      rows.t <- 1:size.T + (t - 1) * size.T; indicators.train <- 1:(size.T - term)
      yt <- y[rows.t, ] %*% t(hierarchy); yt.train <- yt[indicators.train, ]; yt.test <- yt[-indicators.train, ]

      sol.Wt <- CovarianceW(x = yt.train, method.forecast = method.forecast, is.rolling = FALSE)

      sol.base <- BaseAuto(x = yt.train, method.forecast = method.forecast, term = term)
      write.table(sol.base$base - yt.test,
                  paste('LargeScale_n', i, method.forecast, 'base_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      sol.BU <- BottomUp(base = sol.base$base, hierarchy = hierarchy)
      write.table(sol.BU - yt.test,
                  paste('LargeScale_n', i, method.forecast, 'BU_residuals.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)

      time.LS <- time.LAD <- time.Huber <- NULL
      for (j in design.W) {

        sol.LS <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'),
                       design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LS$RoME - yt.test,
                    paste('LargeScale_n', i, method.forecast, 'LS', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.LAD <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1e-4),
                        design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.LAD$RoME - yt.test,
                    paste('LargeScale_n', i, method.forecast, 'LAD', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        sol.Huber <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'Huber', k.Huber = sol.Wt$sd * 1.345),
                          design.W = j, matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
        write.table(sol.Huber$RoME - yt.test,
                    paste('LargeScale_n', i, method.forecast, 'Huber', j, 'residuals.txt', sep = '_'),
                    append = TRUE, row.names = FALSE, col.names = FALSE)

        time.LS <- c(time.LS, sol.LS$time); time.LAD <- c(time.LAD, sol.LAD$time); time.Huber <- c(time.Huber, sol.Huber$time)
      }

      write.table(t(time.LS),
                  paste('LargeScale_n', i, method.forecast, 'LS_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.LAD),
                  paste('LargeScale_n', i, method.forecast, 'LAD_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(time.Huber),
                  paste('LargeScale_n', i, method.forecast, 'Huber_time.txt', sep = '_'),
                  append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}
