####==== Global Environment ====####

rm(list = ls())
setwd('./RoME')
library(MASS)
source('RoME_Basic.R')

####==== Experiment 1:  Non-Gaussian Series ====####

Generate_NonGaussianSeries <- function(size.n = 9, size.c = 3, size.T = 192, alpha = .8, rho = .4, type.epsilon, times.repeat = 1e+3) {
  
  size.s <- size.n - size.c; mu.error <- rep(0, size.s); Sigma.error <- matrix(nrow = size.s, ncol = size.s)
  for (i in 1:size.s) {
    for (j in 1:size.s) {
      Sigma.error[i, j] <- rho^abs(i - j)
    }
  }
  for (t in 1:times.repeat) { 
    cat('Non-Gaussian series. Replication: ', t, '.\n', sep = '')
    for (idx_i in seq_along(type.epsilon)) {
      i <- type.epsilon[idx_i]
      # 后续的 GenerateError 和 mvrnorm 都会顺着这个种子自然生成，既独立又可复现。
      set.seed(10000 * t + idx_i)
      errors.c <- matrix(nrow = size.T, ncol = size.c)
      for (j in 1:size.c) {
        errors.c[, j] <- GenerateError(type = i, size = size.T)
      }
      errors.s <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error)
      y <- cbind(errors.c, errors.s) 
      y <- .5 * y
      
      for (j in 2:size.T) {
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]
      }
        
        write.table(y, paste('NonGaussianSeries_epsilon', i, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

#Generate_NonGaussianSeries(type.epsilon = c('Normal','Mixture', 'Tstudent'))

####==== Experiment 2:  Proportion of Bad Forecasts ====####

Generation_ProportionOfBadForecasts <- function(size.n = 30, size.c, size.T = 192, type.epsilon = c('Mixture', 'Tstudent', 'Cauchy'),
                                                range.alpha = c(.6, .8), times.repeat = 1e+3) {
  
  for (i in size.c) {
    
    size.s <- size.n - i; mu.error <- rep(0, size.s); Sigma.error <- diag(1, size.s)
    
    for (t in 1:times.repeat) { cat('Proportion of bad forecasts. size.c: ', i, '. Replication: ', t, '.\n', sep = '')
      set.seed(123+t)
      
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

# Generation_ProportionOfBadForecasts(size.c = 1:9 * 3)

####==== Experiment 3:  Correlation between normal series ====####

Generation_Correlation <- function(size.n = 9, size.T = 192, alpha = .6, rho, times.repeat = 1e+3) {
  
  mu.error <- rep(0, size.n)
  
  for (i in rho) {
    
    Sigma.error <- matrix(i, nrow = size.n, ncol = size.n); diag(Sigma.error) <- 1
    
    for (t in 1:times.repeat) { cat('Correlation. rho: ', i, '. Replication: ', t, '.\n', sep = '')
      set.seed(123+t)
      y <- mvrnorm(n = size.T, mu = mu.error, Sigma = Sigma.error); y = .5 * y
      
      for (j in 2:size.T)
        y[j, ] <- alpha * y[j - 1, ] + y[j, ]
      
      write.table(y, paste('Correlation_rho', i * 10, 'y.txt', sep = '_'), append = TRUE, row.names = FALSE, col.names = FALSE)
    }
  }
}

# Generation_Correlation(rho = 0:9 * 0.1)

####==== Experiment 4:  Complexity of hierarchy ====####

Generation_LargeScale <- function(size.n, ratio.c = .4, size.T = 192, type.epsilon = c('Mixture', 'Tstudent', 'Cauchy'), 
                                  range.alpha = c(.6, .8), times.repeat = 1e+3) {
  
  for (i in size.n) {
    
    size.c <- floor(i * ratio.c); size.s <- i - size.c; mu.error <- rep(0, size.s); Sigma.error <- diag(1, size.s)
    
    for (t in 1:times.repeat) { cat('Large scale. n: ', i, '. Replication: ', t, '.\n', sep = '')
      set.seed(123+t)
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

# Generation_LargeScale(size.n = 1:6 * 20)