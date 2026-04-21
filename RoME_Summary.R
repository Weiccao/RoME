####==== Global Environment ====####
rm(list = ls())
setwd('./RoME')
source('RoME_Basic.R')

library(forecast)

####==== Experiment 1: Proportion of Bad Forecasts ====####

Summary_NonGaussianSeries <- function(times.repeat = 1e+3, type.epsilon = c('Normal','Mixture', 'Tstudent'), method.forecast = 'ARIMA',
                                      term = 12, range.term, range.m,
                                      type.loss = c('LS', 'LAD1', 'LAD2', 'Huber', 'Averaged', 'Single', 'Double'),
                                      design.W = c('OLS', 'WLSv', 'WLSs', 'Sample', 'Shrink')) {
  
  length.term <- length(range.term)
  length.m <- length(range.m)
  number.loss <- length(type.loss)
  number.design <- length(design.W)
 
  Get_MSE_Matrix <- function(res_matrix, is_cauchy) {
    size.m <- ncol(res_matrix)
    err_mat <- matrix(0, nrow = term, ncol = size.m)
    for (l in 1:term) {
      rows_l <- l + (0:(times.repeat - 1)) * term
      sq_err <- res_matrix[rows_l, , drop = FALSE]^2
      if (is_cauchy) {
        err_mat[l, ] <- apply(sq_err, 2, median) # Median
      } else {
        err_mat[l, ] <- colMeans(sq_err)         
      }
    }
    return(err_mat)
  }
  
  Get_Loss_Cube <- function(res_matrix, is_cauchy) {
    cube <- array(0, dim = c(times.repeat, length.term, length.m))
    for (h in 1:length.term) {
      for (m_idx in 1:length.m) {
        for (t in 1:times.repeat) {
          rows_t <- range.term[[h]] + (t - 1) * term
          sub_err <- res_matrix[rows_t, range.m[[m_idx]], drop = FALSE]^2
          if (is_cauchy) {
            cube[t, h, m_idx] <- median(sub_err)
          } else {
            cube[t, h, m_idx] <- mean(sub_err)
          }
        }
      }
    }
    return(cube)
  }
  
  method_names <- c(paste(rep(type.loss, each = number.design), rep(design.W, times = number.loss), sep = '-'), 'BU', 'base')
  num_methods <- length(method_names)
  
  for (i in type.epsilon) {
    cat("\n================ Processing Epsilon:", i, "================\n")
    is_cauchy <- (i == 'Cauchy')
    
    All_MSE_Mats <- list()
    All_Loss_Cubes <- list()
    matrix.time <- NULL
    
    # 1. load Base-forecast 和 Bottom-Up
    res.base <- as.matrix(read.table(paste('NonGaussianSeries_epsilon', i, method.forecast, 'base_residuals.txt', sep = '_')))
    All_MSE_Mats[["base"]] <- Get_MSE_Matrix(res.base, is_cauchy)
    All_Loss_Cubes[["base"]] <- Get_Loss_Cube(res.base, is_cauchy)
    
    res.BU <- as.matrix(read.table(paste('NonGaussianSeries_epsilon', i, method.forecast, 'BU_residuals.txt', sep = '_')))
    All_MSE_Mats[["BU"]] <- Get_MSE_Matrix(res.BU, is_cauchy)
    All_Loss_Cubes[["BU"]] <- Get_Loss_Cube(res.BU, is_cauchy)
    
    # 2. load competitive methods
    for (j in 1:number.loss) {
      for (k in 1:number.design) {
        m_name <- paste(type.loss[j], design.W[k], sep = '-')
        res.jk <- as.matrix(read.table(paste('NonGaussianSeries_epsilon', i, method.forecast, type.loss[j], design.W[k], 'residuals.txt', sep = '_')))
        
        All_MSE_Mats[[m_name]] <- Get_MSE_Matrix(res.jk, is_cauchy)
        All_Loss_Cubes[[m_name]] <- Get_Loss_Cube(res.jk, is_cauchy)
      }
      time.j <- as.matrix(read.table(paste('NonGaussianSeries_epsilon', i, method.forecast, type.loss[j], 'time.txt', sep = '_')))
      matrix.time <- rbind(matrix.time, colMeans(time.j))
    }
    
    # 3. Initialize metrics matrix
    mat.rmse <- matrix(0, nrow = num_methods, ncol = length.term * length.m)
    mat.avgrelmse <- matrix(0, nrow = num_methods, ncol = length.term * length.m)
    mat.dm_p <- matrix(0, nrow = num_methods, ncol = length.term * length.m)
    
    rownames(mat.rmse) <- method_names; rownames(mat.avgrelmse) <- method_names
    rownames(mat.dm_p) <- method_names
    
    # 4. Estimation and testing
    col_idx <- 1
    for (h in 1:length.term) {
      for (m_idx in 1:length.m) {
        
        for (row_idx in 1:num_methods) {
          m_name <- method_names[row_idx]
          
          # 【A. RMSE & AvgRelMSE】
          sub_err_m    <- All_MSE_Mats[[m_name]][range.term[[h]], range.m[[m_idx]], drop = FALSE]
          sub_err_base <- All_MSE_Mats[["base"]][range.term[[h]], range.m[[m_idx]], drop = FALSE]
          ratio_mat    <- sub_err_m / sub_err_base
          
          if (m_name == "base") {
            if (is_cauchy) {
              mat.rmse[row_idx, col_idx] <- sqrt(median(sub_err_base))
              mat.avgrelmse[row_idx, col_idx] <- median(sub_err_base)
            } else {
              mat.rmse[row_idx, col_idx] <- sqrt(mean(sub_err_base))
              mat.avgrelmse[row_idx, col_idx] <- mean(sub_err_base)
            }
          } else {
            if (is_cauchy) {
              rmse_m_val <- sqrt(median(sub_err_m))
              rmse_b_val <- sqrt(median(sub_err_base))
              mat.rmse[row_idx, col_idx] <- (1 - (rmse_m_val / rmse_b_val)) * 100
              mat.avgrelmse[row_idx, col_idx] <- (1 - median(ratio_mat)) * 100
            } else {
              rmse_m_val <- sqrt(mean(sub_err_m))
              rmse_b_val <- sqrt(mean(sub_err_base))
              mat.rmse[row_idx, col_idx] <- (1 - (rmse_m_val / rmse_b_val)) * 100
              mat.avgrelmse[row_idx, col_idx] <- (1 - mean(ratio_mat)) * 100
            }
          }
          
          # 【B. DM test】
          loss_seq_m <- All_Loss_Cubes[[m_name]][, h, m_idx]
          loss_seq_base <- All_Loss_Cubes[["base"]][, h, m_idx]
          
          if (m_name == "base" || isTRUE(all.equal(loss_seq_m, loss_seq_base))) {
            mat.dm_p[row_idx, col_idx] <- 1.0
          } else {
            dm_res <- forecast::dm.test(loss_seq_base, loss_seq_m, h = 1, power = 1)
            mat.dm_p[row_idx, col_idx] <- dm_res$p.value
          }
        }
        
        col_idx <- col_idx + 1
      }
    }
    
    # 5. Export results
    mat.rmse <- round(mat.rmse, digits = 3)
    mat.avgrelmse <- round(mat.avgrelmse, digits = 3)
    mat.dm_p <- round(mat.dm_p, digits = 4) 
    
    write.csv(mat.rmse, paste('NonGaussianSeries_epsilon', i, method.forecast, 'RMSE.csv', sep = '_'))
    write.csv(mat.avgrelmse, paste('NonGaussianSeries_epsilon', i, method.forecast, 'AvgRelMSE.csv', sep = '_'))
    write.csv(mat.dm_p, paste('NonGaussianSeries_epsilon', i, method.forecast, 'DM_pvalue.csv', sep = '_'))
    matrix.time <- round(matrix.time * 1e+3, digits = 3)
    rownames(matrix.time) <- type.loss; colnames(matrix.time) <- design.W
    write.csv(matrix.time, paste('NonGaussianSeries_epsilon', i, method.forecast, 'time.csv', sep = '_'))
  }
}

#Summary_NonGaussianSeries(range.term = list(1, 1:6, 1:12), range.m = list(c(2, 5:7), c(3:4, 8:13), 1:13))

####==== Experiment 2: Proportion of Bad Forecasts ====####

Summary_ProportionOfBadForecasts <- function(times.repeat = 1e+3, size.c = 1:9 * 3, method.forecast = 'ARIMA', term = 12, range.term, range.m,
                                             type.loss = c('LS', 'LAD', 'Huber'), design.W = c('OLS', 'WLSv', 'WLSs', 'Sample', 'Shrink')) {
  
  length.term <- length(range.term); length.m <- length(range.m); number.loss <- length(type.loss); number.design <- length(design.W)
  
  for (i in size.c) {
    
    residuals.base <- as.matrix(read.table(paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'base_residuals.txt', sep = '_')))
    residuals.BU <- as.matrix(read.table(paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'BU_residuals.txt', sep = '_')))
    
    rmse.base <- rmse.BU <- matrix(nrow = times.repeat, ncol = length.term * length.m)
    for (t in 1:times.repeat) {
      
      rows.t <- 1:term + (t - 1) * term
      
      residuals.base.t <- residuals.base[rows.t, ]
      rmse.base.t <- RMSE(x = residuals.base.t, range.term = range.term, range.m = range.m)
      rmse.base[t, ] <- as.vector(rmse.base.t)
      
      residuals.BU.t <- residuals.BU[rows.t, ]
      rmse.BU.t <- RMSE(x = residuals.BU.t, range.term = range.term, range.m = range.m)
      rmse.BU[t, ] <- as.vector(rmse.BU.t)
    }
    
    rmse.base <- colMeans(rmse.base); rmse.BU <- colMeans(rmse.BU)

    matrix.rmse <- NULL
    for (j in type.loss) {
      
      for (k in design.W) {
        
        residuals.jk <- as.matrix(read.table(paste('ProportionOfBadForecasts_changeable', i, method.forecast, j, k, 'residuals.txt', sep = '_')))
        
        rmse.jk <- matrix(nrow = times.repeat, ncol = length.term * length.m)
        for (t in 1:times.repeat) {
          
          indicators.row <- 1:term + (t - 1) * term; residuals.jk.t <- residuals.jk[indicators.row, ]
          
          rmse.jk.t <- RMSE(x = residuals.jk.t, range.term = range.term, range.m = range.m)
          rmse.jk[t, ] <- as.vector(rmse.jk.t)
        }

        matrix.rmse <- rbind(matrix.rmse, colMeans(rmse.jk))
      }
    }
    
    matrix.rmse <- rbind(matrix.rmse, rmse.BU)
    for (j in 1:nrow(matrix.rmse))
      matrix.rmse[j, ] <- (matrix.rmse[j, ] / rmse.base - 1) * 100
    
    matrix.rmse <- rbind(matrix.rmse, rmse.base); matrix.rmse <- round(matrix.rmse,  digits = 3)
    
    rownames(matrix.rmse) <- c(paste(rep(type.loss, each = number.design), rep(design.W, times = number.loss), sep = '-'), 'BU', 'base')
    write.csv(matrix.rmse, paste('ProportionOfBadForecasts_changeable', i, method.forecast, 'RMSE.csv', sep = '_'))
  }
}

# Summary_ProportionOfBadForecasts(range.term = list(1, 1:6, 1:12), range.m = list(7:36, 1:6, 1:36))

####==== Experiment 3: Correlation between Normal Series ====####

Summary_Correlation <- function(times.repeat = 1e+3, rho = 0:9 * .1, method.forecast, term = 12, range.term, range.m,
                                type.loss = c('LS', 'LAD', 'Huber'), design.W = c('OLS', 'WLSv', 'WLSs', 'Sample', 'Shrink')) {
  
  length.term <- length(range.term); length.m <- length(range.m); number.loss <- length(type.loss); number.design <- length(design.W)
  
  for (i in rho) {
    
    residuals.base <- as.matrix(read.table(paste('Correlation_rho', i * 10, method.forecast, 'base_residuals.txt', sep = '_')))
    residuals.BU <- as.matrix(read.table(paste('Correlation_rho', i * 10, method.forecast, 'BU_residuals.txt', sep = '_')))
    
    rmse.base <- rmse.BU <- matrix(nrow = times.repeat, ncol = length.term * length.m)
    for (t in 1:times.repeat) {
      
      rows.t <- 1:term + (t - 1) * term
      
      residuals.base.t <- residuals.base[rows.t, ]
      rmse.base.t <- RMSE(x = residuals.base.t, range.term = range.term, range.m = range.m)
      rmse.base[t, ] <- as.vector(rmse.base.t)
      
      residuals.BU.t <- residuals.BU[rows.t, ]
      rmse.BU.t <- RMSE(x = residuals.BU.t, range.term = range.term, range.m = range.m)
      rmse.BU[t, ] <- as.vector(rmse.BU.t)
    }
    
    rmse.base <- colMeans(rmse.base); rmse.BU <- colMeans(rmse.BU)
    
    matrix.rmse <- NULL
    for (j in type.loss) {
      
      for (k in design.W) {
        
        residuals.jk <- as.matrix(read.table(paste('Correlation_rho', i * 10, method.forecast, j, k, 'residuals.txt', sep = '_')))
        
        rmse.jk <- matrix(nrow = times.repeat, ncol = length.term * length.m)
        for (t in 1:times.repeat) {
          
          rows.t <- 1:term + (t - 1) * term; residuals.jk.t <- residuals.jk[rows.t, ]
          rmse.jk.t <- RMSE(x = residuals.jk.t, range.term = range.term, range.m = range.m)
          rmse.jk[t, ] <- as.vector(rmse.jk.t)
        }

        matrix.rmse <- rbind(matrix.rmse, colMeans(rmse.jk))
      }
    }
    
    matrix.rmse <- rbind(matrix.rmse, rmse.BU)
    for (j in 1:nrow(matrix.rmse))
      matrix.rmse[j, ] <- (matrix.rmse[j, ] / rmse.base - 1) * 100
    
    matrix.rmse <- rbind(matrix.rmse, rmse.base); matrix.rmse <- round(matrix.rmse,  digits = 3)

    rownames(matrix.rmse) <- c(paste(rep(type.loss, each = number.design), rep(design.W, times = number.loss), sep = '-'), 'BU', 'base')
    write.csv(matrix.rmse, paste('Correlation_rho', method.forecast, i * 10, 'RMSE.csv', sep = '_'))
  }
}

# Summary_Correlation(method.forecast = 'ARIMA', range.term = list(1, 1:6, 1:12), range.m = list(5:13, 1:4, 1:13))
# Summary_Correlation(method.forecast = 'ETS', range.term = list(1, 1:6, 1:12), range.m = list(5:13, 1:4, 1:13))

####==== Experiment 4: Complexity of Hierarchy ====####

Summary_LargeScale <- function(times.repeat = 1e+3, size.n = 1:6 * 20, size.each = 10, method.forecast = 'ARIMA', term = 12, range.term,
                               type.loss = c('LS', 'LAD', 'Huber'), design.W = c('OLS', 'WLSv', 'WLSs', 'Shrink')) {
  
  length.term <- length(range.term); number.loss <- length(type.loss); number.design <- length(design.W)
  
  for (i in size.n) {
    
    number.middle <- i / size.each
    range.m <- list(1:i + 1 + number.middle, 1:(1 + number.middle), 1:(1 + number.middle + i)); length.m <- length(range.m)
    
    residuals.base <- as.matrix(read.table(paste('LargeScale_n', i, method.forecast, 'base_residuals.txt', sep = '_')))
    residuals.BU <- as.matrix(read.table(paste('LargeScale_n', i, method.forecast, 'BU_residuals.txt', sep = '_')))
    
    rmse.base <- rmse.BU <- matrix(nrow = times.repeat, ncol = length.term * length.m)
    for (t in 1:times.repeat) {
      
      rows.t <- 1:term + (t - 1) * term
      
      residuals.base.t <- residuals.base[rows.t, ]
      rmse.base.t <- RMSE(x = residuals.base.t, range.term = range.term, range.m = range.m)
      rmse.base[t, ] <- as.vector(rmse.base.t)
      
      residuals.BU.t <- residuals.BU[rows.t, ]
      rmse.BU.t <- RMSE(x = residuals.BU.t, range.term = range.term, range.m = range.m)
      rmse.BU[t, ] <- as.vector(rmse.BU.t)
    }

        rmse.base <- colMeans(rmse.base); rmse.BU <- colMeans(rmse.BU)
    
    matrix.rmse <- matrix.time <- NULL
    for (j in type.loss) {
      
      for (k in design.W) {
        
        residuals.jk <- as.matrix(read.table(paste('LargeScale_n', i, method.forecast, j, k, 'residuals.txt', sep = '_')))
        
        rmse.jk <- matrix(nrow = times.repeat, ncol = length.term * length.m)
        for (t in 1:times.repeat) {
          
          rows.t <- 1:term + (t - 1) * term; residuals.jk.t <- residuals.jk[rows.t, ]
          rmse.jk.t <- RMSE(x = residuals.jk.t, range.term = range.term, range.m = range.m)
          rmse.jk[t, ] <- as.vector(rmse.jk.t)
        }

                matrix.rmse <- rbind(matrix.rmse, colMeans(rmse.jk))
      }
      
      time.j <- as.matrix(read.table(paste('LargeScale_n', i, method.forecast, j, 'time.txt', sep = '_')))
      matrix.time <- rbind(matrix.time, colMeans(time.j))
    }
    
    matrix.time <- round(matrix.time * 1e+3, digits = 3)
    rownames(matrix.time) <- type.loss; colnames(matrix.time) <- design.W
    
    write.csv(matrix.time, paste('LargeScale_n', i, method.forecast, 'time.csv', sep = '_'))

    matrix.rmse <- rbind(matrix.rmse, rmse.BU)
    for (j in 1:nrow(matrix.rmse))
      matrix.rmse[j, ] <- (matrix.rmse[j, ] / rmse.base - 1) * 100

    matrix.rmse <- rbind(matrix.rmse, rmse.base); matrix.rmse <- round(matrix.rmse,  digits = 3)
    
    rownames(matrix.rmse) <- c(paste(rep(type.loss, each = number.design), rep(design.W, times = number.loss), sep = '-'), 'BU', 'base')
    
    write.csv(matrix.rmse, paste('LargeScale_n', i, method.forecast, 'RMSE.csv', sep = '_'))
  }
}

# Summary_LargeScale(range.term = list(1, 1:6, 1:12))