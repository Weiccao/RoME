####==== Loss Function ====####

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

DiagMatrix_DiffLoss <- function(residual, loss.diff, varsigma = 1e-8, tau = 1e-16) {
  
  residuals <- abs(residual)
  
  vector.diff <- (residuals + varsigma) / (loss.diff(residuals) + tau)
  diag.diff <- diag(as.numeric(vector.diff))
  
  return(diag.diff)
}

SqrtMatrix <- function(x, tau = 1e-16) {
  
  sol.eigen <- eigen(x)
  
  eigen.min <- min(sol.eigen$values)
  if (eigen.min <= 0)
    sol.eigen$values <- sol.eigen$values - eigen.min + tau
  
  eigen.sqrt <- sqrt(sol.eigen$values)
  x.sqrt <- sol.eigen$vectors %*% diag(eigen.sqrt) %*% t(sol.eigen$vectors)
  
  return(x.sqrt)
}

####==== Basic Forecast ====####

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

####==== Reconciliation Forecast ====####

Hierarchy2Nodes <- function(S) {
  
  size.m <- nrow(S); size.n <- ncol(S)
  
  indicators.level <- cumsum(rowSums(S)) %% size.n; ends.level <- (1:size.m)[!indicators.level]
  
  nodes.all <- list(); number.level <- length(ends.level)
  for (i in 2:number.level) {
    
    if (i == 2)
      S.parent <- t(rbind(S[1:ends.level[i - 1], ]))
    else
      S.parent <- t(S[(ends.level[i - 2] + 1):ends.level[i - 1], ])
    
    S.child <- t(S[(ends.level[i - 1] + 1):ends.level[i], ])
    
    nodes.j <- NULL; count.tmp <- 0; S.tmp <- rep(0, size.n)
    for (j in 1:ncol(S.child)) {
      
      count.tmp <- count.tmp + 1; S.tmp <- S.tmp + S.child[, j]
      
      test.new <- sum(apply(S.parent == S.tmp, 2, prod))
      if (test.new) {
        nodes.j <- c(nodes.j, count.tmp)
        count.tmp <- 0; S.tmp <- rep(0, size.n)
      }
    }
    
    nodes.all <- c(nodes.all, list(nodes.j))
  }
  
  return(nodes.all)
}

Nodes2Hierarchy <- function(nodes) {
  
  number.level <- length(nodes); size.n <- sum(nodes[[number.level]])
  
  S <- S.child <- diag(1, size.n)
  for (i in number.level:1) {
    
    S.parent <- NULL
    for (j in nodes[[i]]) {
      S.parent <- rbind(S.parent, colSums(rbind(S.child[1:j, ])))
      S.child <- S.child[-(1:j), ]
    }
    
    S <- rbind(S.parent, S); S.child <- S.parent
  }
  
  return(S)
}

BottomUp <- function(base, hierarchy) {
  
  size.mstar <- nrow(hierarchy) - ncol(hierarchy)
  
  sol.bu <- base[, -(1:size.mstar)] %*% t(hierarchy)
  
  return(sol.bu)
}

RoME <- function(base, hierarchy, loss.diff, design.W = 'OLS', matrix.W = NULL, lambda.shrink = NULL,
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

####==== Generation & Summary ====####

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