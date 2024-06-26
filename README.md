# RoME
 Robust optimal reconciliation for hierarchical time series forecasting

 
## Installation
You can install the **development** version from
[Github](https://github.com/Weiccao/RoME)

```s
# install.packages("remotes")
remotes::install_github("Weiccao/RoME")
```

## Usage

```s

# generate hts
y <- Generate_NonGaussianSeries(type.epsilon = 'Mixture',times.repeat = 1)
hierarchy <- rbind(rep(1, 9), c(rep(1, 3), rep(0, 6)),
     c(rep(0, 3), rep(1, 3), rep(0, 3)), c(rep(0, 6), rep(1, 3)), diag(1, 9))
data.train <- as.matrix(y[1:180, ]) %*% t(hierarchy)
data.test <- as.matrix(y[-c(1:180), ]) %*% t(hierarchy)
sol.base <- BaseAuto(data.train, term = 12)

# forecast with RoME
W <- t(sol.base$residual) %*% sol.base$residual / 180
lambda <- LambdaShrink(sol.base$residual)
sol.RoME <- RoME(base = sol.base$base, hierarchy = hierarchy, loss.diff = DiffLoss(type = 'LS'))
res <- sol.RoME$RoME - data.test
rmse <- RMSE(res, range.term = list(1:12),range.m = list(1,2:4,5:13))
```

## License

This package is free and open source software, licensed under GPL-3.
