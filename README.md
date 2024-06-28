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
# load dara
data("Tourism")
data <- Tourism$data
size.n = ncol(data)
H1 <- Tourism$H1
diag1<-diag(1,size.n)
S <- as.matrix(rbind(H1,diag1))
weight = 1 / rowSums(S)

# split to train and test data set
data.all = as.matrix(data) %*% t(S)
size.T = nrow(data)
size.m = ncol(data)
order.train = 1:96
order.test = (1:12)+96
data.train = data.all[order.train, ]
data.test = data.all[order.test, ]

# forecat with RoME
sol.base = BaseAuto(x = data.train, method.forecast = 'ETS', term = 12)
sol.Wt = CovarianceW(data.train, method.forecast = 'ETS', is.rolling = FALSE)
sol.LAD <- RoME_fore(base = sol.base$base, hierarchy = S, loss.diff = DiffLoss(type = 'LAD'),
design.W = 'OLS', matrix.W = sol.Wt$W, lambda.shrink = sol.Wt$lambda)
residuals = sol.LAD$RoME - data.test
rmse <- RMSE(residuals, range.term = list(1:12), range.m = list(1,2:8,9:35,36:111), weight.m = weight)

```

## License

This package is free and open source software, licensed under GPL-3.
