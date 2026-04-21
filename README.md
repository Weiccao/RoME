We provide code to reproduce the results of the Numerical Experiments in our work Robust Optimal Reconciliation for Hierarchical Time Series Forecasting with M-estimation. Four experimental settings are considered:

* Experiment 1: Non-Gaussian Series
    This experiment investigates three different error distributions for a subset of bottom-level series.
* Experiment 2: Proportion of Bad Forecasts
    This experiment examines the impact of varying proportions of poor base forecasts on reconciliation performance in hierarchical time series (HTS).
* Experiment 3: Correlation Between Series
    This experiment studies the effect of different cross-series dependence structures on the reconciled results.
* Experiment 4: Hierarchical Complexity
    This experiment evaluates the impact of varying levels of hierarchical complexity on reconciliation performance.

All scripts should be executed from the project root directory to ensure that relative paths are correctly resolved.

⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻

Description of Main Simulation Code

* RoME_Basic.R
    Provides core functions for the RoME algorithm and base forecasting methods.
* RoME_Generation.R
    Contains functions for synthetic data generation.
* RoME_Process.R
    Implements the main simulation procedures to reproduce the results in Robust Optimal Reconciliation for Hierarchical Time Series Forecasting with M-estimation.
* RoME_Summary.R
    Provides functions to summarize simulation results, including AvlRelMSE and relative RMSE improvement.

⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻

Example Usage
```
# Generate HTS data
Generate_NonGaussianSeries(type.epsilon = c("Normal", "Mixture", "Tstudent"))
# Run simulation and perform RoME forecasting
Process_NonGaussianSeries(
  method.forecast = "ARIMA",
  design.W = c("OLS", "WLSv", "WLSs", "Sample", "Shrink")
)
# Export and summarize results
Summary_NonGaussianSeries(
  range.term = list(1, 1:6, 1:12),
  range.m = list(c(2, 5:7), c(3:4, 8:13), 1:13)
)
```
