Main Execution Scripts

The folder includes four main driver scripts, each corresponding to a specific empirical or simulation study in the manuscript:
- main_simulation.R  
  This script runs the main simulation study across multiple scenarios, sample sizes, and
  censoring rates. It implements DAer, IPW, and full-data estimations; saves raw simulation
  outputs; summarizes performance metrics (bias, MSE, computation time); and generates
  the tables and figures reported in Section 3/Supplementary of the manuscript.

All scripts listed above should be executed from the project root directory to ensure that relative paths are resolved correctly.

----------------------------------------------------------------------------------------------------
The code folder contains five main subdirectories, which store supporting code modules, datasets, figures, and output results.
These subdirectories are described in detail below.

----------------------------------------------------------------------------------------------------

- datagen_module.R  
  This script generates covariates, true event times, and true expectile regression
  coefficients under various data-generating scenarios. It also provides functions to
  impose right-, left-, or interval-censoring, either adaptively or using fixed,
  scenario-specific parameters. The output is a list of datasets suitable for subsequent
  DAer or IPW analyses.

- estimation_module.R  
  This script provides functions to estimate expectile regression coefficients under
  different censoring mechanisms, including DAer (iterative imputation), IPW (inverse
  probability weighting), and full-data regression (ignoring censoring). It also includes
  auxiliary functions for initial coefficient estimation and censored outcome imputation.

- plot_module.R  
  This script provides functions for loading and summarizing raw simulation results,
  computing MSE and average computation time, and generating visualizations, including
  stacked MSE plots, survival curves, and coefficient plots with confidence intervals.

- sim_module.R  
  This script contains core functions for running simulation studies and aggregating
  results across repeated replications.
