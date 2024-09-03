# spatial-flu-forecasting

This code was used for the analyses in the paper "Does Spatial Information Improve Flu Forecasting?"

The file final_models.csv contains a list of models to be run with their parameters.

The model_pipeline.R file, which can be run from the commmand line, takes arguments for the file name and number of cores to use for parallel processing. Optionally, if you wish to run a subset of models from the CSV file, it takes arguments for start row and end row. If these are omitted, it will run all lines in the CSV file.

To run with 6 cores, for lines 1-4:
Rscript model_pipeline.R "final_models.csv" 6 1 4

This file will run the specified models and output the results in a CSV file with the date and time of the run, using the name of the input file.

**Data files:**

states_mat_adj.csv: contains a matrix indicating which states are neighbors

flu_data_all.csv: contains ILI data at the state level for October 2010 - May 2019; ILI data comes from CDC ILINet, accessed via the Delphi EpiData API

**Function files**

calc_fns.R: Contains helper functions for manipulating dataframes and implementing analysis. Called by model_pipeline.R 

quantile_tracker.R: Contains functions to implement quantile tracker to compute interval forecasts, as in Anastasios N. Angelopoulos, Emmanuel J. Candès, and Ryan J. Tibshirani. Conformal PID Control for Time Series Prediction. Advances in Neural Information Processing Systems, 36, 7 2023. 
Called by model_pipeline.R

do_model_fit_fn.R: Contains main function to perform model fitting, and calls functions from reg_fns.R. File also contains helper functions for arranging data. Called by functions in calc_fns.R 

reg_fns.R: Contains functions to fit last-value-carried-forward, quantile, Poisson, and linear regression models. Called by fit_models function in do_model_fit_fn.R

interval_fns.R: Contains functions to compute interval score and weighted interval score. Called by fit_models function in do_model_fit_fn.R





