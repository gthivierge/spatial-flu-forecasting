# spatial-flu-forecasting

This code was used for the analyses in the paper "Does Spatial Information Improve Flu Forecasting?"

The file final_models.csv contains a list of models to be run with their parameters.

The model_pipeline.R file, which can be run from the commmand line, takes arguments for the file name and number of cores to use for parallel processing. Optionally, if you wish to run a subset of models from the CSV file, it takes arguments for start row and end row. If these are omitted, it will run all lines in the CSV file.

To run with 6 cores, for lines 1-4:
Rscript model_pipeline.R "final_models.csv" 6 1 4

This file will run the specified models and output the results in a CSV file with the date and time of the run, using the name of the input file.

Data files:

states_mat_adj.csv: contains a matrix indicating which states are neighbors
flu_data_all.csv: contains ILI data at the state level for October 2010 - May 2019; ILI data comes from CDC ILINet, accessed via the Delphi EpiData API


