# Non-Linear Time Series Modelling: A Comparative Study for Crude Oil Price Forecasting

**Environment**

*main.Rmd*: Main compiler, returns the modelling outcomes of the selected models. The script is run through the following helper functions:
- Kernels.R: Includes the setups of the squared exponential, Matern and Brownian motion kernel and functions to sample from those distributions in 2D and 3D.
- Data-processing: reads in the data and returns a cleaned, transformed version of the data from January 2003 until December 2022. Moreover, the train-test split function is written in this file.

*Neural Networks.ipynb*: Construction of feed-forward neural networks and LSTM models


*other*: Other, unrelated codes are:
- 3D samples.R: plot draws from a multivariate normal distribution in 3D, mainly used as debug file for Kernels.R 
- GPR.R: provides multiple kernel functions, mainly used as a debug file for MCMC.stan

**Data**

Currently contains the Crude Oil Prices: West Texas Intermediate (WTI) - Cushing, Oklahoma MCOILWTICO series from the Federal Reserve Economic Research (FRED) database. The original time series starts in January 1986 and can be derived from https://fred.stlouisfed.org/series/MCOILWTICO. The explanatory variables are the consumer price index (FRED), Kilian index (FED Dallas), global production (JODI-Oil) and stock change (JODI-Oil).


