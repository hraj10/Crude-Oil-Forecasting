# ANOVA-kernel

**Environment**

*main.Rmd*: Main compiler, returns the modelling outcomes of the selected models. The script is run through the following helper functions:
- Kernels.R: Includes the setups of the squared exponential, Matern and Brownian motion kernel and functions to sample from those distributions in 2D and 3D.
- Data-processing: reads in the data and returns a cleaned, transformed version of the data from January 1990 to the most recent month in 2023. Moreover, the train-test-val split is conducted in this file.

*other*: Other, currently drafted and unused, codes are:
- Rasmussen sample code.R
- ST499 code.R
- attempt two.R

**Data**

Currently contains the Crude Oil Prices: West Texas Intermediate (WTI) - Cushing, Oklahoma MCOILWTICO series from the Federal Reserve Economic Research (FRED) database. The original time series starts in January 1986 and can be derived from https://fred.stlouisfed.org/series/MCOILWTICO.
