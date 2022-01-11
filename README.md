# growth-rate-estim
Reproducible model for: "Bayesian Estimation of real-time Epidemic Growth Rates using Gaussian Processes: local dynamics of SARS-CoV-2 in England” (2022).

Link: https://doi.org/10.1101/2022.01.01.21268131.

Language: `R`.

To see an example to run the model, go to `main.R` (it runs the model for the validation section).

Main functions:
* `setParametersFn()`: sets up the main parameters of the model (e.g. priors).
* `runModelGrowthRate()`: runs the model for the data in `countTable`.

Requires the following packages:
* `data.table`.
* `INLA`.
* `ggplot2`.
