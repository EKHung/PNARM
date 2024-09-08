# PNARM

The repository is organised into the following `R` scripts: 
- `Dissertation.R`: contains the code for generating the results in this dissertation by sourcing the functions in the scripts below.
- `Model_fitting.R`: contains the MCMC algorithms for fitting PNARM models and a function for fitting Poisson regression models.
- `Analysis.R`: contains functions for computing posterior distributions, least-squares partitions, assessing model fit, etc.
- `Stacking.R`: contains functions related to a chain-stacking procedure similar to [1].
- `Simulations_and_figures.R`: contains functions for simulating data from a Poisson vector autoregression model given cluster memberships and coefficients, and for creating figures for easy visualisation of model outputs.
- `Helper_functions.R` contains helper functions used in the functions of the other scripts.
- `gagnar.R`: contains an MCMC algorithm for fitting GAGNAR models [2]. The code is almost entirely the same as from the supplementary material of their paper.
- `COVID-19.csv`: is the raw data set of the COVID-19 cases from the Republic of Ireland.
- `eco_hubs.Rdata`: contains the igraph object for the economic hubs network. This was taken from the GitHub repository associated with [3]. 

References:

[1]  Yao, Y., Vehtari, A., and Gelman, A. (2022). Stacking for non-mixing Bayesian computations: The curse and blessing of multimodal posteriors. Journal of Machine Learning Research, 23(79):1–45.

[2]  Ren, Y., Zhu, X., Lu, X., and Hu, G. (2024). Graphical assistant grouped network autoregression model: A Bayesian nonparametric recourse. Journal of Business and Economic Statistics, 42:49–63.

[3] Armbruster, S., Reinert, G. Network-based time series modeling for COVID-19 incidence in the Republic of Ireland. Appl Netw Sci 9, 23 (2024). https://doi.org/10.1007/s41109-024-00634-2
