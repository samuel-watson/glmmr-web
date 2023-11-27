---
title: Model Class Functions
weight: 1
---

# `Model` class functions
A full list of functions to be added here.

| Method | | Description |
|--------|-|-------------|
|`covariance` | `D` | The matrix D, the covariance matrix of the random effects
|             | `Z` | The matrix Z, the design matrix of the random effects
|             | `chol_D()` | Generates the Cholesky decomposition of D 
|             | `log_likelihood` | For a given vector $u$ calculates the multivariate Gaussian log-likelihood with zero mean and covariance D |
|             | `simulate_re()`  | Simulates a vector $u$
|             | `sparse()`       | Choose whether to use sparse matrix methods
|             | `parameters`     | The parameters $\theta$
|             | `formula`        | The random effects formula used to create D
|             | `update_parameters()` | Updates $\theta$ and related matrices
|             | `nngp()` | Sets or returns options for the nearest neighbour Gaussian process approximation
|             | `hsgp()` | Sets or returns options for the Hilbert space Gaussian process approximation
| `mean`      | `X`  | The matrix X, the design matrix of fixed effects 
|             | `parameters` | The parameters $\beta$
|             | `offset`     | The optional model offset
|             | `formula`    | The fixed effects formula used to create X
|             | `linear_predictor()` | Generates the linear predictor 
|             | `update_parameters()` | Updates $\beta$ and related matrices
| `family` | | An R family object
| `var_par` | | The scale parameter
| `fitted()` | | Full linear predictor including random effects (either simulated or provided)
| `predict()` | | Predictions from the model at new data values
| `sim_data()` | | Simulates data from the model
| `Sigma()` | | Generates $\Sigma$ (or an approximation)
| `information_matrix()` | | Generates the information matrix
| `sandwich()` | | Generates the robust sandwich matrix
| `kenward_roger()` | | Small sample bias-corrected variance matrix of $\hat{\beta}$
| `partial_sigma()` | | Generates matrices $\partial \Sigma / \partial \theta$ and $\partial^2 \Sigma / \partial \theta_i \partial \theta_j$
| `use_attenuation()` | | Option for improving the approximation of $\Sigma$
| `power()` | | Estimates the power
| `MCML()`  | | MCMC Maximum likelihood model fitting
| `LA()`    | | Laplace approximation model fitting
| `mcmc_sample()` | | Sample $u$ using MCMC
| `w_matrix()` | | Returns $diag(W)$
| `dh_deta()`  | | Returns $\partial h^{-1}(\eta) / \partial \eta$
| `log_gradient()` | | Returns the gradient of the log likelihood with respect to either the random effects or $\beta$
| `marginal()`  | | Calculates marginal effects with different types of conditioning or averaging