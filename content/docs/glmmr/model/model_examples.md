---
title: Coding the model
weight: 1
---

# A `Model` object
The `glmmrBase` package defines the `Model` class, which is built on the [R6](https://r6.r-lib.org/articles/Introduction.html) class system. R6 classes are object orientated, which means that when you create a new object of an R6 class, it has a set of methods and objects that "belong" to it. In the case of the `Model` class, it generates and stores all the relevant matrices for the GLMM internally. One can then call functions that belong the object to access or use the stored data. We will demonstrate this below.

The arguments to create a new `Model` object are:
* `formula`: A fomrula specifying the GLMM, see [Formula specification](model_specification). This can be omitted and separate formula specified to the arguments to `covariance` and `mean` below.
* `covariance`: An optional list. The list can contain a formula, parameter values, and data. 
* `mean`: An optional list. The list can contain a formula, parameter values, and data.
* `var_par`: Optionally, the (starting) value of the scale parameter.
* `data`: A data frame with the data for the model. This can be omitted if data are provided separately to the `mean` and `covariance` arguments via their lists.
* `family`: An R family.
We discuss each of these below and then provide some examples.

## Formula
The formula for the model. The formula can either be a single formula including both fixed and random effects, or each component can be specified separately in the argument for the mean and covariance. Separate specifications may be desirable if one wants to quick change the mean or covariance separately in an existing `Model` object. See [Formula specification](model_specification) for more details.

## Family and link function
The package currently supports the following families and link functions, which are specified using the standard R families in the `stats` package. Link functions are specified in the same way, for example `binomial(link="logit")`.

| Family | Link functions          | code         | 
|--------|-------------------------|--------------|
| Gaussian | Identity, log         | `gaussian()` |
| Binomial | Logit, log, identity  | `binomial()` |
| Poisson  | Log, Identity         | `poisson()`  |
| Gamma    | Log, Inverse, Identity| `gamma()`    |
| Beta     | Logit                 | `Beta()`     |


The Beta family is provided by the package function `Beta()`, which generates a barebones list specifying the family and link. We use a mean-variance parameterisation of the Beta family. The likelihood is:

$$
f(y_i | \mu_i, \phi) = \frac{y_i^{\mu_i\phi - 1}(1-y_i)^{(1-\mu_i)\phi - 1}}{B(\mu_i\phi, (1-\mu_i)\phi)}
$$

where $B()$ is the Beta function, and we use logit link

$$
\log\left( \frac{\mu_i}{1-\mu_i} \right) = \mathbf{x}_i\beta + \mathbf{z}_i \mathbf{u}
$$

We similarly use a mean-variance parameterisation for the Gamma regression function:

$$
f(y_i | \mu_i, \nu) = \frac{1}{\Gamma(\nu)}\left( \frac{\nu y_i}{\mu_i} \right)^\nu \exp \left( -\frac{\nu y_i}{\mu_i} \right) \frac{1}{y}
$$

where we also provide logit, inverse, and identity link functions for the specification of $\mu_i$.

## Data
Data are required to build the matrices and perform the calculations. When we want to analyse a study design prior to data being collected we will have to create our own covariate data for this purpose. We provide several useful functions to create `data.frame`s with different structures to support generating dummy data, see [Generating data](../creating_data).

## Parameter values



# Examples

