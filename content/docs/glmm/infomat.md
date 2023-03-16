---
title: Coding the model
weight: 1
---

# Completing the model specification

## Supported Families
The package currently supports the following families and link functions, which are specified using the standard R families in the `stats` package.

| Family | Link functions          | code |
|--------|-------------------------|------|
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
