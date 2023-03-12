---
weight: 1
bookFlatSection: true
title: "What is a generalised linear mixed model?"
---

# What is a generalised linear mixed model?
A generalised linear mixed model (GLMM) is a flexible statistical model that allows for correlation between observations through the incoporation of "random effects" into the model. There may be different reasons for including the random effects in a statistical model. In some cases, we are interested in estimating the effect of a covariate on an outcome. Ignoring the correlation between observations can result in standard errors that are too small as we do not take into account the loss of information due to the correlation. In other cases, the random effects are of interest in and of themselves. These random effects might represent the latent, unobserved variation in an outcome of interest and can identify where average outcomes may be higher or lower than expected given the observed covariates. For example, GLMMs are used in geospatial modelling work to estimate where risk of a disease may be higher or lower than expected, and to quantify the uncertainty surrounding these estimates.

# A statistical model
A general specification of a GLMM is:
$$
y \sim F(\mu,\phi) 
$$
$$
\mu = h^{-1}(X\beta + Z\mathbf{u}) \\
$$
$$
\mathbf{u} \sim N(0,D)
$$
where $y$ is the vector of outcome data, $F$ is a distribution such as Gaussian or Binomial, with mean $\mu$, and $\phi$ is a (optional) scale parameter. The mean is specified using a link function (such as logit or log) $h^{-1}$, a matrix of covaraite data $X$ with parameters $\beta$, and the random effects, which are defined using a "design matrix" $Z$ and the random effects themselves $\mathbf{u}$. Finally, the random effects are assumed to be distributed according to a multivariate Gaussian distribution with mean zero and covariance matrix $D$.

It is the random effects that differentiate these models from standard generalised linear models. 