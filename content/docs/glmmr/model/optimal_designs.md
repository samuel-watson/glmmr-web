---
title: Optimal Experimental Designs
weight: 1
---

# c-Optimal experimental designs
This page describes the combinatorial algorithms described in [Watson, Girling, and Hemming (2023)](https://arxiv.org/abs/2303.07953) and [Watson and Pan (2023)](https://arxiv.org/abs/2207.09183) used to identify c-Optimal experimental designs for experiments with correlated observations. We consider in this example the identification of optimal cluster randomised trial designs, as given in those articles, however more examples will be added later. [Watson, Girling, and Hemming (2023)](https://arxiv.org/abs/2303.07953) summarise and describe a range of other approaches for cluster randomised trials as well.

## What is c-optimality?
Put simply, a c-optimal experimental design is one that minimises the variance of a linear combination of model parameters. Typically, we are only interested in one parameter that represents a treatment effect of interest. As we are considering experiments with correlated observations, we will use a [generalised linear mixed model](../../glmm/_index.md) and consider the generalised least squares estimator of the model fixed effect parameters $\beta$. The variance matrix of the estimates $\hat{\beta}$ is $V = X^T \Sigma^{-1} X$, where $\Sigma^{-1}$ is the covariance matrix of the observations, and so the information matrix is $M = V^{-1}$. The c-optimality criterion for a design $d$ is
$$
g(d) = c^T M_d c
$$
where $c$ is a vector picking out the linear combination of parameters of interest $c^T\beta$.

## What is "a design"?