---
title: Optimal Examples
weight: 1
---

# Introduction
This web page accompanies the article "Optimal Study Designs for Cluster Randomised Trials: An Overview of Methods and Results" forthcoming in Statistical Methods in Medical Research. There is an arXiv version of the paper [here](https://arxiv.org/abs/2303.07953). There is more information about model specification and analysis on the other pages on this website, although not all of them are complete yet. The examples below cover combinatorial optimisation algorithms, the conic optimisation methods, and the [Girling algorithm](girlingalgo), which can all be run using the R package `glmmrOptim`. The package can be used to find optimal designs for any study design that can be represented by a generalised linear mixed model. Here, we focus specifically on cluster randomised trials, but other examples are presented elsewhere on this website.

# Data and model specification
The examples in the article all use the "stepped-wedge" design space, which includes $J$ cluster sequences, $J-1$ time periods and where the intervention roll out is staggered over time. The figure below (reproduced from Watson, Girling, and Hemming) shows this design space a "Design Space A".

![Design spaces](images/design_spaces.jpg)

We must first create a data frame that describes the design space or model of interest. We will use seven cluster sequences:
```
require(glmmrOptim)
df <- nelder(formula(~ (cl(7) * t(6)) > ind(10)))
df$int <- 0
df[df$t >= df$cl,'int'] <- 1
```
See [Data specification](../docs/glmmr/creating_data.md) for more information on creating data sets. We can then specify a generalised linear mixed model using the `glmmrBase` model class. For this example we shall use a linear mixed model with nested exchangeable covariance function (EXC2 in the paper), an ICC of 0.05 and a CAC of 0.8. Assuming an individual level variance of $\sigma^2 = 1$, these values give parameter values of $\tau^2 = 0.05*0.8/0.95$ (cluster-level variance) and $\omega^2 = 0.05*0.2/0.95$ (cluster-period level variance). We can then specify:
```
model <- Model$new(
  formula = ~ factor(t) + int - 1 + (1|gr(cl)) + (1|gr(cl,t)),
  covariance = list(parameters = c(0.05*0.8/0.95,0.05*0.2/0.95)),
  family = gaussian(),
  data = df
)
```
The default values is $\sigma^2 = 1$ and $\beta$ is initialised to random values, so these are not directly set here. 

# Clusters as experimental conditions
All of the `glmmrOptim` algorithms are accessed through a `DesignSpace` class, which holds any number of `Model` objects. When a new `DesignSpace` is created, we provide one or more `Model` objects and the identifiers of the experimental conditions (if this is not provided, it is assumed that each row/observation is an experimental condition). For these examples we will consider that the cluster sequences are the experimental condition
```
ds <- DesignSpace$new(model,experimental_condition = df$cl)
```

## Conic optimisation 
As the experimental conditions are uncorrelated, we can use the conic optimisation methods described by [Holland-Letz T, Dette H, and Pepelyshev (2011)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00757.x) and [Sagnol (2011)](https://linkinghub.elsevier.com/retrieve/pii/S0378375810005318). The function will automatically detect that the experimental conditions are uncorrelated and will default to using this method; to use the combinatorial algorithms instead we can set the argument `use_combin=TRUE`. See later sections for details on combinatorial algorithms. The function will also return the set of optimal apportionments using a range of methods, the value of `m` specifies the total number of clusters for these calculations.
```
outw <- ds$optimal(m = 10,C = list(c(rep(0,6),1)),verbose = TRUE)

> outw
$weights
[1] 0.2101400 0.1159425 0.1159449 0.1159452 0.1159449 0.1159425 0.2101400

$designs
$designs$Pukelsheim_1
[1] 2 1 1 2 1 1 2

$designs$Hamilton
[1] 2 1 1 2 1 1 2

$designs$Jefferson
[1] 2 1 1 1 1 1 3

$designs$Webster
[1] 2 1 1 1 1 1 3
```

# Observations as experimental conditions
We now set up a design space where each observation is an experimental condition:
```
ds <- DesignSpace$new(model)
```
The same function `optimal` is used to run the other algorithms. The argument `algo` specifies which algorithm to run, the options are:
 - `1` Local search
 - `2` Greedy search 
 - `3` Reverse greedy search
 - `"girling"` Girling's algorithm
For the combinatorial algorithms, we can also string together numbers to run one algorithm after another, such as `c(3,1)`. 

## Combinatorial algorithm

## Girling algorithm

# Non-linear models

# Small sample bias corrections

# Robust optimal designs