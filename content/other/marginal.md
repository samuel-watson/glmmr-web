---
title: Marginal standardisation with mixed models
weight: 1
---

# Marginal standardisation with mixed models
Note: this is preliminary work and may contain errors.

The basic presentation of marginal standardisation in the epidemiological literature is of averaging the effect of an exposure over a discrete confounder. For example, if $Z$ is our confounder $Pr(outcome \vert exposure) = \sum_{Z}Pr(outcome \vert exposure, Z=z)Pr(Z=z)$. The idea is to obtain population averaged effects. Marginal standardisation is also used in the trials literature, often to obtain relative risks or risk differences from models like a Binomial-logistic that does not naturally have that interpretation and where there may be covariates. We can calculated fitted probabilities for trial participants under both treatment and control conditions and then take their ratio or difference. Things become more complex though when we have clustering or correlation between the observations as we need to also average over the random effects in a mixed model. This document is a proposal for how one could do that. I'll provide a statistical explanation of the idea then some code to implement it in R.

## Statistical model
### Relative risk
To generalise the description above we assume that the covariates are stochastic and drawn from some set of possible value $X \in \mathcal{X}$. We assume a [generalise linear mixed model](../docs/glmm/_index.md), which have the random effects $\mathbf{u} \sim N(0,D)$. We use $d$ to represent our treatment effect indicator. The linear predictor for an individual $i=1,...,n$ is then:
$$
\eta_i = d_i \alpha + x_i \gamma + z_i \mathbf{u}
$$
We also represent all the model parameters as $\Theta$ and notate $E(y|d,\mathbf{u},X,\Theta) = \mu_d(\mathbf{u},X,\Theta) = h^{-1}(\eta_i)$ where $h^{-1}$ is the link function.

We focus on the relative risk here, but the analysis carries over to risk difference as well. The marginal standardisation relative risk is
$$
RR = \frac{E[y|d=1,\Theta]}{E[y|d=0,\Theta]} = \frac{E_X[E_{\mathbf{u}}[\mu_1(\mathbf{u},X,\Theta)]]}{E_X[E_{\mathbf{u}}[\mu_0(\mathbf{u},X,\Theta)]]}
$$
and equivalently the log relative risk is:
$$
\log(RR) = \log(E_X[E_{\mathbf{u}}[\mu_1(\mathbf{u},X,\Theta)]]) - \log(E_X[E_{\mathbf{u}}[\mu_0(\mathbf{u},X,\Theta)]])
$$

The main issue we have here is the random effects, since:
$$
E_{\mathbf{u}}[\mu_d(\mathbf{u},X,\Theta)] = \int \mu_d(\mathbf{u},X,\Theta) f_{\mathbf{u}}(\mathbf{u}) d\mathbf{u}
$$
and the integral does not generally have a closed form solution. There are different strategies for handling this integral and expectation. One option is to approximate the integral using, for example, a Laplace approximation or Gaussian quadrature. An alternative option, and one I'll pursue here, is to instead treat the random effects as missing data. [McCulloch (1997)](https://doi.org/10.1080/01621459.1997.10473613) proposes a strategy for full maximum likelihood estimation that involves sampling the random effects using Markov Chain Monte Carlo (MCMC). MCMC can be used to obtaining samples of $\mathbf{u}$ since $f_{\mathbf{u}\vert y}(\mathbf{u}\vert y, X,\Theta) \propto f_{y \vert \mathbf{u}}(y \vert \mathbf{u}, X,\Theta)f_{\mathbf{u}}(\mathbf{u})$. So if we draw $M$ samples $\hat{\mathbf{u}}^{m}$ we can use the approximation:
$$
E_{\mathbf{u}}[\mu_d(\mathbf{u},X,\Theta)] \approx \frac{1}{M}\sum_{m=1}^M \mu_d(\hat{\mathbf{u}}^m,X,\Theta)
$$
To average over $X$ to estimate the outermost expectation, we just plug in the whole matrix $X$ and average the values. So, our estimator is:
$$
\hat{RR} = \frac{\frac{1}{NM}\sum_{i=1}^N \sum_{m=1}^M \mu_1(\hat{\mathbf{u}}^m,x_i,\Theta)}{\frac{1}{NM}\sum_{i=1}^N \sum_{m=1}^M \mu_0(\hat{\mathbf{u}}^m,x_i,\Theta)}
$$

### Standard errors
We will focus on the log relative risk and use the delta method to obtain our standard errors. To do this we need the gradient of the log relative risk with respect to the model parameters $\beta = [\alpha, \gamma]$:
$$
\nabla_\beta \log(RR) = \nabla_\beta \log(E_X[E_{\mathbf{u}}[\mu_1(\mathbf{u},X,\Theta)]]) - \nabla_\beta \log(E_X[E_{\mathbf{u}}[\mu_0(\mathbf{u},X,\Theta)]])
$$
$$
 = \frac{\nabla_\beta E_X[E_{\mathbf{u}}[\mu_1(\mathbf{u},X,\Theta)]]}{E_X[E_{\mathbf{u}}[\mu_1(\mathbf{u},X,\Theta)]]} - \frac{\nabla_\beta E_X[E_{\mathbf{u}}[\mu_0(\mathbf{u},X,\Theta)]]}{E_X[E_{\mathbf{u}}[\mu_0(\mathbf{u},X,\Theta)]]}
$$
$$
= \frac{ E_X[E_{\mathbf{u}}[\nabla_\beta \mu_1(\mathbf{u},X,\Theta)]]}{E_X[E_{\mathbf{u}}[\mu_1(\mathbf{u},X,\Theta)]]} - \frac{ E_X[E_{\mathbf{u}}[\nabla_\beta \mu_0(\mathbf{u},X,\Theta)]]}{E_X[E_{\mathbf{u}}[\mu_0(\mathbf{u},X,\Theta)]]}
$$
where we make sufficient assumptions to take the derivative inside the integral ($\mu$ is smooth etc.). 

The denominators of the fractions we can obtain using the same estimator as for the log relative risk. For the numerator we need the gradient of the mean. We will give results for a Binomial-logistic model here with link function $\mu(x) = exp(x)/(1+exp(x))$. In this case we have:
$$
E_X[E_{\mathbf{u}}[\nabla_\beta \mu_d(\mathbf{u},X,\Theta)]]= \frac{1}{N}E_{\mathbf{u}}[ X^{*T} \mu_d(\mathbf{u},X,\Theta) (1-\mu_d(\mathbf{u},X,\Theta)) ]
$$
$$
 \approx \frac{1}{NM} \sum_{m=1}^M X^{*T} \mu_d(\hat{\mathbf{u}}^m,X,\Theta)
$$
where $X^{*}$ just means the matrix $X$ with the treatment indicator added as an additional column.

The standard error for the log relative risk is then given by:
$$
SE_{logRR} = \sqrt{\nabla_\beta^T M_\beta^{-1} \nabla_\beta}
$$
where $\nabla_\beta = \nabla_\beta \log(RR)$ and $M_\beta$ is the informaion matrix of the model parameters. For the GLS estimator of the GLMM this is $M = X^{*T} \Sigma^{-1} X$.

## Coding the marignal standardisation
I will use our statistical package `glmmrBase` as it has most of the relevant functions. You will need the development version and to install it from github, see [here](../docs/glmmr/_index.md) for more details. I would also strongly recommend using `cmdstan` as the MCMC sampler. While there is a sampler in the package that you can use if you can't install Stan, Stan is just better! See [the Stan webpage](https://mc-stan.org/cmdstanr/) for information on installing the R package `cmdstanr`, it also has a built in function to install Stan for you as well.

The links on the left provide information on model specification and other features of `glmmrBase` so I won't provide too much information here. I will specify a stepped-wedge cluster randomised trial design with six clusters and seven time periods with ten individuals per cluster-period cross-sectionally sampled. I will also specify for this example an autoregressive/exponential decay random effect; if $\psi_{kt}$ is the random effect term for cluster $k$ at time period $t$ then I use $Cov(\psi_{kt},\psi_{kt'})= \tau^2 \lambda^{\vert t - t' \vert}$. I'll set $\tau = 0.4$ and $\lambda = 0.8$. I'll randomly generate time period fixed effect parameters from a $N(0,0.1^2)$ distribution and choose a treatment effect parameter $\alpha = 0.3$. We first generate the data frame:
```
require(glmmrBase)
df <- nelder(~(cl(6)*t(7))>ind(10))
df$int <- 0
df[df$cl > df$t, 'int'] <- 1
```
and then the model is specified as:
```
mod <- Model$new(
  formula = ~ factor(t) + int -1 + (1|gr(cl)*ar1(t)),
  mean = list(parameters = c(rnorm(7,0,0.1),0.3)),
  covariance = list(parameters = c(0.4,0.8)),
  family = binomial(),
  data = df
)
```
In `glmmrBase` the model is an object of class `Model` and contains all the relevant functions for that specific GLMM. For example, to simulate outcome data we use:
```
y <- mod$sim_data()
```
I'll also extract the matrix $X$ and create two copies that either have all ones or all zeros in the treatment effect indicator:
```
X1 <- X0 <- mod$mean$X
X0[,8] <- 0
X1[,8] <- 1
```
Now, to draw the random effects using MCMC we just do
```
u <- mod$mcmc_sample(y$y,usestan = TRUE)
```
You can change `usestan=TRUE` to `usestan=FALSE` if you can't use Stan. You can also change the MCMC options by changing values in the list `mod$mcmc_options`.

We now have all the data we need to estimate the relative risk and its standard error.The way I've done this is to create vectors to hold the mean fitted value for $\mu_d$ and another vector to hold the mean values for $X^{*T} \mu_d(\hat{\mathbf{u}}^m,X,\Theta)$. We fill up these vectors and then divide by $M$:
```
y1 <- y0 <- rep(0,length(y))
m1 <- m0 <- rep(0,ncol(X1))

for(i in 1:ncol(u)){
  mu1 <- mod$fitted(type="response",X1,u[,i])
  mu0 <- mod$fitted(type="response",X0,u[,i])
  y1 <- y1 + mu1
  y0 <- y0 + mu0
  m1 <- m1 + t(X1)%*%(mu1*(1-mu1))/nrow(X1)
  m0 <- m0 + t(X0)%*%(mu0*(1-mu0))/nrow(X0)
}

y1 <- y1/ncol(u)
y0 <- y0/ncol(u)
m1 <- m1/ncol(u)
m0 <- m0/ncol(u)
```
The estimate of the relative risk is then:
```
mean(y1)/mean(y0)
> [1] 1.157976
```
The gradient is 
```
grad1 <- m1/mean(y1) - m0/mean(y0)
```
The information matrix for the model can be obtained with `mod$information_matrix()`, so the standard error is:
```
sqrt(t(grad1)%*%solve(mod$information_matrix())%*%grad1)
>           [,1]
> [1,] 0.1572163
```
This would give a relative risk and 95\% confidence interval of 1.15 [0.85, 1.58].

Disclaimer: I don't know if this is the correct value or not! I haven't found anything to compare it against. So I will update this page when we can generate some comparisons. 