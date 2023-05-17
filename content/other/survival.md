---
title: Bayesian survival analysis
weight: 1
---

# Introduction
Survival analysis is a common and widely used set of methods for analysing time to event data. Here, we describe implementing a Bayesian survival analysis including specification of multiple models and model checking. The models are implemented in Stan. The best R package to use Stan is `cmdstanr`. See [the Stan webpage](https://mc-stan.org/cmdstanr/) for information on installing the R package `cmdstanr`, it also has a built in function to install Stan for you as well. If you cannot install this package then `rstan` should work and is available on CRAN. [Watson et al (2021)](https://academic.oup.com/jrsssc/article/70/5/1164/7033963) provides an extended Bayesian Survival analysis example.

# Data
We assume we have patients $i=1,...,n$ and data $\mathbf{D} = [T_i, y_i]$ where $T_i$ is the time to the event occurring or the time at which the invidiual was right-censored, and $y_i$ is an indicator equal to one if the event occurred and zero if they were censored. We also have observed covariates for the patients $\mathbf{x}_i$. 

# Statistical methods
The model is parameterised by a set of parameters $\theta$, so the density function for $T_i$ is $f(T_i | \mathbf{x}_i, \theta) = f_i(t) = h_i(t)S_i(t)$ where $h_i(t)$ is the hazard function, and $S_i(t) = \exp(-\int_0^t h_i(u) du)$ is the survival function. For a proportional hazards model we have:

$$
h_i(t) = h_0(t)\exp(\mathbf{x}_i \beta)
$$

where $\beta$ are the parameters of the linear predictor. We assume here that censoring is not informative; but the model can be easily extended to include a mortality component.

## Log-likelihood
The log likelihood of our model is

$$
\log L(D | \theta) = \prod_{i=1}^n y_i \log(f_i(t)) + (1-y_i)\log(S_i(t))
$$

## Hazard Functions
A key part of Bayesian modelling is to critically assess and compare different possible models to make sure our model is appropriate. We will specify three hazard functions for comparison here, although there are of course many more. The three functions we consider are:
* Exponential: $h_0(t) = \lambda_1$.
* Weibull: $h_0(t) = \lambda_2 t^{\lambda_2 - 1}$ where $\lambda_2$ is the shape parameter.
* Semi-parametric: There are lots of ways of implementing splines. Following Watson et al (2021) we will use $h_0 = \sum_{b=1}^B \lambda_{3,b} M_b(t; d=3, k)$ where $M_b(t; d = 3, k)$ is the $b$th basis term of a degree 3 M-spline with knot locations $k$, which we set at the upper and low tertiles, and $\lambda_{3,b}$ are M-spline coefficients. We include an intercept and constrain the coefficients to a simplex to ensure identifiability of the coefficients and the intercept in the linear predictor.

## Prior distributions
My preferred approach is to use weakly informative prior distributions. These priors constrain the parameter to a plausible range but do not provide information about the location of the parameter within this range. The parameters $\beta$ are the log hazard ratios, which we specify as $\beta \sim N(0,1)$. The parameters $\lambda_r$ are strictly positive, so we give them half-normal priors $\lambda_r \sim N(0,1)[0,\infty)$.

## Model comparisons
We want to choose the model that best fits our data, without overfitting. We use two methods for comparing the models: graphical posterior predictive model checks, and model fit criteria statistics. For the posterior predictive model checks we generate replicate data sets from the posterior predictive distribution of the model parameters and graphically compare these dataset to the actual data. A well fitting model should be able to produce comparable data sets. We use two statistics to compare the models. First, the survival function can be compared to the Kaplan-Meyer survival curve from the actual data. The posterior predictive distribution of the survival function at $t$ is:

$$
p(S_i(t)|\mathbf{D}) = Pr(T_i < t | \mathbf{D}) = \int Pr(T_i < t | \mathbf{D})p(\theta|D) d\theta
$$

and we can also compare the density of the event times

$$
p(T_i|\mathbf{D},\mathbf{x}_i) = \int f(T_i|\theta,\mathbf{x}_i)p(\theta|D) d\theta
$$

We use two statistics to compare the models: the WAIC and the Leave-one-out cross validation criterion (LOO-CV).

# Survival probabilities by age and sex
One of the motivations for this page is an analysis that aims to predict age and sex expected time to death and compare it to a national average. If we aim to predict the time to death for covariate values $\mathbf{x}_i = x$, and $u(x)$ is the national average time to death for that group, then the posterior predictive density of time to death is:

$$
p(T_i|\mathbf{D},\mathbf{x}_i = x) = \int f(T_i|\theta,\mathbf{x}_i = x)p(\theta|D) d\theta
$$

from which we can derive expected values $E(T_i|\mathbf{D},\mathbf{x}_i = x)$ and the probability

$$
Pr(T_i < u(x) | \mathbf{D},\mathbf{x}_i = x)
$$

# Code and implementation
## Stan model

Below we give the Stan model for the Weibull hazard function. I've also written it using the most up-to-date Stan programming syntax. The R package rstan is a few versions behind the main Stan program. The biggest difference if using R stan is that we will have to change the `array[N] real y` to `vector[N] y` or `real y[N]` and similarly with the other arrays.

```
data {
  int N; // number of observations in the data
  int Npred; // number of values to predict at
  int P; // number of covariates
  array[N] real time; // time to event
  array[N] int y; // indicator for death or censoring
  matrix[N,P] x; //covariates
  matrix[Npred,P] xpred; //covariates for new predictions
  array[2] int startend; // start and end values for time range to look at
}
parameters{
  vector[P] beta;
  real mu;
  real<lower=0> lambda;
}
transformed data {
  int ntime = startend[2] - startend[1]; // number of time points
}
model{
  //priors
  mu ~ normal(0,5);
  beta ~ normal(0,1);
  lambda ~ normal(0,1);
  
  //model
  for(n in 1:N){
    if(y[n]==0){
      target += weibull_lccdf(time[n] | lambda,exp(-(mu + x[n]*beta)/lambda));
    }
    if(y[n]==1){
      target += weibull_lpdf(time[n] | lambda,exp(-(mu + x[n]*beta)/lambda));
    }
  }
}
generated quantities{
  vector[N] log_lik;
  vector[N] pred_times;
  vector[Npred] pred_times_new;
  array[Npred,ntime] int surv_prob;
  
  for(n in 1:N){
    log_lik[n] = (1-y[n])*(weibull_lccdf(time[n] | lambda, exp(-(mu + x[n]*beta)/lambda)));
    log_lik[n] += y[n]*(weibull_lpdf(time[n] | lambda, exp(-(mu + x[n]*beta)/lambda)));
  }
  for(n in 1:N){
    pred_times[n] = weibull_rng(lambda,exp(-(mu + x[n]*beta)/lambda));
    for(i in startend[1]:startend[2]){
      surv_prob[n, i - startend[1] + 1] = pred_times[n] > i;
    }
  }
  for(n in 1:Npred){
    pred_times_new[n] = weibull_rng(lambda,exp(-(mu + xpred[n]*beta)/lambda));
  }
}
```

The model should hopefully be self-explanatory: the target is incremented by either the survival function or log-density depending on whether the individual was censored. The generated quantities block includes a calculation for the log-likelihood (`log_lik`). It also includes a vector for sampling from the posterior predictive distribution of surivival times for the data (`pred_times`), and for a matrix of pre-specified covariate values (`pred_times_new`), if needed. The final part is a 2d-array for estimating the survival function for each individual, which just includes an indicator for whether the survival time is greater than each value between the start and end. 

For the exponential hazard function it is just a case of changing the word `weibull` in the above program to `exponential`. For the splines, we just have to make some custom hazard and survival function log probability density functions. At the top of the program we can add a `functions` block. We also need to add the parameters of the spline as a simplex and include the spline data. The Stan program is then
```
functions {
  real ms_surv_lpdf(real T, real ispline, real eta){
    real lprob = -ispline*exp(eta);
    return(lprob);
  }
  real ms_h_lpdf(real T, real mspline, real eta){
    real lprob = log(mspline) + eta;
    return(lprob);
  }
}
data {
  int N; // number of observations in the data
  int Npred; // number of values to predict at
  int P; // number of covariates
  int S; // number of splines
  array[N] real time; // time to event
  array[N] int y; // indicator for death or censoring
  matrix[N,P] x; //covariates
  matrix[Npred,P] xpred; //covariates for new predictions
  array[2] int startend; // start and end values for time range to look at
  matrix[N,S] msp;
  matrix[N,S] isp;
}
parameters{
  vector[P] beta;
  real mu;
  simplex[S] coef_s; // spline parameters
}
transformed data {
  int ntime = startend[2] - startend[1] + 1; // number of time points
}
transformed parameters {
  vector[N] msp1; 
  vector[N] isp1;
  vector[N] eta;

  //linear predictors
  for(n in 1:N){
    eta[n] = mu + x[n]*beta;
    msp1[n] = msp[n] * coef_s;
    isp1[n] = isp[n] * coef_s;
  }
}
model{
  //priors
  mu ~ normal(0,5);
  beta ~ normal(0,1);
  coef_s ~ normal(0,1);

  //model
  for(n in 1:N){
    if(y[n]==0){
      target += ms_surv_lpdf(time[n] | isp1[n], eta[n]);
    }
    if(y[n]==1){
      target += ms_surv_lpdf(time[n] | isp1[n], eta[n]);
      target += ms_h_lpdf(time[n] | msp1[n], eta[n]); 
    }
  }
}
generated quantities{
  vector[N] log_lik;
  vector[N] runif; // vector of random numbers
  vector[N] pred_times;
  vector[Npred] pred_times_new;
  array[N,ntime] int surv_prob;
  
  for(n in 1:N){
    log_lik[n] = (1-y[n])*ms_surv_lpdf(time[n] | isp1[n], eta[n]);
    log_lik[n] += y[n]*(ms_surv_lpdf(time[n] | isp1[n], eta[n]) + ms_h_lpdf(time[n] | msp1[n], eta[n]));
  }
  for(n in 1:N){
    runif[n] = uniform_rng(0,1);
    pred_times[n] = -log(1-runif[n])/exp(mu + x[n]*beta);
    for(i in startend[1]:startend[2]){
      surv_prob[n, i - startend[1] + 1] = pred_times[n] > i;
    }
  }
  for(n in 1:Npred){
    pred_times_new[n] = -log(1-runif[n])/exp(mu + xpred[n]*beta);;
  }
}
```

## Estimating the model in R
To run these models in R, I would typically use cmdstanr, so I will provide instructions for that here. I will update to provide rstan instructions if needed. You will also need the package `loo`, which calculates the model checking statistics. The Stan model should be saved as a `.stan` file somewhere -- in RStudio you can create a new stan file under File -> New File -> Stan File. I will assume here it is saved at `C:/model.stan`. 

We first compile the model
```
require(cmdstanr)
model <- cmdstan_model("C:/model.stan")
```
We will need to put all the data together that appears in the `data` block for the program above as a list. For example,
```
dat <- list(
  N = nrow(df),
  Npred = nrow(Xpred),
  P = ncol(X),
  ...
)
```
For the splines, we will need to calculate the basis functions, which you can do with the `splines2` package. To use the upper and lower tertiles of the data, if the time to event data is a variable named `time` then we can simply use:
```
require(splines2)
knots <- quantile(time,c(0.33,0.67))
ms <- mSpline(time, knots=knots, degree=3, intercept = T)
is <- iSpline(time, knots=knots, degree = 3, intercept = T)
```
And then `ms` and `is` will be included in the data list e.g. `isp = is` and `msp = ms`. The `is` component is the basis functions for the integrated splines - the survival function is the integral of the hazard function, so we also need to integrate the splines.

We also need to include a matrix with all the values we want to predict the time to event for (i.e. the age and sex values), `xpred`. This should be the same structure as `x`.

To draw samples from the model we can then do:
```
fit.exp.1 <- mod.exp.1$sample(
  data = dat,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4,
  refresh = 50
)
```
changing the options as needed. If there are problems with the model fit, then Stan will print warnings to the console after the fitting is complete, such as `Warning: 4 transitions ended with a divergence.` If any warnings or errors appear then we will have to refine the model a bit. If the model is also taking a very long time then we can look at implementing other changes to try to speed it up.

## Model comparison
### Posterior predictive model checks
First, we can compare the density of the time to event outcome. We will pick a random sample of 100 draws to compare to the actual data. We first extract the predicted outcomes and take a random subset:
```
tpred <- fit$draws("pred_times",format="draws_matrix")
tpred <- tpred[,sort(sample(1:ncol(tpred),100,replace=FALSE))]
tpreddata <- reshape2::melt(tpred)
```
Assuming the actual time to event data are in a data frame called `df` and a column called `time`, then we can plot the actual data density versus the predicted data:
```
require(ggplot2)
ggplot()+
  stat_density(data=tpreddata,aes(x=value,group=Var2),position="identity",geom="line",color="blue",alpha=0.1)+
  geom_density(data=df,aes(x=time),size=0.5,color="red")
```
adding any extra options to make it pretty. This should be repeated for each model.

The survival curve is much the same. We just need to produce a data frame that has the predicted probabilities for each draw. In the model above we have included a matrix in the generated quantities `surv_prob` that generates predicted probabilities for each observation that the survival time will be longer than $t$. We can extract and plot this to generate a predicted survival function. For example,
```
probs <- fit.exp.1$draws("surv_prob",format="matrix")
iters <- sample(1:4000,100) 
dfp <- data.frame(t = rep(1:19,100), iter = rep(1:100,each=19), prob = NA) 
for(i in 1:100){ 
  for(j in 1:19){ 
    dfp[dfp$iter == i & dfp$t == j, 'prob'] <- mean(probs[iters[i], (1+1203*(j-1)):(1203*j)]) 
   } 
}
```
Then we could add a sample of curves to a Kaplan-Meier plot, for example
```
km1 + geom_line(data=dfp,aes(x=t,y=prob,group=iter),colour="blue",alpha=0.2)
```

### Model fit statistics
To get the WAIC we can use:
```
loglik <- fit$draws("log_lik")
w1 <- loo::waic(loglik)
```
and for the LOO-CV we have
```
loo <- fit$loo()
```
The lower the number, the better the fit!

## Required probabilities
Here we just need to look at the probabilities for each of the predicted values that it is less than the national average, or other relevant comparator. We will do this with our best fitting model. The predicted values are in `pred_times_new`, and the ordering of the columns will be the same as the rows of `xpred`. To extract the data we can use `out <- fit$draws("pred_times_new",format="draws_matrix")`. The mean values will be `colMeans(out)` and the 95\% credible intervals `apply(out,2,function(x)quantile(x,0.025))` for the lower interval and similarly for the upper interval. To get the probability the $i$th value is below a threshold $a$, we can do `length(out[out[,i]<a,i])/nrow(out)`. 
