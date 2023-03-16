---
title: Estimating power
weight: 1
---

# Estimating power

## Calculating the power
The power for a model parameter $\beta_1$ for a two-sided null hypothesis test at level $\alpha$ is given by
$$
Pow(\beta_1) = \Phi \left(\frac{\vert \beta_1 \vert}{SE_{\beta_1}} - \Phi^{-1}\left(1-\frac{\alpha}{2}\right)\right)
$$
where $\Phi$ is the Gaussian distribution function and $SE_{\beta_1}$ is the standard error. We can estimate the standard error of the [GLMM](../../glmm/_index.md) with the generalised least squares estimator by calculating the information matrix:
$$
M = (X^T \Sigma^{-1} X)^{-1}
$$
then the standard error for the $p$th parameter is $\sqrt{M_{p,p}}$. For a Gaussian model with identity link, $\Sigma$ can be calculated exactly as $\Sigma = \sigma^2I + ZDZ^T$, and we can then plug the values in we need to calculate the power. However, for non-Gaussian models, obtaining the covariance matrix is much harder. The package `glmmrBase` provides [approximations](../../glmm/infomat.md) in these cases, which we can use to quickly estimate the power.

## Example: Cluster randomised trial with a Solomon design
There are a multitude of packages available for R to estimate power for cluster randomised trials. These are shown in the table below. Many of these are relatively limited and provide functionality for specific models with explicit formulae. Because `glmmrBase` can handle a wide range of model specifications and generates the relevant matrices we can calculate the power easily for any type of trial we want under a range of functions. We just have to build the dataset representing the trial. Here, we give an example of a more complex trial design to illustrate.

### The Solomon Design
A Solomon design trial was developed in the 1940s for use in educational trials. One way of evaluating the performance of an educational intervention is to test students before and after using it. However, the act of completing a test on a particular topic may sensitise students to the material and itself act as an intervention. The Solomon design therefore includes four arms to enable us to account for this interaction. The different arms are shown in the Table. 

| Trial arm | Pre-test | Intervention | Post test |
|-----------|----------|--------------|-----------|
| A         | ✓        |  ✓          | ✓         |
| B         | ✓        |             | ✓         |
| C         |          | ✓           | ✓         |
| D         |          |             |  ✓         |

### Statistical Model
Two arms are intervention and two are control, with one of each receiving a pre-test, and all receiving the post-test. If we specify a model for individual $i$ in cluster $k$ at time $t$ where $t$ is an indicator for a post-test observation ($t=1$) or pre-test $t=0$;. We use $p_{k}$ to represent whether the cluster is in a pre-test arm ($p_k=1$) or not ($p_k=0$); and $d_k$ is an indicator for if the cluster is in an intervention arm ($d_k=1$) or control arm ($d_k=0$). Then the linear predictor, without random effects, is:
$$
\eta_{ikt} = \beta_0 + \beta_1 d_k + \beta_2 t + \beta_3 d_k t + \beta_4 p_k t + \beta_5 d_k p_k t
$$
The parameter of interest giving the treatment effect is $\beta_5$.

We will assume that the test produces a score as the outcome. Many analyses might use a linear model for these scores, however to illustrate a slightly more complex example, we will assume that the score is relatively small, say 10 items, and we assume the participant has an underlying probability of answering a question correctly that we wish to improve. So we specify a Binomial-logistic mixed model with the outcome $y_{ikt}$ being the total number of correct responses:
$$
y_{ikt} \sim \text{Binomial}(\mu,10)
$$
$$
\mu = \frac{1}{1+\exp(-(\eta_{ijk}+\alpha_{kt}))}
$$
where $\alpha_{kt}$ represents our random effect. As there are only two periods we will use a "nested exchangable" structure $\alpha_{kt} = \alpha_{1k} + \alpha_{2kt}$ with $\alpha_{1k} \sim N(0,\tau_1^2)$ and $\alpha_{2k} \sim N(0,\tau_2^2)$. 

### Parameter values
We must select parameter values for the model at which to calculate power. As the model is non-linear, the power depends on the choice of $\beta$ as well as the covaraince parameters. We will initially choose $\beta_0,...,\beta_4 = 0$ so that the mean score is 50\%. We will also choose $\tau_1 = 0.25$ and $\tau_2 = 0.1$, which are relatively high. Finally, we will look at values of $\beta_5$ from zero to two, which gives approximate changes to the mean proportion correctly answers from 0 percentage points (pp) to 40pp. 

In terms of sample size we will have 40 clusters allocated in equal proportion to each trial arm. We will also assume that the number of participants in each cluster is not fixed and is likely to be heterogeneous and betweem two and ten. For the first example, we will generate a random sample for the cluster sizes, but we can average the power over this sample size distribution by re-simulating.

### Setting up the data
We now have to create a data frame representing the trial described above. This is a little more complex than the standard blocked structure that we can create with the [nelder()](../../glmmr/creating_data.md) function, so we will start from scratch. First, we will generate a data frame describing the pre- and post-test variables and then expand each row to the desired number of clusters:
```
df <- expand.grid(pre = c(0,1), treat = c(0,1))
df$t <- 2
df1 <- expand.grid(pre = c(1), treat = c(0,1))
df1$t <- 1
df$arm <- 1:4
df1$arm <- c(2,4)
df <- rbind(df,df1)
```
so `df` now has six rows with the variables `pre`, `treat`, `t` and indicating the trial arm. Now we duplicate each row ten times, one for each cluster, and assign cluster labels:
```
df <- df[rep(1:6,each=10),]
df$cl <- c(1:40,11:20,31:40)
```
we can then create the different interaction variables described in the model:
```
df$intpost <- df$treat*df$post
df$prepost <- df$pre*df$post
df$intprepost <- df$pre*df$treat*df$post
```
To generate a randomly sampled number of individuals per cluster, we can draw uniformly from the numbers between two and ten and then duplicate each cluster's rows that many times:
```
n <- sample(1:10,J*4,replace=TRUE)
df$n <- n[df$cl]
df1 <- df[rep(1:nrow(df),df$n),]
```

### Coding the model
We can create the model in one command now:
```
mod <- Model$new(
  formula = ~ treat + post + intpost + prepost + intprepost + (1|gr(cl)) + (1|gr(cl,t)),
  mean = list(parameters = c(0,0,0,0,0,0.1)),
  covariance = list(parameters = c(0.25,0.10)),
  family = binomial(),
  var_par = 10,
  data = df1
)
```
See [model specification](../../glmmr/model/model_specification.md) for details on writing the formula. The `mean` argument gives the $\beta$ parameters, where we have started with a treatment effect of 0.1. The `covariance` argument gives the parameters $\tau_1$ and $\tau_2$, and we have also specified the family and data. The additional `var_par` here gives the extra scale parameter for the binomial distribution, here representing that there are ten `trials`.

### Estimating power
Estimating power is simple, we can just use the member function `mod$power()`, which gives:
```
mod$power()
>    Parameter Value        SE      Power
>1 [Intercept]   0.0 0.1272784 0.02500000
>2       treat   0.0 0.1796723 0.02500000
>3        post   0.0 0.1851434 0.02500000
>4     intpost   0.0 0.2580505 0.02500000
>5     prepost   0.0 0.1851434 0.02500000
>6  intprepost   0.1 0.2580851 0.05791791
```
So you can see the power for our parameter is a measly 5.8\%. It may be useful to estimate power over a range of values. Rather than re-generating the model each time, we can just update its parameters and re-calculate the power. To do this in a loop and extract the power ready for plotting we can do:
```
power.data <- data.frame(b = seq(0,2,by=0.02), power = NA)
for(i in 1:nrow(power.data)){
    mod$update_parameters(mean = c(0,0,0,0,0,power.data$b[i]))
    pwr <- mod$power()
    power.data$power[i] <- pwr$Power[6]
}
```
and then we can plot these data with `ggplot`:
```
require(ggplot2)
ggplot(data=power.data,aes(x=b,y=power))+
  geom_line()+
  theme_bw()
```
giving the figure below.

![Power](/power.jpeg)