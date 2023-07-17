---
title: The Girling Algorithm
weight: 1
---

# Introduction
The Girling Algorithm is a method to find the optimal design weights for an experiment with correlated observations. (To be completed).

# Statistical Details
A "design space" consists of all the possible observations one could make in an experiment. Each unique design "point" is written as $x_i$ for $i=1,...,N$, which includes both the covariate values and "design matrix" for the random effects. The optimal design problem is to find the optimal $n<N$ observations, where "optimum" here is defined by the c-optimality criterion described below. One way of tackling this problem is to identify an approximate optimal design. An approximate design $d$ can be characterised by the tuple $d:= \\{ (x_i,w_i):i=1,...,N, \sum_i w_i = 1 \\}$. We then base our analysis on the linear mixed model:
$$
y = X\beta + Z\mathbf{u} + e
$$
where $y$ is the $n$-length vector of outcomes, $X$ is the $n \times P$ covariate matrix with parameters $\beta$, and $Z$ is the $n \times Q$ design matrix for the random effects. The random effects are specified as $\mathbf{u}\sim N(0,D)$ where $D$ is the covariance matrix of the random effects and we have the independent error term $e \sim N(0,\sigma^2I)$. The covariance matrix of $y$ can be written as $\Sigma = \frac{\sigma^2}{n}W + ZDZ^T$; $W$ is a diagonal matrix with elements equal to $W_{ii} = \frac{1}{w_i}$ where $w_i$ are the weights from before.

Summarising here the analysis from Girling (2023), we use $c$ to represent a vector such that our objective is to estimate $c^T\beta$. The c-optimality criterion is $c^T M^{-1} c$, where $M = X^T \Sigma^{-1} X$ is the information matrix for $\beta$. For most examples of interest, $c$ is a vector comprised of zeroes except for a one in the place of a treatment effect parameter. We write $c^T \beta = \theta$. The best linear unbiased estimator for $\theta$ is:
$$
\hat{\theta} = a^Ty
$$
We can now use the classic generalised least squares approach by pre-multiplying the model by $L = \Sigma^{-1/2}$ where $\Sigma^{-1/2}$ is a square root of $\Sigma^{-1}$; we use the Cholesky decomposition here. We then set $F=LX$ and rewrite our estimator as:
$$
\hat{\theta} = a_1^TLy
$$
where $a = L^Ta_1$. To ensure unbiasedness of the estimator, and by the Gauss-Markov theorem, we have $F^Ta_1 = c$ for $a_1 = F(F^TF)^{-1}c$ and hence $a = \Sigma^{-1}X(X^T\Sigma^{-1}X)^{-1}c$. This also gives us the GLS estimator.

For many of the types of design problems we are interested in, the observations are grouped into independent clusters. For example, observations in clusters in a cluster trial, or observations within individuals in a cohort design. In these cases the covariance matrices are block diagonal, for example:
$$
D = 
\begin{bmatrix}
D_1 & 0 & \dots & 0 \\\
0 & D_2 & \dots & 0 \\\
\vdots & \vdots & \ddots & 0 \\\
0 & 0 & \dots & D_K 
\end{bmatrix}
$$
where there are $K$ "clusters". In these cases we can simplify somewhat some of the calculations, for example, the information matrix can be written as $M = \sum_{k=1}^K X_k^T \Sigma_k^{-1} X_k$ where $X_k$ and $\Sigma_k$ refer to the sub-matrices relating to the observations in cluster $k$. Using this, we aim to work out the variance of $\hat{\theta}$:
$$
\text{Var}(\hat{\theta}) = a_1^T a_1 = a^T\Sigma a
$$
$$
= \sum_{k=1}^K a_k^T \Sigma_k a_k
$$ 
$$
= \sum_{k=1}^K \sum_{i=1}^{n_k} \sum_{j=1}^{n_k} a_{ki} a_{kj} \Gamma_{ij} + \frac{\sigma^2}{n} \sum_{k=1}^K \sum_{i=1}^{n_k} \frac{a_{ki}^2}{w_{ki}}
$$
where $n_k$ is the number of observation in cluster $k$, $\Gamma_{ij}$ is the covariance between observation $i$ and $j$ defined by some covariance function, and $a_k = [a_{k1},a_{k2},...,a_{kn_k}]^T$.

Now, all this is explained in order to work towards a strategy for identifying the weights of the optimal approximate design. By the Cauchy-Schwarz inequality (and ignoring the first part of the variance expression), we know that:
$$
\frac{\sigma^2}{n} \sum_{k=1}^K \sum_{i=1}^{n_k} \frac{a_{ki}^2}{w_{ki}} \geq \frac{\sigma^2}{n} \left ( \sum_{k=1}^K \sum_{i=1}^{n_k} \vert a_{ki}\vert \right)^2
$$
which gives us a lower bound on the variance of the estimator. For the design to be optimal, we require the variance to be minimised, by our c-optimality criterion, which is achieved if we set the weights as:
$$
w_{ki} = \frac{\vert a_{ki} \vert}{\sum_k \sum_i \vert a_{ki} \vert}
$$

## Optimal design algorithm
A general form of the algorithm can be stated to find the optimal weights $w$. 
  1. Set $a = \Sigma^{-1}X(X^T \Sigma^{-1} X)^{-1} c = \Sigma^{-1}XM^{-1} c$
  2. Update $w_{ki} = \frac{\vert a_{ki} \vert}{\sum_k \sum_i \vert a_{ki} \vert}$
  3. Set $\Sigma = \frac{\sigma^2}{n}\text{diag}(w^{-1})+ ZDZ^T$
And iterate until covergence.

# Implementing the algorithm
We implement the algorithm in C++ for speed efficiencies using the `glmmrBase` package, which handles all the model specification and processing of the model, including identifying the independent clusters (see links in the menu for more details). The Eigen linear algebra library is used. There are several steps we can take to try to improve computational time, such as by pre-calculating certain matrix products.

## Pre-calculation
At the start, we calculate and store several model components:
  * A vector containing the matrices $A_k = Z_k D_k Z_k^T$
  * A vector containing the matrices $X_k$

We also reserve a vector containing matrices of the size of $\Sigma_k$, labelled $S_k$ here, for calculation, store a vector of iterators for the observations in each cluster, and create a $P \times P$ matrix $M$.

## Checking for zero-weights
Many optimal designs will not contain all the design points so the weight will become infinitesimal. While this could cause some floating point errors, it is also inefficient to keep trying to calculate their weights. So on each iteration of the algorithm we first check if any of the weights are below some threshold -- we use $10^{-8}$. Then the check is:
  * If any $w < 10^{-8}$ then,
    * For each $k$ in 1 to $K$
      * For each $i$ in cluster $k$
        * If $w_{ki} < 10^{-8}$ then
          1. Remove row $i$ from $X_k$
          2. Remove row and column $i$ from $A_k$
          3. Resize $S_k$
          4. Remove $i$ from $k$'s iterators

## Algorithm
The full algorithm as implemented is then as follows. We start by setting $\delta = 1$ and we also have our vector $c$ and a vector $w_{new}$ for the new set of weights.
* While $\delta > tol$
  * Check for zero-weights (see above)
  * Set $M = \mathbf{0}$
  * For each $k$ in 1 to $K$
    * Set $S_k = A_k$
    * For each $i$ in cluster $k$
      * Set $S_{k,ii} += \frac{\sigma^2}{nw_{ki}}$
    * Set $S_k = S_k^{-1}$ (we use LDL inverse as it is fastest)
    * Set $M += X_k^T S_k X_k$
  * If $M$ not positive semidefinite then
    * Generate $m$ as the row sums of $M$
    * For each $p$ in 1 to $P$ 
      * If $m_p = 0$ then
        * Remove column $p$ from all $X_k$
        * Remove element $p$ from $c$
    * Recalculate $M$ as $M = \sum_k X_k^T S_k X_k$
  * Generate vector $v = M^{-1}c$
  * Set $w_{new}=0$
  * For each $k$ in 1 to $K$
    * Set $w_{new,k} = S_k X_k v$
  * Set $w_{new} = \vert w_{new} \vert$
  * Update $w_{new} = w_{new}/\sum w_new$
  * Update $\delta = \max (\vert w - w_{new} \vert )$
where we've used the programming syntax $a += b$ to mean $a = a + b$.

# Example with code
We provide an example here on executing the above algorithm using `glmmrBase`. The algorithm isn't exposed to the user in the released package yet, as it hasn't been published and fully tested, but it is still hiding in the package! We will aim to find the optimal weights for cluster-periods in a stepped-wedge cluster trial design with six cluster sequences and seven time periods. The links in the menu to the left provide more information on model specification and other features (which are due to be completed), so the exposition here will be fairly bare bones.

We first generate a dataset with six cluster and seven time periods with ten observation per cluster, adding in the treatment effect
```
df <- nelder(~ (cl(6)*t(7)) > ind(10))
df$int <- 0
df[df$t > df$cl, 'int'] <- 1
```
As the algorithm calculates weights for unique experimental conditions, we could set the number of individuals to one per cluster-period and then indpendently tell the algorithm the total sample size $N$. The alternative, and the approach taken below,
is to generate a design with all the observations as above, and then the function will just use the total number of observations as 
$N$.

Now, we will generate a linear mixed model with autoregressive covariance function (although one can use any covariance function with this package). In particular, the covariance function for cluster $k$ and time $t$, if $z_{kt}\mathbf{u} = \alpha_{kt}$ is the random effect term, is:
$$
\text{Cov}(\alpha_{kt},\alpha_{k',t'}) = \begin{cases}
\tau^2 \lambda^{\vert t - t' \vert} & \text{ if } k=k' \\\
0 & \text{ otherwise}
\end{cases}
$$
We will choose covariance parameter values of $\sigma^2 = 1$, $\tau^2 = 0.05$ and $\lambda = 0.8$. The model is set up as:
```
model <- Model$new(
  formula = ~ factor(t) + int - 1 + (1|gr(cl)*ar1(t)),
  covariance = list(parameters = c(0.05,0.8)),
  family = gaussian(),
  data = df
)
```
The default values are $\sigma^2 = 1$ and $\beta$ is initialised to random values, so these are not directly set here. The Girling Algorithm has now been added to `glmmrOptim`. As with all other algorithms in the package, we first create a design space:
```
ds <- DesignSpace$new(des)
```
And, then specify the option `algo="girling"` in the call to `optim()`:
```
outb <- ds$optimal(m=1,C = list(c(rep(0,5),1)),verbose = TRUE,algo="girling")
```
The first argument `m` is redundant so we just set it to one. This takes about 50 milliseconds to run. And produces the weights shown in the figure.

![Optimal design weights](images/example_girling_plot.jpg)

## Rounding weights
There are several ways of rounding design weights to whole numbers of clusters, which are described in Watson & Pan (2023). The function `apportion` in our `glmmrOptim` package calculates all of the different rounding schemes. It turns out to be very slow for 42 weights because the method proposed by Puckelsheim and Rieder calculates all possible combinations of design points with zero weight and so increases exponentially. For this example, we just use Hamilton's method to round the weights, which seems to perform the best in almost all the tests we have done. This rounding is quite simple; we will create a design with 100 observations:
```
n <- 100
des <- floor(w*100)
if(sum(des)<n){
  while(sum(des)<n){
    rem <- n*w - des
    des[which.max(rem)] <- des[which.max(rem)]+1
  }
}
```
which produces the totals in the figure below.

![Optimal design weights](images/example_girling_rounded.jpg)