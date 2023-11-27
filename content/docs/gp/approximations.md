---
title: Gaussian Process Approximations
weight: 1
---

# Gaussian Processes
A Gaussian process is a useful non-parametric model for smooth functions. A finite set of realisations from a Gaussian process is multivariate Gaussian distributed; the Gaussian process is the distribution of the infinitely many possible such realisations. We can use a Gaussian process model in `glmmrBase`. Our linear mixed model specification, (i.e. a GLMM with Gaussian-identity family and link):
$$
y \sim N(X\beta + Z\mathbf{u},\sigma^2) 
$$
$$
\mathbf{u} \sim N(0,D)
$$
Let's assume that there is some variable $s$, and we want to include a smooth non-parametric function in our specification such that $\mu = X\beta + f(s)$. Assuming all the values of $s$ are unique then $Z$ above is just the identity matrix. The matrix $D$ is then defined by a covariance function $Cov(s,s') = g(| s-s'|;\theta)$ where $|.|$ is often the Euclidean distance. For example, a common covariance function is the exponential function:
$$
g(| s-s'|;\theta) = \tau^2 \exp \left(-\frac{|s-s'|}{\phi} \right)
$$
where $\tau^2$ is the magnitude of the function and $\phi$ the "lengthscale", which controls the smoothness. 

## Specifying the model in `glmmrBase`
We could specify this model in `glmmrBase` and simulate some data from it, assuming no fixed effects, as follows. We'll consider a two dimensional Gaussian process on $[-1,1] \times [-1,1]$ First, we can simulate 200 random locations:
```
df <- data.frame(x = runif(200, -1, 1), y = runif(200, -1, 1))
```
then the model is:
```
mod <- Model$new(
  formula = ~ 1 + (1|fexp(x,y)),
  data = df,
  covariance = list(parameters = c(0.5,1)),
  family = gaussian()
)
```
where we have set $\tau^2 = 0.5$ and $\phi=1.0$. The `sim_data` member function simulates both random effects and fixed effects:
```
y <- mod$sim_data(type = "all")
```
the argument `type="all"` ensures the function returns the simulated random effects along with all the other data. 

## Model fitting
We can fit the model using MCML or Laplace appoximation as standard in the package. However, these models quickly become very slow to fit. In particular, estimating the covariance parameters requires inverting the matrix $D$ (or decomposing it, see model estimation - link to be added) which scales as $O(n^3)$. For even reasonably sized data sets this can be prohibitively slow. As an alternative, we could use an approximation. There are two approximations included in `glmmrBase`: the nearest neighbour Gaussian process and the Hilbert Space Gaussian Process. We describe both of these in the next section and then provide some examples in the following section.

# Gaussian Process Approximations
## Nearest Neighbour Gaussian Process
The multivariate Gaussian likelihood can be rewritten as the product of the conditional densities:
$$
f(s) = f(s_1) \prod_{i=2}^n f(s_i|s_1,...,s_{i-1})
$$

Vecchia (1988) proposed that one can approximate $f(s)$ by limiting the conditioning sets for each $s_j$ to a maximum size of $m$. For geospatial applications, Datta (2016) proposed the nearest neighbour Gaussian process (NNGP), in which the conditioning sets are limited to the $m$ "nearest neighbours" of each observation. 

Let $\mathcal{N}_{ji}$ be the set of up to $m$ nearest neighbours of $j$ with index less than $j$; we return to the question of what "nearest" means in this context shortly. The approximation is:

$$
f(z) \approx f(z_1) \prod_{i=2}^n f(z_i|\mathcal{N}_i)
$$

which leads to:
$$
z_1 = \eta_1 
$$
$$
z_j = \sum_{i=1}^p a_{ji}z_{\mathcal{N}_{ji}}
$$

where $\mathcal{N}_{ji}$ is the $i$th nearest neighbour of $j$.

The equation above can be more compactly written as $z = Az + \eta$ where $A$ is a sparse, strictly lower triangular matrix, $\eta \sim N(0,D)$ with $D$ a diagonal matrix with entries $D_{11} = \text{Var}(z_1)$ and $D_{ii} = \text{Var}(w_i | \mathcal{N}_i)$ (Finley et al (2019)). The approximate covariance matrix can then be written as $\Sigma_0 \approx (I-A)^{-1}D(I-A)^{-T}$. We implement the algorithms described by Finley et al (2019) to generate the matrices $A$ and $D$. The cost of generating these matrices is $O(nm^3)$ owing to the need to solve $n-1$ linear systems of up to size $m$.

To approximate the multivariate Gaussian log likelihood we have 
$$
\log(|\Sigma_0|) \approx \sum_{i=1}^m -\log(D_{ii}) 
$$
and the quadratic form can be calculated using the algorithm in Finley (2019). The log determinant is $log(\vert \Sigma_0 \vert) \approx 2\sum_i D_{ii}$ The NNGP has computational complexity scaling with $O(Tnm^3)$ and so the log-likelhood is now linear in both time and grid size.

 There are two main factors that can affect the accuracy of the NNGP approximation. The ordering of the observations can affect performance of the approximation. Guinness et al (2018) compares the performance of several different ordering schemes. They suggest that a ``minimax'' scheme, in which the next observation in the order is the one which maximises the minimum distance to the previous observations, often performed best, although in several scenarios ordering in terms of the vertical coordinate, or at random provided comparable performance. The number of neighbours did not appear to affect the relative performance of the orderings; they considered both 30 and 60, although Datta (2016) uses 15 nearest neighbours in their simulations.

## Hilbert Space Gaussian Process
Low, or reduced, rank approximations aim to approximate the matrix $\Sigma$ with a matrix $\tilde{\Sigma}$ with rank $m < n$. The optimal low-rank approximation is $\tilde{\Sigma} = \Phi \Lambda \Phi^T$ where $\Lambda$ is a diagonal matrix of the $m$ leading eigenvalues of $\Sigma$ and $\Phi$ the matrix of the corresponding eigenvectors. However, generating the eigendecomposition scales the same as matrix inversion. Solin and Sarkka (2020) propose an efficient method to approximate the eigenvalues and eigenvectors using Hilbert space methods, so we refer to it as a Hilbert Space Gaussian Process (HSGP). Riutort-Mayol et al (2022) provide further discussion of these methods.

Stationary covariance functions, including those in the Matern class like exponential and squared exponential, can be represented in terms of their spectral densities. For example, the spectral density function of the squared exponential function in $D$ dimensions is:
$$
S(\omega) = \sigma^2 (\sqrt{2\pi})^D \phi^D \text{exp}(-\phi^2 \omega^2/2)
$$

Consider first a unidimensional space ($D=1$) with support on $[-L,L]$. The eigenvalues $\lambda_j$ (which are the diagonal elements of $\Lambda$) and eigenvectors $\phi_j$ (which form the columns of $\Phi$) of the Laplacian operator in this domain are:
$$
\lambda_j = \left( \frac{j\pi}{2L} \right)^2 
$$
and
$$
\phi_j(s) = \sqrt{\frac{1}{L}} \text{sin}\left(  \sqrt{\lambda_j}(x+L) \right)
$$
Then the approximation in one dimension is
$$
\mathcal{Z}(s) \approx \sum_{j=1}^m S\left(\sqrt{\lambda_j}\right)^{1/2}\phi_j(s)\beta_j
$$
where $\beta_j \sim N(0,1)$. This result can be generalised to multiple dimensions. The total number of eigenvalues and eigenfunctions in multiple dimensions is the combination of all univariate eigenvalues and eigenfunctions over all dimensions. For example, if there were two dimensions and $m=2$ then there would be four multivariate eigenfunctions and eigenvalues equal to the four combinations over the two dimensions. The approximation is then as described above but with the multivariate equivalents.

The approximation can be used in the full model such that the intensity becomes
$$
\mu = X\beta + \Phi\Lambda^{\frac{1}{2}} v
$$
and where $v \sim N(0,I)$. The matrix $Phi$ does not depend on the covariance parameters and can be pre-computed, so only the product $\Phi\Lambda^{\frac{1}{2}}$ needs to be re-calculated during model fitting, which scales as $O(nm^2)$. Similarly, the values of $\Lambda$ do not depend on the positions of the observations (i.e. the value of $s$) and so prediction at new locations just requires generating a new matrix $\Phi$.

# Examples
We will fit a a model using the HSGP. TO BE COMPLETED...