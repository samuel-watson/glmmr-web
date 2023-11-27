---
title: Specifying a formula
weight: 1
---

# Specifying the formula
The `Model` class uses a similar formula specification to other packages in R. It aims to be fairly flexible and the user can specify separate formulae for mean and covariance functions, or a single formula for both. Specifying the fixed effect component is the same as for other functions in R, for instance `~ x1 + x2`. Specifying the covariance function for the random effects extends the approaches of packages like `lme4` and `glmmTMB`.

## Specifying the covariance formula
A covariance function is specified as an additive formula made up of components with structure `(1|f(j))`. The left side of the vertical bar specifies the parameters of the which covariates in the model have a random effects structure. For example, if we wanted to include a random slope for a covariate `x`, we could write `(x|f(j))`. 

The right side of the vertical bar specify the covariance function `f` for that term using variable named in the data `j`. Covariance functions on the right side of the vertical bar are multiplied together, i.e., `(1|f(j)*g(t))`. The table below shows the currently implemented covariance functions

| Function | $Cov(x_i,x_{i'})$ | $\theta$ | `Code` |
|----------|-------------------|----------|--------|
| Identity/ Group membership | $\theta_1^2 \mathbf{1}(x_i = x_{i'})$ | $\theta_1 > 0$ | `gr(x)` |
| Exponential | $\theta_1 \text{exp}(- \vert x_i - x_{i'}\vert / \theta_2 )$ | $\theta_1,\theta_2 > 0$ | `fexp(x)`|
| | $\text{exp}(- \vert x_i - x_{i'}\vert /\theta_1)$ | $\theta_1 > 0$ | `fexp0(x)` |
| Squared Exponential | $\theta_1 \text{exp}(- (\vert x_i - x_{i'}\vert / \theta_2)^2)$ | $\theta_1,\theta_2 > 0$ | `sqexp(x)` |
| | $\text{exp}(-( \vert x_i - x_{i'}\vert/\theta_1)^2 )$ | $\theta_1 > 0$ | `sqexp0(x)` |
| Autoregressive order 1 | $\theta_1^{\vert x_i - x_{i'} \vert}$ | $0 < \theta_1 < 1$ | `ar1(x)` |
| Bessel | $K_{\theta_1}(x)$ | $\theta_1$ > 0 | `bessel(x)` |
| Matern | $\frac{2^{1-\theta_1}}{\Gamma(\theta_1)}\left( \sqrt{2\theta_1}\frac{x}{\theta_2} \right)^{\theta_1} K_{\theta_1}\left(\sqrt{2\theta_1}\frac{x}{\theta_2})\right)$ | $\theta_1,\theta_2 > 0$ | `matern(x)` |
| Compactly supported* | || |
| Wendland 0 | $\theta_1(1-y)^{\theta_2}, 0 \leq y \leq 1; 0, y \geq 1$ | $\theta_1>0, \theta_2 \geq (d+1)/2$ | `wend0(x)` |
| Wendland 1 | $\theta_1(1+(\theta_2+1)y)(1-y)^{\theta_2+1}, 0 \leq y \leq 1; 0, y \geq 1$ | $\theta_1>0, \theta_2 \geq (d+3)/2$ | `wend1(x)` |
| Wendland 2 | $\theta_1(1+(\theta_2+2)y + \frac{1}{3}((\theta_2+2)^2 - 1)y^2)(1-y)^{\theta_2+2}, 0 \leq y \leq 1 $ | $\theta_1>0,\theta_2 \geq (d+5)/2$ | `wend1(x)` |
| | $0, y \geq 1$ | | |
| Whittle-Matern $\times$ Wendland** | $\theta_1\frac{2^{1-\theta_2}}{\Gamma(\theta_2)}y^{\theta_2}K_{\theta_2}(y)(1+\frac{11}{2}y + \frac{117}{12}y^2)(1-y), 0 \leq y \leq 1; 0, y \geq 1$ | $\theta_1,\theta_2 > 0$ | `prodwm(x)` |
| Cauchy $\times$ Bohman*** | $\theta_1(1-y^{\theta_2})^{-3}\left( (1-y)\text{cos}(\pi y)+\frac{1}{\pi}\text{sin}(\pi y) \right), 0 \leq y \leq 1; 0, y \geq 1$ | $\theta_1>0, 0 \leq \theta_2 \leq 2$ | `prodcb(x)` |
| Exponential $\times$ Kantar**** | $\theta_1\exp{(-y^{\theta_2})}\left( (1-y)\frac{\sin{(2\pi y)}}{2\pi y} + \frac{1}{\pi}\frac{1-\cos{(2\pi y)}}{2\pi y} \right), 0 \leq y \leq 1$ | $\theta_1,\theta_2 > 0$ | `prodek(x)` |
| | $0, y \geq 1$ | | |

$\vert . \vert$ is the Euclidean distance. $K_a$ is the modified Bessel function of the second kind. 
*Variable $y$ is defined as $x/r$ where $r$ is the effective range. For the compactly supported functions $d$ is the number of dimensions in `x`. 
**Permissible in one or two dimensions. ***Only permissible in one dimension. ****Permissible in up to three dimensions.

One combines functions to provide the desired covariance function. For example, for a stepped-wedge cluster trial we could consider the standard specification with an exchangeable random effect for the cluster level, and a separate exchangeable random effect for the cluster-period, which would be `~(1|gr(j))+(1|gr(j,t))`. Alternatively, we could consider an autoregressive cluster-level random effect that decays exponentially over time so that, for persons $i$ in cluster $j$ at time $t$, $Cov(y_{ijt},y_{i'jt}) = \theta_1$, for $i\neq i'$, $Cov(y_{ijt},y_{i'jt'}) = \theta_1 \theta_2^{\vert t-t' \vert}$ for $t \neq t'$, and $Cov(y_{ijt},y_{i'j't}) = 0$ for $j \neq j'$. This function would be specified as `~(1|gr(j)*ar(t))`.

The above elements can be combined to create relatively complex random effects structures such as `(x|gr(j))+(1|fexp(x,y)*ar(t))`. The current version does not yet support correlated random effects terms, but this is planned for a future version. Additional functions can be added on request.

## Adding the fixed effects
The fixed effects are specified in the same was as other R packages and the random effects can just be added onto the end. For example, to continue our stepped-wedge cluster trial example, we may want to include fixed effect indicators for each time period `t`, and indicator for the intervention `int`, and remove the intercept: `~ factor(t)+ int - 1 + (1|gr(j)*ar1(t))`.

We can also specify non-linear fixed effects models and include desired parameter names in the formula. For example, the fixed effect specification
$$
\eta_i = \beta_0 + \beta_1*int_i + \beta_2\exp(\beta_3 x_i)
$$
can be written as `~ int + beta_2*exp(beta_3 * x)`. The interpreter assumes that any token (i.e. string of characters) that isn't the name of a variable or function is the name of a parameter. The formula is stored as a sequence of mathematical operations. By default, any data names will have a parameter attached to them, for example `exp(-x/2)` will be interpreted as $\exp(- \beta_x x / 2)$. To include data without a parameter, wrap it in brackets, for example `exp(-(x)/2)` will be interpreted as $\exp(- x / 2)$. This approach could also be used to add offsets, such as `+ (x) + ...`. One can check how the program has interpreted the formula using `model$print_instructions()`, which will print out the steps used to calculate the mean. An auto-differentiation scheme is used to calculate the relevant derivative for model fitting and other operations.

As an example using the model and printing the instructions gives us:
```
model <- Model$new(
    formula = ~ int + beta_2*exp(beta_3 * x) + (1|gr(j)*ar1(t)),
    data = data,
    family = gaussian()
)
model$print_instructions()

```

