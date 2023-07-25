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

## Combinatorial algorithm
The same function `optimal` is used to run the other algorithms. The argument `algo` specifies which algorithm to run, the options are:
 - `1` Local search
 - `2` Greedy search 
 - `3` Reverse greedy search
 - `"girling"` Girling's algorithm
For the combinatorial algorithms, we can also string together numbers to run one algorithm after another, such as `c(3,1)`. 

Now, to find a design of size 80 using the reverse greedy search:
```
out <- ds$optimal(80,C = list(c(rep(0,6),1)),verbose = TRUE,algo=3)
```
where we have included the vector $C$ as an element of a list. We've also specified `verbose=TRUE` so we can have detailed output of the progress of the algorithm. The reverse greedy search is relatively slow at first, especially for larger designs spaces, and then speeds up as the remaining number of observations declines. To then reproduce the plots in the article with the resulting design we can use the following `ggplot2` command, where we've also made use of the package `colorRamps`.
```
require(ggplot2)
df2 <- df[out[[i]]$rows,]
dfa2 <- aggregate(df2$int,list(cl=df2$cl,t=df2$t),length)
dfa2$int <- aggregate(df2$int,list(cl=df2$cl,t=df2$t),mean)$x
dfa2 <- dfa2[order(dfa2$cl,dfa2$t),]
dfa2$method <- paste0("Combin\n",round(out[[i]]$val,6))
dfa2$total2 <- factor(dfa2$x,levels=1:10,ordered=TRUE)
ggplot(data=dfa2,aes(x=factor(t),y=factor(cl),fill=total2))+
             geom_tile(color="black")+
             geom_text(aes(label=int))+
             scale_fill_manual(name="Count",values=colorRamps::blue2green2red(10),breaks=1:10,drop=FALSE)+
             facet_wrap(~method)+
             theme_bw()+
             theme(panel.grid=element_blank(),axis.text.y = element_blank())+
             labs(x="Time",y="Cluster")+
             scale_x_discrete(expand = c(0,0))+
             scale_y_discrete(expand = c(0,0))

```

## Girling algorithm
The "Girling algorithm" identifies the optimal weights to place in each cluster-period. It is run in the same way as the combinatorial algorithms above, but the function returns a data frame giving all the unique cluster-periods and their respective weights.
```
outg <- ds$optimal(80,C = list(c(rep(0,6),1)),verbose = TRUE,algo="girling")
```
Then we can plot the result directly. In the article we used a colour palette from the package `scico`:
```
ggplot(data=outg,aes(x=factor(t),y=factor(cl),fill=weight))+
             geom_tile(color="black")+
             geom_text(aes(label=int))+
             scale_fill_scico(name = "Weight",palette="lajolla",limits=c(0,0.15))+
             facet_wrap(~method)+
             theme_bw()+
             theme(panel.grid=element_blank(),axis.text.y = element_blank())+
             labs(x="Time",y="Cluster")+
             scale_x_discrete(expand = c(0,0))+
             scale_y_discrete(expand = c(0,0))
```

The complete code to reproduce the figures with all the sub-figures is provided at the end of this page.

# Non-linear models
Non-linear models are run identically to their linear counterparts described above. For example, we could create a binomial-logistic design using the same covariance parameter values as above, and with zero parameters for the fixed effects:
```
model_bin <- Model$new(
  formula = ~ factor(t) + int - 1 + (1|gr(cl)) + (1|gr(cl,t)),
  covariance = list(parameters = c(0.05*0.8/0.95,0.05*0.2/0.95)),
  mean = list(parameters = rep(0,7))
  family = binomial(),
  data = df
)
```
We can then use this model object when creating a new design space. See the links in the menu for more information on the available models in the package.

# Robust optimal designs
Finally, we can extend a design space to include multiple designs and then aim to identify the robust optimal design. Let's create a second design like `model` above, but with different covariance parameters (an ICC of 0.005).
```
model2 <- Model$new(
  formula = ~ factor(t) + int - 1 + (1|gr(cl)) + (1|gr(cl,t)),
  covariance = list(parameters = c(0.005*0.8/0.95,0.005*0.2/0.95)),
  family = gaussian(),
  data = df
)
```
Then, we can create a design space with two designs in it:
```
ds <- DesignSpace$new(model,model2)
```
And then a call to `optimal` will automatically average over the variances of the two designs. We can change the weights on each design by specifying the `weights` when creating the design space.
```
out <- ds$optimal(80,C = list(c(rep(0,6),1),c(rep(0,6),1)),verbose = TRUE,algo=3)
```
Notice that we must now provide two $C$ vectors, one for each design, as each design may include different parameters. The running time of the combinatorial algorithms increases linearly with the number of designs.

# Complete code for the figure
Here, we provide the complete code to reproduce the figures in the article, with individual observations as experimental units.
```
### Optimal design paper models ###
require(glmmrOptim)
require(ggplot2)
require(scico)
require(patchwork)

df <- nelder(formula(~ (cl(7) * t(6)) > ind(10)))
df$int <- 0
df[df$t >= df$cl,'int'] <- 1
des <- Model$new(covariance = list(formula = ~(1|gr(cl)*ar1(t)), 
                                   parameters = c(0.05,0.8)),
                 mean = list(formula = ~ factor(t)  + int - 1,
                             parameters = rep(0,7)),
                 var_par = 0.95,
                 data = df,
                 family=gaussian())
#
# modelling

dfpar <- expand.grid(m=20,icc=c(0.05,0.01,0.1),cac=c(0.8,0.2))

out <- list()
outg <- list()

for(i in 1:nrow(dfpar)){
  sigma <- 1-dfpar$icc[i]
  tau <- dfpar$icc[i]#*dfpar$cac[i]
  #omega <- dfpar$icc[i]*(1-dfpar$cac[i])
  lambda <- dfpar$cac[i]
  des$update_parameters(cov.pars = c(tau,lambda),var.par = sigma)
  
  ds <- DesignSpace$new(des,experimental_condition = 1:nrow(df))
  out[[i]] <- ds$optimal(80,C = list(c(rep(0,6),1)),verbose = TRUE,algo=1)
  outg[[i]] <- ds$optimal(80,C = list(c(rep(0,6),1)),verbose = TRUE,algo="girling")
  
  df2 <- df[out[[i]]$rows,]
  dfa2 <- aggregate(df2$int,list(cl=df2$cl,t=df2$t),length)
  dfa2$int <- aggregate(df2$int,list(cl=df2$cl,t=df2$t),mean)$x
  dfa2 <- dfa2[order(dfa2$cl,dfa2$t),]
  
  dfa2$method <- paste0("Combin\n",round(out[[i]]$val,6))
  dfa2$total2 <- factor(dfa2$x,levels=1:10,ordered=TRUE)
  
  if(i%in%1){
    assign(paste0("p",dfpar$m[i],"_",dfpar$icc[i],"_",dfpar$cac[i]),
           ggplot(data=dfa2,aes(x=factor(t),y=factor(cl),fill=total2))+
             geom_tile(color="black")+
             geom_text(aes(label=int))+
             scale_fill_manual(name="Count",values=colorRamps::blue2green2red(10),breaks=1:10,drop=FALSE)+
             facet_wrap(~method)+
             theme_bw()+
             theme(panel.grid=element_blank(),axis.text.y = element_blank())+
             labs(x="Time",y="Cluster")+
             scale_x_discrete(expand = c(0,0))+
             scale_y_discrete(expand = c(0,0)),envir = .GlobalEnv)
  } else {
    assign(paste0("p",dfpar$m[i],"_",dfpar$icc[i],"_",dfpar$cac[i]),
           ggplot(data=dfa2,aes(x=factor(t),y=factor(cl),fill=total2))+
             geom_tile(color="black")+
             geom_text(aes(label=int))+
             scale_fill_manual(name="Count",values=colorRamps::blue2green2red(10),breaks=1:10,drop=FALSE)+
             facet_wrap(~method)+
             theme_bw()+
             theme(panel.grid=element_blank(),axis.text.y = element_blank(),legend.position = "none")+
             labs(x="Time",y="Cluster")+
             scale_x_discrete(expand = c(0,0))+
             scale_y_discrete(expand = c(0,0)),envir = .GlobalEnv)
  }
  
  outg[[i]]$method <- paste0("Weights\n")
  outg[[i]] <- outg[[i]][outg[[i]]$weight > 1e-5,]
  
  if(i%in%1){
    assign(paste0("q",dfpar$m[i],"_",dfpar$icc[i],"_",dfpar$cac[i]),
           ggplot(data=outg[[i]],aes(x=factor(t),y=factor(cl),fill=weight))+
             geom_tile(color="black")+
             geom_text(aes(label=int))+
             scale_fill_scico(name = "Weight",palette="lajolla",limits=c(0,0.15))+
             facet_wrap(~method)+
             theme_bw()+
             theme(panel.grid=element_blank(),axis.text.y = element_blank())+
             labs(x="Time",y="Cluster")+
             scale_x_discrete(expand = c(0,0))+
             scale_y_discrete(expand = c(0,0)),envir = .GlobalEnv)
  } else {
    assign(paste0("q",dfpar$m[i],"_",dfpar$icc[i],"_",dfpar$cac[i]),
           ggplot(data=outg[[i]],aes(x=factor(t),y=factor(cl),fill=weight))+
             geom_tile(color="black")+
             geom_text(aes(label=int))+
             scale_fill_scico(name="Weight",palette="lajolla",limits=c(0,0.15))+
             facet_wrap(~method)+
             theme_bw()+
             theme(panel.grid=element_blank(),axis.text.y = element_blank(),legend.position = "none")+
             labs(x="Time",y="Cluster")+
             scale_x_discrete(expand = c(0,0))+
             scale_y_discrete(expand = c(0,0)),envir = .GlobalEnv)
  }
  
}

row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="ICC = 0.05", angle = 90) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="ICC = 0.01", angle = 90) + theme_void() 
row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="ICC = 0.1", angle = 90) + theme_void() 
#col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label=	"\u03BB = 0.8") + theme_void() 
#col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="lambda = 0.5") + theme_void()
#col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="\u03BB = 0.2") + theme_void()
col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label=	"CAC = 0.8") + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="CAC = 0.2") + theme_void()

scol1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Combinatorial") + theme_void() 
#col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="lambda = 0.5") + theme_void()
scol2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Weights") + theme_void()

layoutplot <- "
#uuuuuuvvvvvv
#cccdddiiiqqq
aeeefffjjjrrr
aeeefffjjjrrr
bggghhhlllsss
bggghhhlllsss
pmmmnnnooottt
pmmmnnnooottt
"

plotlist <- list(a = row1, b = row2, p=row3, c=col1,d=col3,i=col1, 
                 q = col3, u = scol1, v = scol2,
                 e=p20_0.1_0.8,f=p20_0.1_0.2,j=q20_0.1_0.8,r=q20_0.1_0.2,
                 g=p20_0.05_0.8,h=p20_0.05_0.2,l=q20_0.05_0.8,s=q20_0.05_0.2,
                 m=p20_0.01_0.8,n=p20_0.01_0.2,o=q20_0.01_0.8,t=q20_0.01_0.2)

wrap_plots(plotlist, guides = 'collect', design = layoutplot) # 15*9

```