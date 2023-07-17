---
title: Creating data
weight: 1
---

# Creating data
The package includes the function `nelder()`, which we use to generate data for the examples below. [Nelder (1965)](https://royalsocietypublishing.org/doi/10.1098/rspa.1965.0012) suggested a simple notation that could express a large variety of different blocked designs. The notation was proposed in the context of split-plot experiments for agricultural research, where researchers often split areas of land into blocks, sub-blocks, and other smaller divisions, and apply different combinations of treatments. However, the notation is useful for expressing a large variety of experimental designs with correlation and clustering including cluster trials, cohort studies, and spatial and temporal 
prevalence surveys. We have included the function `nelder()` that generates a data frame of a design using the notation. 

There are two operations:
* `>` (or $\to$ in Nelder's notation) indicates "clustered in".
* `*` (or $\times$ in Nelder's notation) indicates a crossing that generates all combinations of two factors.

The function takes a formula input indicating the name of the variable and a number for the number of levels, such as `abc(12)`. 
So for example `~cl(4) > ind(5)` means in each of five levels of `cl` there are five levels of `ind`, and the individuals are different between clusters. The formula `~cl(4) * t(3)` indicates that each of the four levels of `cl` are observed for each of the three levels of `t`. Brackets are used to indicate the order of evaluation.