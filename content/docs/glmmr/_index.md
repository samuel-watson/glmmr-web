---
weight: 1
bookFlatSection: true
title: "The `glmmr` packages"
---

# `glmmrBase`
The currently available packages are `glmmrBase`, which provides the foundations including model specification, analysis, and estimation, and `glmmrOptim`, which implements a set of optimal design algorithms, with more packages planned in future.

## Installation
Both packages can be installed from CRAN using, for example
```
install.packages("glmmrBase")
```
or the most recent release is available from GitHub as:
```
devtools::install_github("samuel-watson/glmmrBase")
```

A pre-compiled binary is also available with each release on the (GitHub page)[https://github.com/samuel-watson/glmmrBase]. 

## Building from source

It is strongly recommended to build from source with the flags `-fno-math-errno -O3 -g`, this will cut the time to run many functions by as much as 90%. One way to do this is to set CPP_FLAGS in `~/.R/Makevars`. Another alternative is to download the package source `.tar.gz` file and run from the command line 
```
R CMD INSTALL --configure-args="CPPFLAGS=-fno-math-errno -O3 -g" glmmrBase_0.5.2.tar.gz
```
Yet another alternative would be to clone the package from this repository and edit the `makevars` or `makevars.win` file to include the flags. The package can then be installed using `devtools::install`.
