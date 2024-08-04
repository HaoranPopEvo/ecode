# ecode

<!-- badges: start -->
[![R-CMD-check](https://github.com/Asa12138/ecode/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Asa12138/ecode/actions/workflows/R-CMD-check.yaml)
[![](https://img.shields.io/badge/doi-10.1016/j.ecolmodel.2024.110676-yellow.svg)](https://doi.org/10.1016/j.ecolmodel.2024.110676)
[![](http://cranlogs.r-pkg.org/badges/grand-total/ecode)](https://cran.r-project.org/package=ecode)
[![](https://www.r-pkg.org/badges/version/ecode?color=green)](https://cran.r-project.org/package=ecode)
[![](https://img.shields.io/badge/devel%20version-0.1.0-green.svg)](https://github.com/Asa12138/ecode)
<!-- badges: end -->

`ecode`, a novel package for modelling ecological populations and communities using ordinary differential equation systems, designed with a user-friendly framework. 

By following a three-cycle procedure, users can easily construct ecological models and explore their behaviors through a wide range of graphical, analytical, and numerical techniques. 

The package incorporates advanced techniques such as grid search methods and simulated annealing algorithms, enabling users to iteratively refine their models and achieve accurate predictions. 

Notably, ecode minimises external dependencies, ensuring robustness and reducing the risk of package failure caused by updates in dependencies. 

Please find the user manual at https://bookdown.org/Asa12138/ecode_book/.

## Citation

To cite `ecode` in publications use:

Haoran Wu,
ecode: An R package to investigate community dynamics in ordinary differential equation systems,
Ecological Modelling,
2024,
<https://doi.org/10.1016/j.ecolmodel.2024.110676>.


## Installation

You can install the released version of `ecode` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ecode")
```

You can install the development version of `ecode` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("HaoranPopEvo/ecode")
```


Test code:

```{r}
##Example1: Lotka-Volterra competition model
library(ecode)
eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
x <- eode(dxdt = eq1, dydt = eq2)
x
plot(x)
```

