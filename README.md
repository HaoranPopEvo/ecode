# ecode

<!-- badges: start -->
[![](https://img.shields.io/badge/doi-10.1016/j.ecolmodel.2024.110676-yellow.svg)](https://doi.org/10.1016/j.ecolmodel.2024.110676)
<!-- badges: end -->

`ecode`, a novel package for modelling ecological populations and communities using ordinary differential equation systems, designed with a user-friendly framework. 

By following a three-cycle procedure, users can easily construct ecological models and explore their behaviours through a wide range of graphical, analytical, and numerical techniques. 

The package incorporates advanced techniques such as grid search methods and simulated annealing algorithms, enabling users to iteratively refine their models and achieve accurate predictions. 

Notably, ecode minimises external dependencies, ensuring robustness and reducing the risk of package failure caused by updates in dependencies. 

## Citation

To cite ReporterScore in publications use:

Haoran Wu,
ecode: An R package to investigate community dynamics in ordinary differential equation systems,
Ecological Modelling,
2024,
<https://doi.org/10.1016/j.ecolmodel.2024.110676>.


## Installation

The `ecode` package is currently only available at GitHub. Installation can be done by the following code:

```{r}
library(devtools)
install_github("HaoranPopEvo/ecode")
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

