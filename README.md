# ecode

**HOW TO INSTALL**

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

