# ecode
HOW TO INSTALL:
Please download file 'ecode_0.0.0.9000.tar.gz' and install it manually through 'Tools/Install Packages...' in RStudio.

After installation, run the following test code to make sure package ecode is ready for use.

<code>
##Example1: Lotka-Volterra competition model
library(ecode)
eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
x <- eode(dxdt = eq1, dydt = eq2)
x
plot(x)
</code>
