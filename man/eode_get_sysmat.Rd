% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability.R
\name{eode_get_sysmat}
\alias{eode_get_sysmat}
\title{System Matrix}
\usage{
eode_get_sysmat(x, value, delta = 0.00001)
}
\arguments{
\item{x}{object of class "\code{eode}" representing an ODE system.}

\item{value}{an object of class "\code{pp}" representing a phase point in the
ODE system under consideration.}

\item{delta}{Spacing or step size of a numerical grid used for calculating
numerical differentiation.}
}
\value{
An object of class "\code{matrix}". Each matrix entry, (i,j), represents
the partial derivative of the i-th equation of the system with respect to the
j-th variable taking other variables as a constant.
}
\description{
Linearlise an ODE system and calculate its system matrix.
}
\examples{
## Example1: Lotka-Volterra competition model
eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
x <- eode(dxdt = eq1, dydt = eq2)
eode_get_sysmat(x, value = pp(list(x = 1, y = 1)), delta = 10e-6)

## Example2: Susceptible-infected model
dX_Cdt <- function(X_C, Y_C, X_A, Y_A, nu = 0.15, beta = 0.1, mu = 0.15, g = 0.04) {
  nu * (X_A + Y_A) - beta * X_C * (Y_C + Y_A) - (mu + g) * X_C
}

dY_Cdt <- function(X_C, Y_C, Y_A, beta = 0.1, mu = 0.15, g = 0.04, rho = 0.2) {
  beta * X_C * (Y_C + Y_A) - (mu + g + rho) * Y_C
}

dX_Adt <- function(X_C, Y_C, X_A, Y_A, beta = 0.1, g = 0.04) {
  g * X_C - beta * X_A * (Y_C + Y_A)
}

dY_Adt <- function(X_A, Y_C, Y_A, beta = 0.1, g = 0.04, rho = 0.2) {
  beta * X_A * (Y_C + Y_A) + g * Y_C - rho * Y_A
}

x <- eode(
  dX_Cdt = dX_Cdt, dY_Cdt = dY_Cdt, dX_Adt = dX_Adt, dY_Adt = dY_Adt,
  constraint = c("X_C>=0", "Y_C>=0", "X_A>=0", "Y_A>=0")
)
eode_get_sysmat(x, value = pp(list(X_A = 4, Y_A = 4, X_C = 4, Y_C = 4)), delta = 10e-6)

}
