# VELOCITY---------------
#' Velocity Vector
#'
#' Calculates the velocity vector at a given phase point.
#' @param x The ODE system under consideration. An object of class "\code{eode}".
#' @param value an object of "\code{pp}" classs representing a phase point in the
#' ODE system under consideration.
#' @param phase_curve an object of "\code{pc}" class representing a phase curve.
#' Only used when there is delayed items in the ODE system.
#'
#' @import stringr
#' @return an object of class "\code{velocity}" representing a velocity vector at
#' a given phase point.
#' @export
#'
#' @examples
#'
#' ## Example1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_get_velocity(x, value = pp(list(x = 1, y = 1)))
#'
#' ## Example2: Susceptible-infected model
#' dX_Cdt <- function(X_C, Y_C, X_A, Y_A, nu = 0.15, beta = 0.1, mu = 0.15, g = 0.04) {
#'   nu * (X_A + Y_A) - beta * X_C * (Y_C + Y_A) - (mu + g) * X_C
#' }
#'
#' dY_Cdt <- function(X_C, Y_C, Y_A, beta = 0.1, mu = 0.15, g = 0.04, rho = 0.2) {
#'   beta * X_C * (Y_C + Y_A) - (mu + g + rho) * Y_C
#' }
#'
#' dX_Adt <- function(X_C, Y_C, X_A, Y_A, beta = 0.1, g = 0.04) {
#'   g * X_C - beta * X_A * (Y_C + Y_A)
#' }
#'
#' dY_Adt <- function(X_A, Y_C, Y_A, beta = 0.1, g = 0.04, rho = 0.2) {
#'   beta * X_A * (Y_C + Y_A) + g * Y_C - rho * Y_A
#' }
#'
#' x <- eode(
#'   dX_Cdt = dX_Cdt, dY_Cdt = dY_Cdt, dX_Adt = dX_Adt, dY_Adt = dY_Adt,
#'   constraint = c("X_C>=0", "Y_C>=0", "X_A>=0", "Y_A>=0")
#' )
#' eode_get_velocity(x, value = pp(list(X_A = 5, Y_A = 5, X_C = 3, Y_C = 2)))
#'
eode_get_velocity <- function(x, value, phase_curve = NULL) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (any(grepl("\\[", as.character(unlist(lapply(x$system, deparse))))) && is.null(phase_curve)) stop("delayed items found but 'phase_curve' not known. Cannot calculate it.")

  velo_list <- lapply(x$system, function(fun) {
    input_args <- fn_fmls(fun)
    vars_need_value <- names(input_args)[as.character(input_args) == ""]
    if (any(grepl("\\[", deparse(fun)))) {
      for (model_var in x$variables) {
        if (any(grepl(paste0(model_var, "\\[-.+\\]"), deparse(fun)))) {
          all_delayed_model_vars <- do.call(c, str_extract_all(deparse(fun), paste0(model_var, "\\[-.+\\]")))
          for (delayed_model_var in all_delayed_model_vars) {
            delayed_time <- eval(parse(text = paste0("with(input_args,", strsplit(strsplit(delayed_model_var, "\\-")[[1]][2], "\\]")[[1]], ")")))
            delayed_step <- delayed_time / phase_curve$step
            past_series <- phase_curve[[which(names(phase_curve) == model_var)]]
            if (length(past_series) - delayed_step <= 0) {
              delayed_value <- 0
            } else {
              delayed_value <- past_series[length(past_series) - delayed_step]
            }
            fun <- eval(parse(text = gsub(str_escape(delayed_model_var), delayed_value, deparse(fun))))
          }
        }
      }
    }
    eval(parse(text = paste0("fun(", paste0(do.call(c, lapply(vars_need_value, function(var) {
      paste0(var, "=", as.numeric(value[var == names(value)]))
    })), collapse = ","), ")")))
  })
  ret <- t(t(as.numeric(velo_list)))
  row.names(ret) <- names(velo_list)
  colnames(ret) <- "value"
  class(ret) <- "velocity"
  ret
}


#' Test the Validity of a Phase Point
#'
#' Test whether a phase point sits out of the boundary of an ODE system.
#' @param x object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#'
#' @return returns \code{TRUE} or \code{FALSE} depending on whether the values of
#' the system variables are within the boundary or not.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_validval(x, value = pp(list(x = 1, y = 1)))
eode_is_validval <- function(x, value) {
  if (!inherits(value, "pp")) stop("'value' should be an object of class 'pp'")
  if (!all(x$variables %in% names(value))) stop("the value of the variable '", x$variables[!(x$variables %in% names(value))][1], "' not found")
  is_valid <- TRUE
  for (con in x$constraint) {
    if (!eval(parse(text = paste0("value$", con)))) is_valid <- FALSE
  }
  is_valid
}

#' Print Brief Details of a Phase Velocity Vector
#'
#' Prints a very brief description of a phase velocity vector.
#' @param x Object of class "\code{velocity}".
#' @param ... Ignored.
#' @return No value
#' @exportS3Method
#' @method print velocity
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_get_velocity(x, value = pp(list(x = 1, y = 1)))
print.velocity <- function(x, ...) {
  cat("phase velocity vector:\n")
  for (i in 1:nrow(x)) {
    cat(row.names(x)[i], "=", x[i, ], "\n")
  }
}

#' Create a Plot of a Phase Velocity Vector
#'
#' Plot a phase velocity vector
#' @param x Object of class "\code{velocity}".
#' @param ... Ignored.
#'
#' @import ggplot2
#' @return ggplot object
#' @exportS3Method
#' @method plot velocity
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' plot(eode_get_velocity(x, value = pp(list(x = 1, y = 1))))
plot.velocity <- function(x, ...) {
  velo <- x # change variable names
  x_ <- y_ <- x_end_ <- y_end_ <- NULL
  theme1 <- theme_bw() +
    theme(
      axis.text.x = element_text(size = 16, angle = 0, colour = "black"),
      axis.text.y = element_text(size = 16, angle = 0, colour = "black"),
      axis.title = element_text(size = 18),
      axis.line = element_line(linetype = 1, color = "black", linewidth = 0.1),
      axis.ticks = element_line(colour = "black"),
      panel.grid.major = element_blank(), # change the major and minor grid lines,
      panel.grid.minor = element_blank(), # if want to change, check this parameters, I think it's easier to dao that
      # strip.background = element_rect(colour = "black",size = 0.8),
      # panel.background = element_rect(colour="black", fill="white"),
      # panel.border = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      plot.title = element_text(size = 14, angle = 0, colour = "black", face = "italic"),
      plot.tag = element_text(size = 14, angle = 0, colour = "black", face = "bold"),
      plot.caption = element_text(size = 14, angle = 0, colour = "black", face = "italic"),
      axis.title.y = element_text(vjust = 1.9),
      axis.title.x = element_text(vjust = 0.5),
      legend.text = element_text(colour = "black", size = 14),
      legend.background = element_rect(fill = "transparent", color = NA),
      # legend.position = "bottom",
      legend.title = element_blank()
    )

  arrows <- data.frame(
    name = row.names(velo),
    y_end_ = as.numeric(velo),
    y_ = 0,
    x_ = 1:length(velo),
    x_end_ = 1:length(velo)
  )

  ggplot() +
    geom_segment(
      data = arrows, mapping = aes(x = x_, y = y_, xend = x_end_, yend = y_end_), size = 1,
      arrow = arrow(length = unit(1 * 0.15, "inches"), type = "closed", angle = 20)
    ) +
    geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed", color = "blue") +
    xlab(NULL) +
    ylab("velocity") +
    scale_x_continuous(labels = arrows$name, breaks = c(1, 2)) +
    theme1 +
    theme(axis.text.x = element_text(hjust = 0.8))
}

# EQUILBRIUM-----------
#' System Matrix
#'
#' Linearlise an ODE system and calculate its system matrix.
#' @param x object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param delta Spacing or step size of a numerical grid used for calculating
#' numerical differentiation.
#'
#' @return An object of class "\code{matrix}". Each matrix entry, (i,j), represents
#' the partial derivative of the i-th equation of the system with respect to the
#' j-th variable taking other variables as a constant.
#' @export
#'
#' @examples
#' ## Example1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_get_sysmat(x, value = pp(list(x = 1, y = 1)), delta = 10e-6)
#'
#' ## Example2: Susceptible-infected model
#' dX_Cdt <- function(X_C, Y_C, X_A, Y_A, nu = 0.15, beta = 0.1, mu = 0.15, g = 0.04) {
#'   nu * (X_A + Y_A) - beta * X_C * (Y_C + Y_A) - (mu + g) * X_C
#' }
#'
#' dY_Cdt <- function(X_C, Y_C, Y_A, beta = 0.1, mu = 0.15, g = 0.04, rho = 0.2) {
#'   beta * X_C * (Y_C + Y_A) - (mu + g + rho) * Y_C
#' }
#'
#' dX_Adt <- function(X_C, Y_C, X_A, Y_A, beta = 0.1, g = 0.04) {
#'   g * X_C - beta * X_A * (Y_C + Y_A)
#' }
#'
#' dY_Adt <- function(X_A, Y_C, Y_A, beta = 0.1, g = 0.04, rho = 0.2) {
#'   beta * X_A * (Y_C + Y_A) + g * Y_C - rho * Y_A
#' }
#'
#' x <- eode(
#'   dX_Cdt = dX_Cdt, dY_Cdt = dY_Cdt, dX_Adt = dX_Adt, dY_Adt = dY_Adt,
#'   constraint = c("X_C>=0", "Y_C>=0", "X_A>=0", "Y_A>=0")
#' )
#' eode_get_sysmat(x, value = pp(list(X_A = 4, Y_A = 4, X_C = 4, Y_C = 4)), delta = 10e-6)
#'
eode_get_sysmat <- function(x, value, delta = 10e-6) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  sysmat <- do.call(cbind, lapply(x$variables, function(var) {
    value_plus_det <- value
    value_plus_det[[var]] <- value_plus_det[[var]] + delta
    eode_get_velocity(x, value_plus_det) - eode_get_velocity(x, value)
  })) / delta
  colnames(sysmat) <- paste0("\u2202", x$variables)
  sysmat
}

# CRITICAL POINT-----------------------

#' Find Equilibrium Point
#'
#' Finds an equilibrium point (or critical point) of an ODE system based on
#' Newton iteration method.
#'
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param init_value An object of class "\code{pp}" representing a phase point
#' giving start estimates.
#' @param eps Precision for the stopping criterion. Iteration will stop after
#' the movement of the phase point in a single step is smaller than \code{eps}.
#' @param max_step Maximum number
#' @param method one of "Newton", "GradDesc"
#' @param verbose Logical, whether to print the iteration process.
#'
#' @return An object of class "\code{pp}" representing an equilibrium point found.
#' @export
#'
#' @import stringr
#' @examples
#' ## Example 1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 1, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_get_cripoi(x, init_value = pp(list(x = 0.5, y = 0.5)))
#'
#' ## Example 2: Susceptible-infected model
#' dX_Cdt <- function(X_C, Y_C, X_A, Y_A, nu = 0.15, beta = 0.1, mu = 0.15, g = 0.04) {
#'   nu * (X_A + Y_A) - beta * X_C * (Y_C + Y_A) - (mu + g) * X_C
#' }
#'
#' dY_Cdt <- function(X_C, Y_C, Y_A, beta = 0.1, mu = 0.15, g = 0.04, rho = 0.2) {
#'   beta * X_C * (Y_C + Y_A) - (mu + g + rho) * Y_C
#' }
#'
#' dX_Adt <- function(X_C, Y_C, X_A, Y_A, beta = 0.1, g = 0.04) {
#'   g * X_C - beta * X_A * (Y_C + Y_A)
#' }
#'
#' dY_Adt <- function(X_A, Y_C, Y_A, beta = 0.1, g = 0.04, rho = 0.2) {
#'   beta * X_A * (Y_C + Y_A) + g * Y_C - rho * Y_A
#' }
#'
#' x <- eode(
#'   dX_Cdt = dX_Cdt, dY_Cdt = dY_Cdt, dX_Adt = dX_Adt, dY_Adt = dY_Adt,
#'   constraint = c("X_C>=0", "Y_C>=0", "X_A>=0", "Y_A>=0")
#' )
#' eode_get_cripoi(x, init_value = pp(list(X_C = 1, Y_C = 1, X_A = 1, Y_A = 1)))
eode_get_cripoi <- function(x, init_value, eps = 10e-4, max_step = 0.01, method = c("Newton", "GradDesc"), verbose = TRUE) {
  if (!eode_is_validval(x, init_value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (length(method) > 1) method <- method[1]

  var_range <- lapply(x$variables, function(var) {
    c(
      minval = as.numeric(str_extract_all(x$constraint[grepl(var, x$constraint) & grepl(">", x$constraint)], "\\d+\\.?\\d*")[[1]]),
      maxval = as.numeric(str_extract_all(x$constraint[grepl(var, x$constraint) & grepl("<", x$constraint)], "\\d+\\.?\\d*")[[1]])
    )
  })
  names(var_range) <- x$variables

  if (method == "Newton") {
    var_input <- init_value
    track <- list(var_input)
    while (1) {
      last_input <- var_input
      res <- tryCatch(solve(eode_get_sysmat(x, var_input), eode_get_velocity(x, var_input)),
        error = function(e) e
      )
      if ("error" %in% class(res)) {
        message("Fail to find an equilibrium point. The function will return iteration history\n")
        if (verbose) {
          print(data.frame(do.call(rbind, lapply(track, function(xx) {
            xx
          }))))
        }
        stop("Fail to find an equilibrium point. Please try other initial values")
      }
      for (ii in 1:length(var_input)) {
        var_input[[ii]] <- var_input[[ii]] - res[ii, ]
      }
      track <- c(track, list(var_input))
      if (pdist(var_input, last_input) < eps) break
    }
  } else if (method == "GradDesc") {
    scaleFactor <- max_step / sqrt(sum(eode_get_velocity(x, init_value)^2))
    var_input <- init_value
    track <- list(var_input)
    while (1) {
      last_input <- var_input
      pvv <- eode_get_velocity(x, var_input)
      move <- pvv * scaleFactor
      for (ii in 1:length(var_input)) {
        var_input[[ii]] <- var_input[[ii]] + move[ii, ]
      }
      track <- c(track, list(var_input))
      if (!eode_is_validval(x, var_input)) {
        message("Phase points go out of the boundary. Iteration History:\n")
        if (verbose) {
          print(data.frame(do.call(rbind, lapply(track, function(xx) {
            xx
          }))))
        }
        stop("Fail to find an equilibrium point. Please try other initial values")
      }
      if (sqrt(sum(pvv^2)) < eps) break
    }
  } else {
    stop("invalid parameter 'method'.")
  }
  var_input
}


# STABLE-------------------------
#' Stability Analysis
#'
#' Check whether an equilibrium point is stable or not, and give
#' the specific type (i.e. stable node, stable focus, unstable node,
#' unstable focus, saddle, or centre).
#' Check whether an equilibrium point is stable or not.
#'
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a string character indicate the type of the equilbrium point.
#' @export
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_stability_type(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_stability_type <- function(x, value, eps = 10e-4) {
  if (eode_is_stanod(x, value, eps)) {
    return("stable node")
  }
  if (eode_is_stafoc(x, value, eps)) {
    return("stable focus")
  }
  if (eode_is_unsnod(x, value, eps)) {
    return("unstable node")
  }
  if (eode_is_unsfoc(x, value, eps)) {
    return("unstable focus")
  }
  if (eode_is_saddle(x, value, eps)) {
    return("saddle")
  }
  if (eode_is_centre(x, value, eps)) {
    return("centre")
  }
}


#' Stable Equilibrium Point
#'
#' Check whether an equilibrium point is stable or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_stapoi(x, value = pp(list(x = 0.3333, y = 0.3333)))
#'
#' ## Example2: Susceptible-infected model
#' dX_Cdt <- function(X_C, Y_C, X_A, Y_A, nu = 0.15, beta = 0.1, mu = 0.15, g = 0.04) {
#'   nu * (X_A + Y_A) - beta * X_C * (Y_C + Y_A) - (mu + g) * X_C
#' }
#'
#' dY_Cdt <- function(X_C, Y_C, Y_A, beta = 0.1, mu = 0.15, g = 0.04, rho = 0.2) {
#'   beta * X_C * (Y_C + Y_A) - (mu + g + rho) * Y_C
#' }
#'
#' dX_Adt <- function(X_C, Y_C, X_A, Y_A, beta = 0.1, g = 0.04) {
#'   g * X_C - beta * X_A * (Y_C + Y_A)
#' }
#'
#' dY_Adt <- function(X_A, Y_C, Y_A, beta = 0.1, g = 0.04, rho = 0.2) {
#'   beta * X_A * (Y_C + Y_A) + g * Y_C - rho * Y_A
#' }
#'
#' x <- eode(
#'   dX_Cdt = dX_Cdt, dY_Cdt = dY_Cdt, dX_Adt = dX_Adt, dY_Adt = dY_Adt,
#'   constraint = c(
#'     "X_C>=0", "Y_C>=0", "X_A>=0", "Y_A>=0",
#'     "X_C<5", "Y_C<5", "X_A<5", "Y_A<5"
#'   )
#' )
#' eode_is_stapoi(x, value = pp(list(X_A = 1.3443, Y_A = 0.2304, X_C = 1.0655, Y_C = 0.0866)))
eode_is_stapoi <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  all(Re(eigen(eode_get_sysmat(x, value))$value) < 0)
}


#' Stable Node
#'
#' Check whether an equilibrium point is a stable node or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 2, a12 = 1) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 1, a22 = 2) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_stanod(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_is_stanod <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  eigenvalues <- eigen(eode_get_sysmat(x, value))$value
  all(Im(eigenvalues) == 0) & all(Re(eigenvalues) < 0)
}

#' Stable Focus
#'
#' Check whether an equilibrium point is a stable focus or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 2, a12 = 1) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 1, a22 = 2) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_stafoc(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_is_stafoc <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  eigenvalues <- eigen(eode_get_sysmat(x, value))$value
  all(Im(eigenvalues) != 0) & all(Re(eigenvalues) < 0)
}


# UNSTABLE----------------------
#' Unstable Node
#'
#' Check whether an equilibrium point is an unstable node or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 2, a12 = 1) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 1, a22 = 2) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_unsnod(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_is_unsnod <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  eigenvalues <- eigen(eode_get_sysmat(x, value))$value
  all(Im(eigenvalues) == 0) & all(Re(eigenvalues) > 0)
}

#' Unstable Focus
#'
#' Check whether an equilibrium point is an unstable focus or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 2, a12 = 1) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 1, a22 = 2) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_unsfoc(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_is_unsfoc <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  eigenvalues <- eigen(eode_get_sysmat(x, value))$value
  all(Im(eigenvalues) != 0) & all(Re(eigenvalues) > 0)
}

# QUASI-STABLE---------------------
#' Saddle
#'
#' Check whether an equilibrium point is a saddle or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_saddle(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_is_saddle <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  eigenvalues <- eigen(eode_get_sysmat(x, value))$value
  all(Im(eigenvalues) == 0) & any(Re(eigenvalues) > 0) & any(Re(eigenvalues) < 0)
}

#' Centre
#'
#' Check whether an equilibrium point is a neutral centre or not.
#' @param x Object of class "\code{eode}" representing an ODE system.
#' @param value an object of class "\code{pp}" representing a phase point in the
#' ODE system under consideration.
#' @param eps Precision used to check whether the input phase point is an equilibrium
#' point. If the absolute value of any component of the phase velocity vector at a
#' phase point is lower than \code{eps}, an error would be thrown out.
#'
#' @return a \code{TRUE} or \code{FALSE}.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' eode_is_centre(x, value = pp(list(x = 0.3333, y = 0.3333)))
eode_is_centre <- function(x, value, eps = 10e-4) {
  if (!eode_is_validval(x, value)) stop("phase point out of the boundary. Please check the constraints of the ODE system")
  if (!all(abs(eode_get_velocity(x, value)) < eps)) stop("'value' is not an equilibrium point")
  eigenvalues <- eigen(eode_get_sysmat(x, value))$value

  all(Im(eigenvalues) != 0) & all(Re(eigenvalues) == 0)
}
