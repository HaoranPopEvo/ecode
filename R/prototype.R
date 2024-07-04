# CLASS `eode`--------------------
#' Create an ODE System
#'
#' Creates an object of class "\code{eode}" representing an ODE system that
#' describes population or community dynamics.
#' @param ... Objects of class "\code{function}", each representing a single
#' differential equation.
#' @param constraint Character strings that must have a format of
#' "\code{variable>value}" or "\code{variable<value}", such as "x>1", "y>3", etc.
#' If not specified, then all model variables are considered as positive real
#' numbers by default.
#'
#' @return An object of class "\code{eode}" describing population or community
#' dynamics.
#' @import rlang
#' @export
#' @examples
#' ## Example1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' x
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
#' x
eode <- function(..., constraint = NULL) {
  lifun <- list(...)
  all_model_variables <- c()
  all_model_parameters <- list()

  for (fun in lifun) {
    arglist <- fn_fmls(fun)
    model_variables <- names(arglist)[do.call(c, lapply(arglist, function(xx) xx == ""))]
    model_parameters <- arglist[do.call(c, lapply(arglist, function(xx) xx != ""))]

    test_code <- paste0("fun(", paste0(model_variables, "=0", collapse = ","), ")")
    error_message <- as.character(tryCatch(eval(parse(text = test_code)), error = function(e) e))
    if (length(error_message) == 0) error_message <- ""
    if (grepl("Error", error_message) & grepl("not found", error_message) & grepl("Error", error_message)) {
      invalid_variable <- strsplit(error_message, "'")[[1]][2]
      # print(fun)
      stop(paste0("variable '", invalid_variable, "' not found in function arguments. Consider the syntax of the function printed above.\n"))
    }

    all_model_variables <- c(all_model_variables, model_variables)
    all_model_parameters <- c(all_model_parameters, model_parameters)
  }

  all_model_variables <- unique(all_model_variables)
  for (single_variable in all_model_variables) {
    if (!paste0("d", single_variable, "dt") %in% names(lifun)) {
      stop("model variable '", single_variable, "' found, but equation '", paste0("d", single_variable, "dt"), "' not specified")
    }
  }

  for (parameter_name in names(all_model_parameters)) {
    single_model_parameters <- all_model_parameters[names(all_model_parameters) == parameter_name]
    if (length(single_model_parameters) >= 2) {
      if (length(unique(do.call(c, single_model_parameters))) > 1) stop(paste0("model parameter '", unique(names(single_model_parameters)), "' has multiple values"))
      all_model_parameters <- all_model_parameters[-which(names(all_model_parameters) == parameter_name)[-1]]
    }
  }

  if (is.null(constraint)) {
    constraint <- c(paste0(all_model_variables, ">0"), paste0(all_model_variables, "<1000"))
  } else {
    if (!is.character(constraint)) stop("invalid argument 'constraint'")
    if (!all(grepl(">", constraint) | grepl("<", constraint))) stop("invalid argument 'constraint'")
    var_constrainted <- c()
    for (con in constraint) {
      if (length(strsplit(con, ">")[[1]]) == 2) {
        var_constrainted <- c(var_constrainted, strsplit(con, ">")[[1]][1])
      } else if (length(strsplit(con, "<")[[1]]) == 2) {
        var_constrainted <- c(var_constrainted, strsplit(con, "<")[[1]][1])
      } else {
        stop("invalid argument 'constraint'")
      }
    }
    if (!all(all_model_variables %in% var_constrainted)) {
      constraint <- c(constraint, paste0(all_model_variables[!(all_model_variables %in% var_constrainted)], ">0"))
    }

    for (con in constraint) {
      valid_con <- c(paste0(all_model_variables, ">"), paste0(all_model_variables, "<"))
      if (!any(unlist(lapply(valid_con, grepl, con)))) {
        stop("invalid constraint '", con, "' found")
      }
    }

    for (var in all_model_variables) {
      valid_con <- c(paste0(var, ">"), paste0(var, "<"))
      if (sum(grepl(valid_con[1], constraint)) > 1) {
        wordy_index <- which(grepl(valid_con[1], constraint))
        split_index <- strsplit(constraint[wordy_index], ">")
        max_value <- max(as.numeric(lapply(split_index, function(xx) xx[2])))
        constraint <- c(constraint[-wordy_index], paste0(var, ">", max_value))
      }
      if (sum(grepl(valid_con[2], constraint)) > 1) {
        wordy_index <- which(grepl(valid_con[2], constraint))
        split_index <- strsplit(constraint[wordy_index], "<")
        min_value <- min(as.numeric(lapply(split_index, function(xx) xx[2])))
        constraint <- c(constraint[-wordy_index], paste0(var, "<", min_value))
      }
    }

    for (var in all_model_variables) {
      if (!any(grepl(paste0(var, ">"), constraint))) constraint <- c(constraint, paste0(var, ">0"))
    }
    for (var in all_model_variables) {
      if (!any(grepl(paste0(var, "<"), constraint))) constraint <- c(constraint, paste0(var, "<1000"))
    }
  }

  ret <- list(
    system = lifun,
    variables = all_model_variables,
    parameter = all_model_parameters,
    constraint = sort(constraint)
  )
  class(ret) <- "eode"
  ret
}

## S3 FUN--------------

#' Print Brief Details of an ODE System
#'
#' Prints a very brief description of an ODE system.
#' @param x Object of class "\code{eode}".
#' @param ... Ignored.
#' @return No value
#' @exportS3Method
#' @method print eode
#' @examples
#' ## Example 1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' x
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
#' x
print.eode <- function(x, ...) {
  cat(paste("An ODE system:", length(x$system), "equations\n"))
  cat("equations:\n")
  for (i in 1:length(x$system)) {
    fun_text <- paste(x$system)[[i]]
    RHS_text <- strsplit(fun_text, "\n")[[1]][2]
    LHS_text <- names(x$system)[i]
    cat(paste(" ", LHS_text, "=", RHS_text, "\n"))
  }
  cat("variables:", x$variables, "\n")
  cat("parameters:", paste(names(x$parameter), "=", x$parameter, collapse = ", "), "\n")
  cat("constraints:", x$constraint)
}

#' Length of an ODE System
#'
#' Get the number of equations in an ODE system.
#' @param x Object of "\code{eode}" class representing an ODE system under consideration.
#' @param ... Ignored.
#' @exportS3Method
#' @method length eode
#' @return The number of equations in the ODE system.
#' @examples
#' ## Example 1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' length(x)
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
#' length(x)
length.eode <- function(x, ...) {
  length(x$system)
}


#' Plot Phase Velocity Vector Field
#'
#' Creates a plot of phase velocity vector (or its field) of an ODE system. With
#' two free variables the function creates a plot of the phase velocity vector field
#' within the range defined by model constraints. With a single free variable the
#' function creates a plot of one-dimensional phase velocity vector field. With
#' more than two variables the function will automatically set a value (the middle
#' of the lower and upper boundaries) for any extra variable to reduce axes. The
#' values can also be set using parameter "\code{set_covar}". If all model variables
#' have had their values, the function creates a plot to show the exact values of a
#' single phase velocity vector.
#'
#' @param x The ODE system under consideration. An object of class "\code{eode}".
#' @param n.grid number of grids per dimension. A larger number indicates a smaller
#' grid size. Default 20.
#' @param scaleL scale factor of the arrow length.
#' @param scaleAH scale factor of the arrow head.
#' @param scaleS scale factor of the arrow size.
#' @param set_covar give certain values of model variables that are going to be fixed.
#' @param ... ignored arguments
#'
#' @import ggplot2
#' @import stringr
#' @return a graphic object
#' @exportS3Method
#' @method plot eode
#' @examples
#' ## Example 1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' plot(x)
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
#' plot(x)
#' plot(x, set_covar = list(X_A = 60, Y_A = 80))
#' plot(x, set_covar = list(Y_C = 80, X_A = 60, Y_A = 80))
#' plot(x, set_covar = list(X_C = 60, Y_C = 80, X_A = 60, Y_A = 80))
#'
plot.eode <- function(x, n.grid = 20, scaleL = 1, scaleAH = 1, scaleS = 1,
                      set_covar = NULL, ...) {
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

  if (!is.null(set_covar)) {
    if (any(!(names(set_covar) %in% x$variables))) warning(paste0(paste(names(set_covar)[!(names(set_covar) %in% x$variables)], collapse = ", "), " not used in the ODE system."))
    set_covar <- set_covar[names(set_covar) %in% x$variables]
    if (length(set_covar) == 0) set_covar <- NULL
  }

  axis_var_names <- x$variables[!(x$variables %in% names(set_covar))]
  if (length(axis_var_names) > 2) {
    add_covar <- lapply(axis_var_names[-c(1:2)], function(cov_name) {
      mean(as.numeric(unlist(str_extract_all(x$constraint[grepl(cov_name, x$constraint)], "\\d+\\.?\\d*"))))
    })
    names(add_covar) <- axis_var_names[-c(1:2)]
    set_covar <- c(set_covar, add_covar)
    axis_var_names <- axis_var_names[1:2]
    message(paste0("Set ", paste(names(add_covar), as.numeric(add_covar), sep = " = ", collapse = ", "), " for mapping in two axis\n"))
  }

  if (length(axis_var_names) == 2) {
    x_lowerbound <- min(as.numeric(str_extract_all(x$constraint[grepl(axis_var_names[1], x$constraint)], "\\d+\\.?\\d*")))
    x_upperbound <- max(as.numeric(str_extract_all(x$constraint[grepl(axis_var_names[1], x$constraint)], "\\d+\\.?\\d*")))
    y_lowerbound <- min(as.numeric(str_extract_all(x$constraint[grepl(axis_var_names[2], x$constraint)], "\\d+\\.?\\d*")))
    y_upperbound <- max(as.numeric(str_extract_all(x$constraint[grepl(axis_var_names[2], x$constraint)], "\\d+\\.?\\d*")))
    inteval_x <- (x_upperbound - x_lowerbound) / n.grid
    inteval_y <- (y_upperbound - y_lowerbound) / n.grid
    xrange <- seq(x_lowerbound, x_upperbound, inteval_x)
    yrange <- seq(y_lowerbound, y_upperbound, inteval_y)

    dxdt <- x$system[[paste0("d", axis_var_names[1], "dt")]]
    dydt <- x$system[[paste0("d", axis_var_names[2], "dt")]]
    if (!is.null(set_covar)) {
      for (i in 1:length(set_covar)) {
        dxdt <- eval(parse(text = gsub(paste0(names(set_covar[i]), ","), paste0(names(set_covar[i]), "=", set_covar[[i]], ","), paste(trimws(deparse(dxdt)), collapse = " "))))
        dydt <- eval(parse(text = gsub(paste0(names(set_covar[i]), ","), paste0(names(set_covar[i]), "=", set_covar[[i]], ","), paste(trimws(deparse(dydt)), collapse = " "))))
      }
      for (i in 1:length(axis_var_names)) {
        dxdt_string <- paste(trimws(deparse(dxdt)), collapse = " ")
        dydt_string <- paste(trimws(deparse(dydt)), collapse = " ")

        if (!grepl(axis_var_names[i], dxdt_string)) {
          dxdt_string <- paste0("function (", axis_var_names[i], ", ", strsplit(dxdt_string, "function \\(")[[1]][2])
          dxdt <- eval(parse(text = dxdt_string))
        }
        if (!grepl(axis_var_names[i], dydt_string)) {
          dydt_string <- paste0("function (", axis_var_names[i], ", ", strsplit(dydt_string, "function \\(")[[1]][2])
          dydt <- eval(parse(text = dydt_string))
        }
      }
    }
    arrows <- do.call(rbind, apply(expand.grid(xrange, yrange), MARGIN = 1, function(xx) {
      x_ <- xx[1] # start_point.x
      y_ <- xx[2] # start_point.y
      dx <- eval(parse(text = paste0("dxdt(", axis_var_names[1], "=x_, ", axis_var_names[2], "=y_)")))
      # arrow_length.x
      dy <- eval(parse(text = paste0("dydt(", axis_var_names[1], "=x_, ", axis_var_names[2], "=y_)")))
      # arrow_length.y
      data.frame(x_ = x_, y_ = y_, dx = dx, dy = dy, r = sqrt(dx * dx + dy * dy))
    }))

    rescaleL_x <- inteval_x * scaleL * 0.95 / max(abs(arrows$dx)) # scale_arrow_length.x
    rescaleL_y <- inteval_y * scaleL * 0.95 / max(abs(arrows$dy)) # scale_arrow_length.y
    arrows$x_end_ <- with(arrows, rescaleL_x * dx + x_)
    arrows$y_end_ <- with(arrows, rescaleL_y * dy + y_)

    max_arrowhead <- 0.07 * scaleAH
    arrowhead <- arrows$r / max(arrows$r) * max_arrowhead # scale_arrow_head

    arrowsize <- with(arrows, r / max(r) * 0.7 * scaleS) # scale_arrow_size

    ggplot() +
      geom_segment(
        data = arrows, mapping = aes(x = x_, y = y_, xend = x_end_, yend = y_end_, size = arrowsize),
        arrow = arrow(length = unit(arrowhead, "inches"), type = "closed", angle = 20)
      ) +
      scale_size_identity() +
      xlab(axis_var_names[1]) +
      ylab(axis_var_names[2]) +
      theme1
  } else if (length(axis_var_names) == 1) {
    x_lowerbound <- min(as.numeric(str_extract_all(x$constraint[grepl(axis_var_names, x$constraint)], "\\d+\\.?\\d*")))
    x_upperbound <- max(as.numeric(str_extract_all(x$constraint[grepl(axis_var_names, x$constraint)], "\\d+\\.?\\d*")))
    inteval_x <- (x_upperbound - x_lowerbound) / n.grid
    xrange <- seq(x_lowerbound, x_upperbound, inteval_x)

    dxdt <- x$system[[paste0("d", axis_var_names, "dt")]]
    if (!is.null(set_covar)) {
      for (i in 1:length(set_covar)) {
        dxdt <- eval(parse(text = gsub(paste0(names(set_covar[i]), ","), paste0(names(set_covar[i]), "=", set_covar[[i]], ","), paste(trimws(deparse(dxdt)), collapse = " "))))
      }
    }
    arrows <- do.call(rbind, lapply(xrange, function(xx) {
      x_ <- xx # start_point.x
      dx <- eval(parse(text = paste0("dxdt(", axis_var_names[1], "=x_)"))) # arrow_length.x
      data.frame(x_ = x_, dx = dx, r = abs(dx))
    }))

    rescaleL_x <- inteval_x * scaleL * 0.95 / max(abs(arrows$dx)) # scale_arrow_length.x
    arrows$x_end_ <- with(arrows, rescaleL_x * dx + x_)

    max_arrowhead <- 0.07 * scaleAH
    arrowhead <- arrows$r / max(arrows$r) * max_arrowhead # scale_arrow_head
    arrowsize <- with(arrows, r / max(r) * 0.7 * scaleS) # scale_arrow_size

    arrows$y_ <- 0
    arrows$y_end_ <- 0

    ggplot() +
      geom_segment(
        data = arrows, mapping = aes(x = x_, y = y_, xend = x_end_, yend = y_end_, size = arrowsize),
        arrow = arrow(length = unit(arrowhead, "inches"), type = "closed", angle = 20)
      ) +
      scale_size_identity() +
      xlab(axis_var_names) +
      ylab(NULL) +
      theme1 +
      theme(
        axis.text.y = element_blank(),
        axis.line = element_line(size = 1.2),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()
      )
  } else if (length(axis_var_names) == 0) {
    arrows <- lapply(x$system, function(fun) {
      do.call(fun, set_covar[names(set_covar) %in% names(formals(fun))])
    })
    arrows <- data.frame(
      name = names(arrows),
      y_end_ = as.numeric(unlist(arrows)),
      y_ = 0,
      x_ = 1:length(arrows),
      x_end_ = 1:length(arrows)
    )

    ggplot() +
      geom_segment(
        data = arrows, mapping = aes(x = x_, y = y_, xend = x_end_, yend = y_end_), size = scaleS,
        arrow = arrow(length = unit(scaleAH * 0.15, "inches"), type = "closed", angle = 20)
      ) +
      geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed", color = "blue") +
      xlab(NULL) +
      ylab("velocity") +
      scale_x_continuous(labels = arrows$name) +
      theme1 +
      theme(axis.text.x = element_text(hjust = 0.8))
  } else {
    stop("incorrect length of axes.")
  }
}

## NON-S3------------------

#' Set New Constraints
#'
#' Set new constraints for an ODE system.
#' @param x an object of class "\code{eode}".
#' @param new_constraint a vector of characters indicating new constraints
#' to be assigned to the ODE system.
#'
#' @import stringr
#' @return an object of "\code{eode}" class
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' x
#' eode_set_constraint(x, c("x<5", "y<5"))
eode_set_constraint <- function(x, new_constraint) {
  for (con in new_constraint) {
    if (!any(unlist(lapply(c(">", "<"), grepl, con)))) stop(paste0("invalid constraint '", con, "' detected"))
    con <- gsub("\\s+", "", con)
    model_cons <- gsub("\\s+", "", x$constraint)
    if (grepl("<", con)) {
      new_upper_bound <- as.numeric(str_extract_all(con, "\\d+")[[1]])
      var_name <- strsplit(con, "<")[[1]][1]
      old_lower_bound <- as.numeric(str_extract_all(
        model_cons[grepl(">", model_cons) & grepl(var_name, model_cons)],
        "\\d+"
      )[[1]])
      if (new_upper_bound <= old_lower_bound) stop(paste0("conflict between new constraint '", con, "' and old constraint '", model_cons[grepl(">", model_cons) & grepl(var_name, model_cons)], "'"))
      x$constraint[grepl("<", model_cons) & grepl(var_name, model_cons)] <- con
    } else if (grepl(">", con)) {
      new_lower_bound <- as.numeric(str_extract_all(con, "\\d+")[[1]])
      var_name <- strsplit(con, ">")[[1]][1]
      old_upper_bound <- as.numeric(str_extract_all(
        model_cons[grepl("<", model_cons) & grepl(var_name, model_cons)],
        "\\d+"
      )[[1]])
      if (new_lower_bound >= old_upper_bound) stop(paste0("conflict between new constraint '", con, "' and old constraint '", model_cons[grepl("<", model_cons) & grepl(var_name, model_cons)], "'"))
      x$constraint[grepl(">", model_cons) & grepl(var_name, model_cons)] <- con
    }
  }
  x
}

#' Set New Parameters
#'
#' Set new parameters for an ODE system.
#' @param x an object of class "\code{eode}".
#' @param ParaList a list of names of parameters with their corresponding values
#' to be assigned to the ODE system.
#'
#' @return an object of "\code{eode}" class
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' x
#' eode_set_parameter(x, ParaList = list(r1 = 3))
eode_set_parameter <- function(x, ParaList) {
  if (length(ParaList) == 0) {
    return(x)
  }
  for (i in 1:length(ParaList)) {
    var_name <- names(ParaList)[i]
    var_value <- ParaList[[i]]
    if (!any(names(x$parameter) == var_name)) {
      warning(paste0("parameter '", var_name, "' does not exist."))
      next
    }
    for (j in 1:length(x$system)) {
      fun <- x$system[[j]]
      if (any(grepl(paste0(var_name, " = "), deparse(fun)))) {
        x$system[[j]] <- eval(parse(text = gsub(
          paste0(var_name, " = \\d+(\\.\\d+)?"),
          paste0(var_name, " = ", var_value), deparse(fun)
        )))
      }
    }
    x$parameter[[which(names(x$parameter) == var_name)]] <- var_value
  }
  x
}

# CLASS `pp`-------------------
#' Create a Phase Point
#'
#' Creates an object of class "\code{pp}" representing a phase point in the phase
#' space.
#' @param value A list of several elements, each giving a value on one dimension.
#'
#' @return An object of class "\code{pp}" describing a phase point.
#' @export
#'
#' @examples
#' pp(list(x = 1, y = 1))
pp <- function(value) {
  if (!inherits(value, "list")) stop("argument 'value' should be an object of class 'list'")
  if (max(as.numeric(lapply(value, length))) > 1) stop("each component of the list 'value' should be a scalar")
  if (!all(as.logical(lapply(value, is.numeric)))) stop("non-numerical values found")
  if (any(names(value) == "")) stop("component of the list 'value' should be labelled with its name")
  class(value) <- "pp"
  value
}

#' Print Brief Details of a Phase Point
#'
#' Prints a very brief description of a phase point.
#' @param x Object of class "\code{pp}".
#' @param ... Ignored.
#' @return No value
#' @export
#' @examples
#' pp(list(x = 1, y = 1))
print.pp <- function(x, ...) {
  cat("phase point:\n")
  for (i in 1:length(x)) {
    cat(names(x)[i], "=", x[[i]], "\n")
  }
}

#' Add Operator for Phase Points
#'
#' The add operator used to perform arithmetic on phase points.
#' @param x an object of \code{"pp"} class representing a phase point.
#' @param y an object of \code{"pp"} class representing another phase point.
#'
#' @return an object of \code{"pp"} class representing the calculated result.
#' @export
#'
#' @examples
#' pp(list(x = 1, y = 1)) + pp(list(x = 3, y = 4))
"+.pp" <- function(x, y) {
  if (!(all(names(x) %in% names(y)) & all(names(y) %in% names(x)))) stop("phase point 'x' and 'y' are not in the same phase space")
  for (i in 1:length(x)) {
    x[[i]] <- x[[i]] + y[names(x)[i] == names(y)][[1]]
  }
  x
}

#' Subtract Operator for Phase Points
#'
#' The subtract operator used to perform arithmetic on phase points.
#' @param x an object of \code{"pp"} class representing a phase point.
#' @param y an object of \code{"pp"} class representing another phase point.
#'
#' @return an object of \code{"pp"} class representing the calculated result.
#' @export
#'
#' @examples
#' pp(list(x = 1, y = 1)) - pp(list(x = 3, y = 4))
"-.pp" <- function(x, y) {
  if (!(all(names(x) %in% names(y)) & all(names(y) %in% names(x)))) stop("phase point 'x' and 'y' are not in the same phase space")
  for (i in 1:length(x)) {
    x[[i]] <- x[[i]] - y[names(x)[i] == names(y)][[1]]
  }
  x
}

#' Multiply Operator for Phase Points
#'
#' The multiply operator used to perform arithmetic on phase points.
#' @param x an object of \code{"pp"} class representing a phase point.
#' @param y an object of \code{"pp"} class representing another phase point.
#'
#' @return an object of \code{"pp"} class representing the calculated result.
#' @export
#'
#' @examples
#' pp(list(x = 1, y = 1)) * pp(list(x = 3, y = 4))
"*.pp" <- function(x, y) {
  if (!(all(names(x) %in% names(y)) & all(names(y) %in% names(x)))) stop("phase point 'x' and 'y' are not in the same phase space")
  for (i in 1:length(x)) {
    x[[i]] <- x[[i]] * y[names(x)[i] == names(y)][[1]]
  }
  x
}

#' Divide Operator for Phase Points
#'
#' The devide operator used to perform arithmetic on phase points.
#' @param x an object of \code{"pp"} class representing a phase point.
#' @param y an object of \code{"pp"} class representing another phase point.
#'
#' @return an object of \code{"pp"} class representing the calculated result.
#' @export
#'
#' @examples
#' pp(list(x = 1, y = 1)) / pp(list(x = 3, y = 4))
"/.pp" <- function(x, y) {
  if (!(all(names(x) %in% names(y)) & all(names(y) %in% names(x)))) stop("phase point 'x' and 'y' are not in the same phase space")
  for (i in 1:length(x)) {
    x[[i]] <- x[[i]] / y[names(x)[i] == names(y)][[1]]
  }
  x
}

#' Distance between Phase Points
#'
#' Computes and returns the distance between two phase points.
#' @param x an object of \code{"pp"} class representing a phase point.
#' @param y an object of \code{"pp"} class representing another phase point.
#'
#' @return A scalar.
#' @export
#'
#' @examples
#' a <- pp(list(x = 1, y = 1))
#' b <- pp(list(x = 3, y = 4))
#' pdist(a, b)
pdist <- function(x, y) {
  if (!(all(names(x) %in% names(y)) & all(names(y) %in% names(x)))) stop("phase point 'x' and 'y' are not in the same phase space")
  distance <- 0
  for (i in 1:length(x)) {
    distance <- distance + (x[[i]] - y[names(x)[i] == names(y)][[1]])^2
  }
  sqrt(distance)
}


# CLASS `pdata`-------------------
#' Create Population Dynamics Data
#'
#' Create a new data set to describe population dynamics with initial conditions.
#' The data is mainly used to train the ODE model described by an object of "\code{eode}"
#' class.
#' @param x an object of class "\code{eode}" describing the ODE system in focus.
#' @param init a data frame describing the initial conditions of the population.
#' Each column should correspond to a variable in the ODE system, hence the name of
#' column will be automatically checked to make sure it matches a variable in the ODE
#' system. Each row represents an independent observation
#' @param t a numerical vector that describes the time for each observation
#' @param lambda a data frame describing the observation values of population
#' dynamics. Each column is an variable and each row represents an independent
#' observation.
#' @param formula a character that specifies how the observed variable can be
#' predicted by the ODE system.
#'
#' @return an object of "\code{pdata}" class
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' pdata(x,
#'   init = data.frame(x = c(10, 20), y = c(5, 15)),
#'   t = c(3, 3),
#'   lambda = data.frame(z = c(15, 30)),
#'   formula = "z = x + y"
#' )
pdata <- function(x, init, t, lambda, formula) {
  if (!is.data.frame(init)) stop("'init' should be a data frame.")
  if (!is.data.frame(lambda)) stop("'lambda' should be a data frame.")
  if (nrow(init) != nrow(lambda) || nrow(init) != length(t) || length(t) != nrow(lambda)) {
    stop("unequal length of input data. Cannot decide the number of observations.")
  }
  if (ncol(lambda) != length(formula)) {
    stop(paste0(ncol(lambda), " observed variables in 'lambda' but ", length(formula), " formula", ifelse(length(formula) > 1, "s", ""), " found."))
  }
  if (!all(colnames(init) %in% x$variables)) {
    stop(paste0("variable '", colnames(init)[!(colnames(init) %in% x$variables)][1], "' is not found in the ODE system."))
  }
  if (!all(x$variables %in% colnames(init))) {
    stop(paste0("the initial value of variable '", x$variables[!(x$variables %in% colnames(init))][1], "' in the ODE system is not specified."))
  }
  ret <- list(init = init, t = t, lambda = lambda, formula = formula)
  class(ret) <- "pdata"
  attr(ret, "ode") <- x
  ret
}

#' Print Brief Details of Population Dynamics Data
#'
#' Prints a very brief description of population dynamics data.
#' @param x Object of class "\code{pdata}".
#' @param ... Ignored.
#' @return No value
#' @exportS3Method
#' @method print pdata
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' pdata(x,
#'   init = data.frame(x = c(10, 20), y = c(5, 15)),
#'   t = c(3, 3),
#'   lambda = data.frame(z = c(15, 30)),
#'   formula = "z = x + y"
#' )
print.pdata <- function(x, ...) {
  cat("Population Dynamics Data.\n")
  colnames(x$init) <- paste0(colnames(x$init), "_0")
  print(data.frame(x$init, observation_time = x$t, x$lambda))
  cat(paste0("Formula", ifelse(length(x$formula) > 1, "s", ""), ":\n"))
  for (i in 1:length(x$formula)) {
    cat(x$formula[i], "\n")
  }
}

#' Get Length Of Population Dynamics Data
#'
#' Get the number of observations of population dynamics data.
#' @param x Object of class "\code{pdata}".
#' @param ... Ignored.
#' @return A scalar.
#' @exportS3Method
#' @method length pdata
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' pdata(x,
#'   init = data.frame(x = c(10, 20), y = c(5, 15)),
#'   t = c(3, 3),
#'   lambda = data.frame(z = c(15, 30)),
#'   formula = "z = x + y"
#' )
length.pdata <- function(x, ...) {
  nrow(x$init)
}

#' Extract Parts of an Population Dynamics Data
#'
#' Operators acting on population dynamics data.
#' @param x Object of class "\code{pdata}".
#' @param i indice specifying elements to extract.
#' @return pdata object
#' @export
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2)
#' pdat <- pdata(x,
#'   init = data.frame(x = c(10, 20), y = c(5, 15)),
#'   t = c(3, 3),
#'   lambda = data.frame(z = c(15, 30)),
#'   formula = "z = x + y"
#' )
#' pdat[1]
"[.pdata" <- function(x, i) {
  pdat <- x # change variable name
  init <- pdat$init[i, ]
  if (!is.data.frame(init)) {
    init
  }
  lambda <- pdat$lambda[i, ]
  if (!is.data.frame(lambda)) {
    lambda <- data.frame(lambda)
    colnames(lambda) <- colnames(pdat$lambda)
  }

  pdata(attr(pdat, "ode"), init, t = pdat$t[i], lambda, formula = pdat$formula)
}
