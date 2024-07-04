# PHASE CURVE--------------
#' Solve ODEs
#'
#' Provides the numerical solution of an ODE under consideration given an initial
#' condition.
#' @param x the ODE system under consideration. An object of class "\code{eode}".
#' @param value0 value an object of class "\code{pp}" representing the phase point
#' at the initial state.
#' @param N number of iterations
#' @param step interval of time for each step
#'
#' @import stringr
#' @return an object of class "\code{pc}" that represents a phase curve
#' @export
#'
#' @examples
#' ## Example1: Lotka-Volterra competition model
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<1", "y<1"))
#' eode_proj(x, value0 = pp(list(x = 0.2, y = 0.1)), N = 100)
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
#' eode_proj(x, value0 = pp(list(X_A = 5, Y_A = 5, X_C = 3, Y_C = 2)), N = 100)
eode_proj <- function(x, value0, N, step = 0.01) {
  trajectory <- as.data.frame(t(as.numeric(value0))) # initial state
  colnames(trajectory) <- names(value0)

  value <- value0
  for (i in 1:N) {
    pc_now <- trajectory
    pc_now$t <- 0:(nrow(trajectory) - 1)
    class(pc_now) <- "pc"
    pc_now$step <- step

    movement <- eode_get_velocity(x, value, phase_curve = pc_now) * step
    for (name in names(value)) {
      value[[name]] <- value[[name]] + movement[paste0("d", name, "dt"), ]
    }
    if (!eode_is_validval(x, value)) {
      warning("Phase point reaches boundary at iteration ", i)
      # options(warn = -1)
      suppressWarnings({
        for (model_var in names(value)) {
          for (focus_constraint in x$constraint[grepl(model_var, x$constraint)]) {
            if (!eval(parse(text = paste0("with(value,", focus_constraint, ")")))) {
              if (grepl("=", focus_constraint)) {
                value[[which(names(value) == model_var)]] <- as.numeric(str_extract_all(focus_constraint, "\\d+")[[1]])
              } else {
                if (grepl(">", focus_constraint)) {
                  value[[which(names(value) == model_var)]] <- as.numeric(str_extract_all(focus_constraint, "\\d+")[[1]]) + 0.001
                } else {
                  value[[which(names(value) == model_var)]] <- as.numeric(str_extract_all(focus_constraint, "\\d+")[[1]]) - 0.001
                }
              }
            }
          }
        }
      })
    }
    traj <- as.data.frame(t(as.numeric(value)))
    colnames(traj) <- names(value)

    trajectory <- rbind(trajectory, traj)
  }
  # options(warn = 0)

  trajectory$t <- 0:(nrow(trajectory) - 1)
  class(trajectory) <- "pc"
  trajectory$step <- step
  attr(trajectory, "ode") <- x
  trajectory
}

## S3 FUN-----------

#' Print Brief Details of a Phase Curve
#'
#' Prints a very brief description of a phase curve.
#' @param x Object of class "\code{pc}".
#' @param ... Ignored.
#' @return No value
#' @exportS3Method
#' @method print pc
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<1", "y<1"))
#' eode_proj(x, value0 = pp(list(x = 0.9, y = 0.9)), N = 8)
print.pc <- function(x, ...) {
  cat("phase curve:\n")
  x$step <- NULL

  for (i in 1:length(x$t)) {
    if (i <= 6 || i == length(x$t)) {
      cat(
        paste0("t = ", x$t[i]),
        paste0("   "),
        paste0("(", paste0(names(x)[names(x) != "t"], collapse = ", "), ") = ", paste0("(", paste0(round(as.numeric(lapply(x[names(x) != "t"], function(xx) xx[i])), 3), collapse = ", "), ")")),
        "\n"
      )
    } else if (i == 7) {
      cat("...\n")
    }
  }
}

#' Plot a Phase Curve With Time
#'
#' Creates a plot of a phase curve
#' @param x object of class "\code{pc}" that represents a phase curve.
#' @param model_var_label a list indicating labels of model variables. Name should
#' be old variable names and values should be names of labels.
#' @param model_var_color a list indicating colors of model variable. Name should
#' be old variable names and values should be 16-bit color codes.
#' @param ... additional parameters accepted by the function "\code{geom_line}".
#'
#' @import ggplot2
#' @return a graphic object
#' @exportS3Method
#' @method plot pc
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<1", "y<1"))
#' pc <- eode_proj(x, value0 = pp(list(x = 0.2, y = 0.1)), N = 100)
#' plot(pc)
plot.pc <- function(x, model_var_label = NULL, model_var_color = NULL, ...) {
  Time <- Value <- Varible <- NULL
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

  if (!is.null(model_var_color)) {
    color_palette <- c()
    for (model_var in names(x)) {
      if (model_var %in% names(model_var_color)) {
        palette_name <- model_var
        if (!is.null(model_var_label)) {
          if (model_var %in% names(model_var_label)) {
            palette_name <- model_var_label[[which(names(model_var_label) == model_var)]]
          }
        }
        palette_value <- model_var_color[[which(names(model_var_color) == model_var)]]
        color_palette <- c(color_palette, palette_value)
        names(color_palette)[length(color_palette)] <- palette_name
      }
    }
  }
  if (!is.null(model_var_label)) {
    for (model_var in names(x)) {
      if (model_var %in% names(model_var_label)) {
        names(x)[names(x) == model_var] <- model_var_label[[which(names(model_var_label) == model_var)]]
      }
    }
  }
  df <- do.call(rbind, lapply(names(x)[names(x) != "t" & names(x) != "step"], function(name) {
    data.frame(
      Time = x$t * x$step,
      Value = x[[name]],
      Varible = name
    )
  }))
  if (is.null(model_var_color)) {
    ggplot(df, aes(x = Time, y = Value, color = Varible)) +
      geom_line(size = 0.9) +
      theme1
  } else {
    ggplot(df, aes(x = Time, y = Value, color = Varible)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = color_palette) +
      theme1
  }
}

## NON-S3------------------
#' Plot a Phase Curve in Phase Plane
#'
#' Creates a plot of a phase curve in its phase plane
#' @param x object of class "\code{pc}" that represents a phase curve.
#' @return a graphic object
#' @import ggplot2
#' @import cowplot
#' @export
#' @examples
#' eq1 <- function(x, y, r1 = 1, a11 = 2, a12 = 1) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 1, a22 = 2) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<100", "y<100"))
#' pc <- eode_proj(x, value0 = pp(list(x = 0.2, y = 0.1)), N = 50, step = 0.1)
#' pc_plot(pc)
pc_plot <- function(x) {
  if (!inherits(x, "pc")) stop("`x` should be an object of 'pc' class.")
  x$t <- NULL
  x$step <- NULL
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

  plot_li <- apply(t(combn(names(x), 2)), MARGIN = 1, function(xx) {
    ggplot() +
      geom_point(aes(x = x[[xx[1]]][1], y = x[[xx[2]]][1]), size = 3) +
      geom_line(aes(x = x[[xx[1]]], y = x[[xx[2]]]), linewidth = 0.9) +
      xlab(xx[1]) +
      ylab(xx[2]) +
      theme1
  })

  plot_grid(plotlist = plot_li)
}

#' Variable Calculator In ODE Systems
#'
#' Calculate one or several variables during the simulation of an ODE system
#' @param x the ODE system under consideration. An object of "\code{eode}" class.
#' @param formula formula that specifies how to calculate the variable to be traced.
#'
#' @import stringr
#' @return an object of "\code{pc}" class with extra variables being attached.
#' @export
#'
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<100", "y<100"))
#' pc <- eode_proj(x, value0 = pp(list(x = 0.2, y = 0.1)), N = 100)
#' pc_calculator(pc, formula = "z = x + y")
pc_calculator <- function(x, formula) {
  phase_curve <- x # change variable name
  x <- attr(phase_curve, "ode") # get original ode system

  ret_phase_curve <- phase_curve
  for (i in 1:length(formula)) {
    formula_now <- formula[i]
    for (para in names(x$parameter)) {
      if (grepl(para, formula_now)) formula_now <- gsub(para, x$parameter[[which(names(x$parameter) == para)]], formula_now)
    }
    calc_var_name <- strsplit(formula_now, " = ")[[1]][1]
    formula_now <- strsplit(formula_now, " = ")[[1]][2]
    for (model_var in x$variables) {
      if (grepl(paste0(model_var, "\\[\\-\\d+\\]"), formula_now)) {
        delayed_terms <- str_extract_all(formula_now, paste0(model_var, "\\[\\-\\d+\\]"))[[1]]
        for (delayed_term in delayed_terms) {
          delayed_time <- as.numeric(str_extract_all(delayed_term, "\\d")[[1]])
          delayed_step <- delayed_time / phase_curve$step
          value_series <- phase_curve[[which(names(phase_curve) == model_var)]]
          index_series <- phase_curve$t - delayed_step + 1
          delayed_value <- rep(0, length(index_series))
          delayed_value[index_series > 0] <- value_series[index_series > 0]
          phase_curve <- c(phase_curve, list(delayed_value))
          names(phase_curve)[length(phase_curve)] <- paste0(model_var, "_delayed_", delayed_time)
          formula_now <- gsub(paste0(model_var, "\\[\\-", delayed_time, "\\]"), paste0(model_var, "_delayed_", delayed_time), formula_now)
        }
      }
    }
    if (grepl("dt", formula_now)) {
      if (grepl("\\[", formula_now)) {
        time_interval_text <- strsplit(formula_now, " dt")[[1]][2]
        time_interval_split <- unlist(str_extract_all(time_interval_text, "[\\d+t]"))
        time_interval_start <- time_interval_split[1]
        time_interval_end <- time_interval_split[2]
        formula_now <- strsplit(formula_now, " dt")[[1]][1]
        res <- as.numeric(lapply(0:(length(phase_curve$t) - 1), function(step_now) {
          start_step <- eval(parse(text = gsub("t", "step_now", time_interval_start)))
          end_step <- eval(parse(text = gsub("t", "step_now", time_interval_end)))
          if (start_step < 0 || end_step > (length(phase_curve$t) - 1)) {
            return(NA)
          }
          sum(eval(parse(text = paste0("with(phase_curve,", formula_now, ")")))[start_step:end_step]) * phase_curve$step
        }))
        res <- c(res, rep(NA, length(phase_curve$t) - length(res)))
      } else {
        formula_now <- gsub("dt", paste0("* ", phase_curve$step), formula_now)
        res <- eval(parse(text = paste0("with(phase_curve,", formula_now, ")")))
      }
    } else {
      res <- eval(parse(text = paste0("with(phase_curve,", formula_now, ")")))
    }
    ret_phase_curve <- c(ret_phase_curve, list(res))
    names(ret_phase_curve)[length(ret_phase_curve)] <- calc_var_name
  }
  class(ret_phase_curve) <- "pc"
  ret_phase_curve
}

# TRANING------------------
#' Calculate Loss Function
#'
#' Calculate the quadratic loss function given an ODE system and observed population dynamics.
#' @param x the ODE system under consideration. An object of "\code{eode}" class.
#' @param pdat observed population dynamics. An object of "\code{pdata}" class.
#' @param step interval of time for each step. Parameter of the function "\code{eode_proj()}".
#'
#' @return a value of loss function
#' @export
#'
#' @examples
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
#' training_data <- pdata(x,
#'   init = data.frame(
#'     X_A = c(9, 19, 29, 39),
#'     Y_A = c(1, 1, 1, 1),
#'     X_C = c(5, 5, 5, 5),
#'     Y_C = c(0, 0, 0, 0)
#'   ),
#'   t = c(3, 3, 3, 3),
#'   lambda = data.frame(incidence = c(0.4, 0.8, 0.9, 0.95)),
#'   formula = "incidence = (Y_A + Y_C)/(X_A + X_C + Y_A + Y_C)"
#' )
#' eode_lossf(x, pdat = training_data)
eode_lossf <- function(x, pdat, step = 0.01) {
  loss_fun <- 0
  for (i in 1:length.pdata(pdat)) {
    if (any(grepl("dt\\[", pdat$formula))) {
      interval_max_t <- max(unlist(lapply(pdat$formula, function(formula_text) {
        time_interval_text <- strsplit(strsplit(strsplit(formula_text, "dt\\[")[[1]][2], "\\]")[[1]], "-")[[1]]
        max(c(
          eval(parse(text = paste0("with(pdat,", time_interval_text[1], ")"))),
          eval(parse(text = paste0("with(pdat,", time_interval_text[2], ")")))
        ))
      })))
      iteration_times <- round(interval_max_t / step)
    } else {
      iteration_times <- round(pdat$t[i] / step)
    }
    phase_curve <- eode_proj(x, value0 = pp(as.list(pdat$init[1, ])), N = iteration_times, step = step)
    phase_curve$step <- NULL
    for (j in 1:ncol(pdat$lambda)) {
      pred_var_name <- colnames(pdat$lambda)[j]
      pred_formula <- pdat$formula[grepl(pred_var_name, pdat$formula)]
      if (length(pred_formula) == 0) stop(paste0("no formula for predicting varible '", pred_var_name, "'"))
      if (length(pred_formula) > 1) stop(paste0("multiple formulas for predicting varible '", pred_var_name, "'. Confused."))
      observed_data <- pdat$lambda[i, j] # get observed data
      if (any(unlist(lapply(names(x$parameter), function(para) grepl(para, pred_formula))))) { # in case the formula has model parameters
        mentioned_model_paras <- names(x$parameter)[unlist(lapply(names(x$parameter), function(para) grepl(para, pred_formula)))]
        for (mentioned_model_para in mentioned_model_paras) {
          pred_formula <- gsub(mentioned_model_para, x$parameter[[which(names(x$parameter) == mentioned_model_para)]], pred_formula)
        }
      }
      if (!grepl("dt\\[", pred_formula)) { # in case we want the value of a point
        pred <- lapply(phase_curve, function(component) component[phase_curve$t == round(pdat$t[i] / step)])
        predicted_data <- eval(parse(text = paste0("with(pred,", trimws(strsplit(pred_formula, "=")[[1]][2]), ")")))
      } else { # in case we want integration over a period
        pred_time_interval_text <- gsub("\\]", "", gsub("\\[", "", strsplit(pred_formula, "dt")[[1]][length(strsplit(pred_formula, "dt")[[1]])]))
        pred_start_t <- eval(parse(text = gsub("t", "pdat$t[i]", strsplit(pred_time_interval_text, "-")[[1]][1])))
        pred_end_t <- eval(parse(text = gsub("t", "pdat$t[i]", strsplit(pred_time_interval_text, "-")[[1]][2])))
        pred <- lapply(phase_curve, function(component) {
          component[phase_curve$t > round(pred_start_t / step) &
            phase_curve$t <= round(pred_end_t / step)]
        })
        pred_formula <- trimws(strsplit(pred_formula, "dt")[[1]][1])
        predicted_data <- sum(eval(parse(text = paste0("with(pred,", trimws(strsplit(pred_formula, "=")[[1]][2]), ")")))) * step
      }
      loss_fun <- loss_fun + (observed_data - predicted_data)^2
    }
  }
  loss_fun
}

#' Grid Search For Optimal Parameters
#'
#' Find optimal parameters in the ODE system using grid search method.
#'
#' @param x the ODE system under consideration. An object of "\code{eode}" class.
#' @param pdat observed population dynamics. An object of "\code{pdata}" class.
#' @param step interval of time for running simulations. Parameter of the function "\code{eode_proj()}".
#' @param space space
#'
#' @return a data frame showing attempted parameters along with the corresponding values of loss function.
#' @export
#'
#' @examples
#' \donttest{
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
#' training_data <- pdata(x,
#'   init = data.frame(
#'     X_A = c(9, 19, 29, 39),
#'     Y_A = c(1, 1, 1, 1),
#'     X_C = c(5, 5, 5, 5),
#'     Y_C = c(0, 0, 0, 0)
#'   ),
#'   t = c(3, 3, 3, 3),
#'   lambda = data.frame(incidence = c(0.4, 0.8, 0.9, 0.95)),
#'   formula = "incidence = (Y_A + Y_C)/(X_A + X_C + Y_A + Y_C)"
#' )
#' res <- eode_gridSearch(x, pdat = training_data, space = list(beta = seq(0.05, 0.15, 0.01)))
#' res
#' }
eode_gridSearch <- function(x, pdat, space, step = 0.01) {
  grids <- expand.grid(space)
  lossf <- data.frame(grids, lossf = 0)
  message("Start Grid Search...\n")
  for (i in 1:nrow(grids)) {
    paras <- data.frame(grids[i, ])
    colnames(paras) <- colnames(grids)
    try_ode <- eode_set_parameter(x, ParaList = as.list(paras))
    try_lossf <- eode_lossf(try_ode, pdat, step = step)
    message(paste0("Parameter: ", paste(colnames(paras), paras, sep = " = ", collapse = ", "), ", Loss Function: ", round(try_lossf, 3), "\n"))
    lossf[i, "lossf"] <- try_lossf
  }
  message("Finished.\n")

  optPara <- data.frame(lossf[which.min(lossf$lossf)[1], -ncol(lossf)])
  colnames(optPara) <- colnames(grids)
  message("Optimal Parameters:", paste(colnames(optPara), optPara, sep = " = ", collapse = ", "), "\n") #
  message("Loss Function:", min(lossf$lossf), "\n")
  lossf
}

#' Simulated Annealing For Optimal Parameters
#'
#' Find optimal parameters in the ODE system using simulated annealing.
#' @param x the ODE system under consideration. An object of "\code{eode}" class.
#' @param pdat observed population dynamics. An object of "\code{pdata}" class.
#' @param paras parameters to be optimised. A character vector. If multiple parameters
#' are specified, the simulation annealing process will proceed by altering multiple
#' parameters at the same time, and accept an alteration if it achieves a lower value of
#' the loss function. Default is "ALL", which means to choose all the parameters.
#' @param max_disturb maximum disturbance in proportion. The biggest disturbance acts on
#' parameters at the beginning of the simulated annealing process.
#' @param AnnN steps of simulated annealing.
#' @param step interval of time for running simulations. Parameter of the function "\code{eode_proj()}".
#' @param prop.train proportion of training data set. In each step of annealing, a proportion
#' of dataset will be randomly decided for model training, and the rest for model validation.
#'
#' @return a data frame showing attempted parameters along with the corresponding values of loss function.
#' @export
#'
#' @examples
#' \donttest{
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
#' training_data <- pdata(x,
#'   init = data.frame(
#'     X_A = c(9, 19, 29, 39),
#'     Y_A = c(1, 1, 1, 1),
#'     X_C = c(5, 5, 5, 5),
#'     Y_C = c(0, 0, 0, 0)
#'   ),
#'   t = c(3, 3, 3, 3),
#'   lambda = data.frame(incidence = c(0.4, 0.8, 0.9, 0.95)),
#'   formula = "incidence = (Y_A + Y_C)/(X_A + X_C + Y_A + Y_C)"
#' )
#' res <- eode_simuAnnealing(x, pdat = training_data, paras = "beta", max_disturb = 0.05, AnnN = 20)
#' res
#' }
eode_simuAnnealing <- function(x, pdat, paras = "ALL", max_disturb = 0.2, AnnN = 100, step = 0.01, prop.train = 1) {
  whole_pdat <- pdat
  message("Start Simulated Annealing...\n")
  if (prop.train != 1) {
    sample_index <- sample(1:length.pdata(whole_pdat), length.pdata(whole_pdat))
    sample_selected <- sample_index > length.pdata(whole_pdat) * prop.train
    pdat <- whole_pdat[sample_selected]
    pvali <- whole_pdat[!sample_selected]
  }
  lossf_now <- eode_lossf(x, pdat, step = step)
  if (prop.train != 1) {
    if (length(pvali) != 0) {
      lossf_validation <- eode_lossf(x, pvali, step = step)
      sa_recording <- data.frame(x$parameter, lossf = lossf_now, lossf_validation = lossf_validation)
    } else {
      warning("Data set too small. No validation data available.")
      sa_recording <- data.frame(x$parameter, lossf = lossf_now, lossf_validation = NA)
    }
  } else {
    sa_recording <- data.frame(x$parameter, lossf = lossf_now)
  }
  message("Step 0 - loss function =", lossf_now, "\n")
  for (i in 1:AnnN) {
    if (prop.train != 1) {
      sample_index <- sample(1:length.pdata(whole_pdat), length.pdata(whole_pdat))
      sample_selected <- sample_index > length.pdata(whole_pdat) * prop.train
      pdat <- whole_pdat[sample_selected]
      pvali <- whole_pdat[!sample_selected]
    }
    disturb_i <- max_disturb - max_disturb / AnnN * (i - 1)
    message("Disturbance:", paste0(disturb_i * 100, "%"), "...\n")
    if (paras[1] == "ALL") {
      parameter_sets <- x$parameter
    } else {
      parameter_sets <- x$parameter[which(names(x$parameter) %in% paras)]
    }
    paras_new <- lapply(parameter_sets, function(val) eval(parse(text = deparse(val))) * runif(1, min = 1 - disturb_i, max = 1 + disturb_i))
    x_new <- eode_set_parameter(x, ParaList = paras_new)
    lossf_new <- eode_lossf(x_new, pdat, step = step)
    if (lossf_new < lossf_now) {
      x <- x_new
      lossf_now <- lossf_new
      message("Step", i, "- parameter resetted:", paste(names(parameter_sets), round(as.numeric(parameter_sets), 3), sep = "=", collapse = ", "), "- loss function =", lossf_now, "\n")
    } else {
      message("Step", i, "- parameters remain unchanged.\n")
    }

    if (prop.train != 1) {
      if (length(pvali) != 0) {
        lossf_validation <- eode_lossf(x, pvali, step = step)
        sa_recording <- rbind(sa_recording, data.frame(x$parameter, lossf = lossf_now, lossf_validation = lossf_validation))
      } else {
        sa_recording <- rbind(sa_recording, data.frame(x$parameter, lossf = lossf_now, lossf_validation = NA))
      }
    } else {
      sa_recording <- rbind(sa_recording, data.frame(x$parameter, lossf = lossf_now))
    }
  }
  message("Finished.\n")
  message("Parameters:", paste(names(x$parameter), round(as.numeric(x$parameter), 3), sep = "=", collapse = ", "), "- loss function =", lossf_now, "\n")

  sa_recording
}

# SENSITIVITY-------------------
#' Sensitivity Analysis
#'
#' Run a sensitivity analysis on an ODE system.
#' @param x object of class "\code{pc}" that represents the ODE system under
#' consideration.
#' @param valueSpace a list indicating initial conditions and parameters. Model
#' variables must be included to specify initial values of each variable. Values
#' can be a vector, indicating all the potential values to be considered in the
#' sensitivity analysis.
#' @param N number of iterations
#' @param step interval of time for running simulations. Parameter of the function "\code{eode_proj()}".
#'
#' @return an object of "\code{pcfamily}" class, each component having three
#' sub-components:
#' \code{$grid_var}: variables or parameters whose values are changed throughout
#' the sensitivity analysis.
#' \code{$fixed_var}: variables whose values are not changed.
#' \code{$pc}: phase curve. An object of "\code{pc}" class.
#' @export
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<100", "y<100"))
#' eode_sensitivity_proj(x, valueSpace = list(x = c(0.2, 0.3, 0.4), y = 0.1, a11 = 1:3), N = 100)
eode_sensitivity_proj <- function(x, valueSpace, N, step = 0.01) {
  if (!all(x$variables %in% names(valueSpace))) {
    stop(paste0("variable '", x$variables[!(x$variables %in% names(valueSpace))], "' not specified by `valueSpace`."))
  }
  grid_var_name <- names(valueSpace)[which(as.logical(lapply(valueSpace, function(xx) length(xx) > 1)))]
  grids <- expand.grid(valueSpace)
  res <- list()
  for (i in 1:nrow(grids)) {
    message(paste0("Parameter Settings: ", paste(grid_var_name, grids[i, grid_var_name], sep = " = ", collapse = ", "), "\n"))
    inits <- as.list(grids[i, colnames(grids) %in% x$variables])
    if (length(inits) == 1) names(inits) <- colnames(grids)[colnames(grids) %in% x$variables]
    paras <- as.list(grids[i, !(colnames(grids) %in% x$variables)])
    if (length(paras) == 1) names(paras) <- colnames(grids)[!(colnames(grids) %in% x$variables)]
    res <- c(res, list(list(
      grid_var = unlist(c(inits, paras)[grid_var_name]),
      fixed_var = names(grids)[!(names(grids) %in% grid_var_name)],
      pc = eode_proj(eode_set_parameter(x, paras), value0 = pp(inits), N, step = step)
    )))
  }
  class(res) <- "pcfamily"
  res
}

#' Print Brief Details of a Phase Curve Family
#'
#' Prints a very brief description of a phase curve family.
#' @param x Object of class "\code{pcfamily}".
#' @param ... Ignored.
#' @return No value
#' @exportS3Method
#' @method print pcfamily
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<100", "y<100"))
#' eode_sensitivity_proj(x, valueSpace = list(x = c(0.2, 0.3, 0.4), y = 0.1, a11 = 1:3), N = 100)
print.pcfamily <- function(x, ...) {
  pcf <- x # change variable name
  cat("Phase Curve Family:\n")
  texts <- as.data.frame(do.call(rbind, lapply(pcf, function(xx) {
    c(
      paras = paste(names(xx$grid_var), xx$grid_var, sep = " = ", collapse = ", "),
      time = paste0("time = 0 ~ ", xx$pc$t[length(xx$pc$t)] * xx$pc$step)
    )
  })))
  texts$paras <- do.call(c, lapply(texts$paras, function(xx) {
    paste0(xx, paste(rep(" ", nchar(xx) - max(nchar(texts$paras)) + 4), collapse = ""))
  }))
  for (i in 1:nrow(texts)) {
    cat(paste0("Phase Curve ", i, ": "))
    cat(texts$paras[i], texts$time[i], "\n")
  }
}

#' Plot Phase Curve Family
#'
#' Creates a plot of a phase curve family.
#'
#' @param x object of class "\code{pcfamily}" that represents a phase curve family.
#' @param model_var_label a list indicating labels of model variables. Name should
#' be old variable names and values should be names of labels.
#' @param model_var_color a list indicating colors of model variable. Name should
#' be old variable names and values should be 16-bit color codes.
#' @param facet_grid_label a list indicating facet labels. Name should be old
#' variable names and values should be facet labels.
#' @param facet_paras a logical vector indicating whether facet?
#' @param ... other parameters
#'
#' @import ggplot2
#' @return a graphic object
#' @exportS3Method
#' @method plot pcfamily
#' @examples
#' eq1 <- function(x, y, r1 = 4, a11 = 1, a12 = 2) (r1 - a11 * x - a12 * y) * x
#' eq2 <- function(x, y, r2 = 1, a21 = 2, a22 = 1) (r2 - a21 * x - a22 * y) * y
#' x <- eode(dxdt = eq1, dydt = eq2, constraint = c("x<1", "y<1"))
#' value_space <- list(x = c(0.2, 0.3, 0.4), y = 0.1, a11 = 1:3)
#' plot(eode_sensitivity_proj(x, valueSpace = value_space, N = 100))
plot.pcfamily <- function(x, model_var_label = NULL, model_var_color = NULL,
                          facet_grid_label = NULL, facet_paras = TRUE, ...) {
  Time <- Value <- Varible <- NULL
  if (!is.null(model_var_color)) {
    color_palette <- c()
    for (model_var in names(x[[1]]$pc)) {
      if (model_var %in% names(model_var_color)) {
        palette_name <- model_var
        if (!is.null(model_var_label)) {
          if (model_var %in% names(model_var_label)) {
            palette_name <- model_var_label[[which(names(model_var_label) == model_var)]]
          }
        }
        palette_value <- model_var_color[[which(names(model_var_color) == model_var)]]
        color_palette <- c(color_palette, palette_value)
        names(color_palette)[length(color_palette)] <- palette_name
      }
    }
  }
  if (!is.null(model_var_label)) {
    for (model_var in names(x[[1]]$pc)) {
      if (model_var %in% names(model_var_label)) {
        for (i in 1:length(x)) {
          names(x[[i]]$pc)[names(x[[i]]$pc) == model_var] <- model_var_label[[which(names(model_var_label) == model_var)]]
        }
      }
    }
  }
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
      strip.text = element_text(size = 14),
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

  df <- do.call(rbind, lapply(x, function(pc_component) {
    data.frame(do.call(rbind, lapply(names(pc_component$pc)[names(pc_component$pc) != "t" & names(pc_component$pc) != "step"], function(name) {
      data.frame(
        Time = pc_component$pc$t * pc_component$pc$step,
        Value = pc_component$pc[[name]],
        Varible = name
      )
    })), as.list(pc_component$grid_var))
  }))

  if (facet_paras) {
    if (!is.null(facet_grid_label)) {
      for (grid_var_name in names(x[[1]]$grid_var)) {
        if (grid_var_name %in% names(facet_grid_label)) {
          df[, grid_var_name] <- factor(paste(
            facet_grid_label[[which(names(facet_grid_label) == grid_var_name)]],
            df[, grid_var_name],
            sep = " = "
          ), levels = paste(facet_grid_label[[which(names(facet_grid_label) == grid_var_name)]], unique(df[, grid_var_name]), sep = " = "))
        }
      }
    } else {
      for (grid_var_name in names(x[[1]]$grid_var)) {
        df[, grid_var_name] <- factor(paste(
          grid_var_name, df[, grid_var_name],
          sep = " = "
        ), levels = paste(grid_var_name, unique(df[, grid_var_name]), sep = " = "))
      }
    }

    g <- ggplot(df, aes(x = Time, y = Value, color = Varible)) +
      geom_line(size = 0.9) +
      theme1
    if (!is.null(model_var_color)) g <- g + scale_color_manual(values = color_palette)

    if (length(x[[1]]$grid_var) == 2) {
      g + eval(parse(text = paste0("facet_grid(", paste(names(x[[1]]$grid_var), collapse = " ~ "), ")")))
    } else if (length(x[[1]]$grid_var) == 1) {
      g + eval(parse(text = paste0("facet_wrap(~ ", names(x[[1]]$grid_var), ")")))
    } else {
      g + eval(parse(text = paste0("facet_wrap(~ ", paste(names(x[[1]]$grid_var), collapse = " + "), ")")))
    }
  } else {
    if (length(x[[1]]$grid_var) > 1) {
      stop("Cannot create a plot of multiple grid varibles facet by model variables. Please use `facet_paras = TRUE`.")
    }
    focus_var_name <- names(x[[1]]$grid_var)

    alpha_gradients <- eval(parse(text = paste0("unique(df$", focus_var_name, ")")))
    alpha_index_steps <- floor(length(alpha_gradients) / 4)
    alpha_breaks <- alpha_gradients[c(
      1, 1 + alpha_index_steps,
      1 + 2 * alpha_index_steps,
      1 + 3 * alpha_index_steps,
      1 + 4 * alpha_index_steps
    )]


    g <- eval(parse(text = paste0("ggplot(df, aes(x = Time, y = Value, color = Varible, group = ", focus_var_name, ", alpha = ", focus_var_name, "))"))) +
      geom_line(size = 0.9) + facet_wrap(~ df$Varible, nrow = 1) + guides(color = "none") +
      scale_alpha_continuous(breaks = alpha_breaks) + theme1
    if (!is.null(model_var_color)) g <- g + scale_color_manual(values = color_palette)
    g
  }
}
