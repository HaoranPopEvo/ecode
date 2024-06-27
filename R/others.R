# Load the value of the option on package startup
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\033[38;5;246mIf you use the ReporterScore package in published research, please cite:\n\nHaoran Wu, ecode: An R package to investigate community dynamics in ordinary differential equation systems, Ecological Modelling, 2024, https://doi.org/10.1016/j.ecolmodel.2024.110676 \033[39m")
}
