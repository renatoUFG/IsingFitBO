#' Plot Method for IsingFit Objects
#'
#' @description
#' Visualizes the results of an Ising model estimation with Bayesian optimization.
#'
#' @param x An object of class `IsingFit`
#' @param type Plot type: "network" (default), "convergence", or "parameters"
#' @param main Plot title (optional)
#' @param ... Additional arguments passed to plot functions
#'
#' @return
#' Invisibly returns the plot object. For "network" type, returns a `qgraph` object.
#'
#' @seealso
#' \code{\link{qgraph}} for network visualization options
#'
#' @export
plot.IsingFitBO <-
function(x,...) qgraph(x$q,DoNotPlot = FALSE, ...)
