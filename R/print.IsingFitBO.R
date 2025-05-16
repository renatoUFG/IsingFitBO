#' Print Method for IsingFit Objects
#'
#' @description
#' Displays concise output of Ising model estimation results.
#'
#' @param x An object of class `IsingFit`
#' @param ... Additional arguments
#'
#' @return
#' Invisibly returns the input object.
#'
#' @export
print.IsingFitBO <-
function(x, ...)
{
  cat("Estimated network:\n")

  print(round(x$weiadj,2))

  cat("\n\nEstimated Thresholds:\n")

  print(x$thresholds)
}
