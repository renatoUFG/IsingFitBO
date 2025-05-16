
#' Summary Method for IsingFitBO Objects
#'
#' @description
#' Provides comprehensive summary statistics for Ising model estimation results.
#'
#' @param object An object of class `IsingFit`
#' @param digits Number of significant digits (default=3)
#' @param ... Additional arguments
#'
#' @return
#' An object of class `summary.IsingFitBO` containing network properties and optimization details.
#'
#' @export
summary.IsingFitBO <-
function(object, ...)
{
  cat("\tNetwork Density:\t\t", round(mean(object$weiadj[upper.tri(object$weiadj)]!=0),2),"\n",
      "Gamma:\t\t\t",round(object$gamma,2),"\n",
      "Rule used:\t\t",ifelse(object$AND,"And-rule","Or-rule"),"\n",
      "Analysis took:\t\t",format(object$time,format="%s"),"\n"
      )
}
