\name{IsingBO-methods}
\alias{print.IsingFitBO}
\alias{plot.IsingFitBO}
\alias{summary.IsingFitBO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for IsingFitBO objects
}
\description{
Print method prints the IsingFitBO output , plot method plots the estimated network (with the \code{qgraph} package), and summary method returns density of the network, the value of gamma used, the rule used, and the time the analysis took.
}
\usage{
\method{print}{IsingFitBO}(x, \dots)
\method{summary}{IsingFitBO}(object, \dots)
\method{plot}{IsingFitBO}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
output of \code{\link{IsingFitBO}}
}
  \item{object}{
output of \code{\link{IsingFitBO}}
}
  \item{\dots}{
  Arguments sent to qgraph. Only used in plot method.
}
}
\author{
Renato Rodrigues Silva and Claudia van Borkulo
}