\name{IsingFit}
\alias{IsingFit}
\title{
Network estimation using the eLasso method with Bayesian Optimization
}
\description{
An extension of the IsingFit package (van Borkulo et al., 2014) that implements Bayesian Optimization for lambda hyperparameter tuning. This adaptation retains the original network estimation logic while adding automated hyperparameter selection. Licensed under GPL-2 as a derivative work. 
}
\usage{
IsingFitBO(x, method="BayesOpt",family = "binomial", AND = TRUE,
          niter=20, plot = TRUE, gamma_hyp = 0.25, ...)
}

\arguments{
  \item{x}{
Input matrix. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}

\item{method}{
The default method is 'BayesOpt', which uses Bayesian Optimization for lambda hyperparameter tuning. If the user sets the method argument to 'Grid', the function will instead use a predefined grid of nine values for hyperparameter tuning.
}

\item{family}{
The default is 'binomial', treating the data as binary. Currently, this procedure is only supported for binary data.
}
  \item{AND}{
Logical. Can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network. Defaults to TRUE.
}

\item{niter}{
The number of iterations for the Bayesian Optimization procedure. Default is 20.
}

\item{gamma}{
A value of hyperparameter gamma in the extended BIC. Can be anything between 0 and 1. Defaults to .25.
}
  \item{plot}{
Logical. Should the resulting network be plotted?
}

\item{\dots}{
Arguments sent to \code{qgraph}.
}

}

\value{
IsingFit returns (invisibly) a 'IsingFit' object that contains the following items:
\item{weiadj }{The weighted adjacency matrix.}
\item{thresholds }{Thresholds of the variables.}
\item{q }{The object that is returned by qgraph (class 'qgraph').}
\item{gamma }{The value of hyperparameter gamma.}
\item{AND }{A logical indicating whether the AND-rule is used or not. If not, the OR-rule is used.}
\item{time }{The time it took to estimate the network.}
\item{asymm.weights }{The (asymmetrical) weighted adjacency matrix before applying the AND/OR rule.}
\item{lambda.values }{The values of the tuning parameter per node that ensured the best fitting set of neighbors.}
}

\references{
Chen, J., & Chen, Z. (2008). Extended bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Foygel, R., & Drton, M. (2011). Bayesian model choice and information criteria in sparse generalized linear models. arXiv preprint arXiv:1112.5635.

Ravikumar, P., Wainwright, M. J., & Lafferty, J. D. (2010). High-dimensional Ising model selection using l1-regularized logistic regression. The Annals of Statistics, 38, 1287 - 1319.

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific Reports 4, 5918; DOI:10.1038/srep05918. 
}
\author{
Claudia D. van Borkulo, Sacha Epskamp; with contributions from Alexander Robitzsch and Mihai Alexandru Constantin

Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
}
\note{
This function extends the original \code{IsingFit} package (van Borkulo et al., 2014). 
The Bayesian Optimization feature was added by Renato Rodrigues Silva (2025).
}


\examples{
library("IsingSampler")

### Simulate dataset ###
# Input:
N <- 6 # Number of nodes
nSample <- 1000 # Number of samples

# Ising parameters:
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.8, 0.2)),N,N) * runif(N^2,0.5,2)
Graph <- pmax(Graph,t(Graph))
diag(Graph) <- 0
Thresh <- -rowSums(Graph) / 2

# Simulate:
Data <- IsingSampler(nSample, Graph, Thresh)

### Fit using IsingFitBO ###
Res <- IsingFitBO(Data, method="BayesOpt",
family='binomial', AND = TRUE, niter=20, plot = TRUE, 
gamma_hyp = 0.25)

# Plot results:
library("qgraph")
layout(t(1:2))
qgraph(Res$weiadj,fade = FALSE)
title("Estimated network")
qgraph(Graph,fade = FALSE)
title("Original network")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
