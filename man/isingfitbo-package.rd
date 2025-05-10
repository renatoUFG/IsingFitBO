\name{IsingFitBO-package}
\alias{IsingFitBO-package}
\docType{package}
\title{
Network estimation using the eLasso method with an Bayesian Optimization for lambda parameters
}
\description{
An extension of the IsingFit package (van Borkulo et al., 2014) that implements Bayesian Optimization for lambda hyperparameter tuning. This adaptation retains the original network estimation logic while adding automated hyperparameter selection. Licensed under GPL-2 as a derivative work. 
}
\details{
\tabular{ll}{
Package: \tab IsingFitBO\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2025-05-12\cr

License: \tab GPL-2\cr 
}
}
\author{
Renato Rodrigues Silva, Claudia D. van Borkulo, Sacha Epskamp;
with contributions from Alexander Robitzsch and Mihai Alexandru Constantin

Maintainer: Renato Rodrigues Silva <renato.rrsilva@ufg.br>
}
\references{
van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific Reports 4, 5918; DOI:10.1038/srep05918. 

Wang, J. (2023). An intuitive tutorial to Gaussian process regression. IEEE. https://doi.org/10.1109/XXX.0000.0000000

}
% ~~ Optionally other standard keywords, one per line, from file ~~
% ~~ KEYWORDS in the R documentation directory ~~

\note{
This package extends the original \code{IsingFit} (van Borkulo et al., 2014). 
The Bayesian Optimization feature was added by Renato Rodrigues Silva (2025).
}