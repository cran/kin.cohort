\name{internal}
\alias{pyear}
\alias{pwexp}
\alias{cp1exact}
\alias{cp32exact}
\alias{tpgg00exact}
\alias{mendelian}
\alias{mendelian.combine}

\title{ Internal functions for marginal method }
\description{
  \item{pyear}{ calculates number of events and person years}
   \item{pwexp}{ estimates survival and hazard for piece-wise exponential model}
   \item{cp1exact}{ calculates the mendelian probabilities of carrying the mutation conditional on the proband genotype. Case of 2 genotypes (noncarriere/carrier).}
   \item{cp32exact}{ calculates the mendelian probabilities of carrying the mutation conditional on the proband genotype. Case of 3 genotypes.}
   \item{cp1exact}{ calculates the mendelian probabilities of carrying the mutation conditional on the proband genotype for the combinations of 2 genes.}
   \item{mendelian}{ calculates the mendelian probabilities of carrying the mutation conditional on the proband genotype for 1 gene.}
   \item{mendelian.combine}{ combines mendelian probabilities of carrying the mutation conditional on the proband genotype for 2 genes.}
}

\seealso{ \code{\link{kc.marginal}} }
\keyword{internal}