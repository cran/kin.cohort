\name{kc.marginal}
\alias{kc.marginal}
\title{ Marginal Maximum Likelihood estimation of kin-cohort data }
\description{
This function estimates cumulative risk and hazard at given ages for carriers and noncarriers of a mutation based on the probands genotypes. It uses the Marginal Maximum Likelihood estimation method (Chatterjee and Wacholder, 2001). Piece-wise exponential distribution is assumed for the survival function.
}
\usage{
kc.marginal(t, delta, genes, r, knots, f, pw = rep(1,length(t)), 
            set = NULL, B = 1, maxit = 1000, tol = 1e-5, subset,
            logrank=TRUE, trace=FALSE)
}
\arguments{
  \item{t}{ time variable. Usually age at diagnosis or at last follow-up }
  \item{delta}{ disease status (1: event, 0: no event }
  \item{genes}{ factor or numeric vector (1 gene), matrix or dataframe (2 genes) with genotypes of proband numeric. factors and data.frame with factors are prefered in order to use user-defined labels. Otherwise use codes (1:noncarrier, 2: carrier, 3: homozygous carrier) }
  \item{r}{ relationship with proband 1:parent, 2:sibling 3:offspring 0:proband. Probands will be excluded from analysis and offspring will be recoded 1 internally.}
  \item{knots}{ time points (ages) for cumulative risk and hazard estimates }
  \item{f}{ vector of mutation allele frequencies in the population }
  \item{pw}{ prior weights, if needed }
  \item{set}{ family id (only needed for bootstrap)   }
  \item{B}{ number of boostrap samples (only needed for bootstrap) }
  \item{maxit}{ max number of iterations for the EM algorithm }
  \item{tol}{ convergence tolerance }
  \item{subset}{ logical condition to subset data }
  \item{logrank}{ Perform a logrank test }
  \item{trace}{ Show iterations for bootstrap }
}

\value{
  object of classes "kin.cohort" and "chatterjee". 
  \item{cumrisk }{matrix with cumulative risk estimates for noncarriers, carriers and the cumulative risk ratio. Estimates are given for the times indicated in the knot vector}
  \item{hazard }{matrix with hazard estimates for noncarriers, carriers and the hazard ratio. Estimates are given for the times indicated in the knot vector}
  \item{knots}{vector of knots}
  \item{conv}{if the EM algorithm converged}
  \item{niter}{number of iterations needed for convergence}
  \item{ngeno.rel}{number of combinations of genotypes in the relatives}
  \item{events}{matrix with number of events and person years per each knot}
  \item{logHR}{mean log hazard ratio estimate (unweighted)}
  \item{logrank}{logrank test. If 2 genes, for the main effects, the cross-classification and the stratified tests}
  \item{call}{copy of call}

  if bootstrap confidence intervals are requested (B>1) then the returned object is of classes "kin.cohort.boot" and "chatterjee" with previous items packed in value estimate and each bootstrap sample packed in matrices. 
}
\references{ 
Chatterjee N and Wacholder S. A Marginal Likelihood Approach for Estimating Penetrance from Kin-Cohort Designs. Biometrics. 2001; 57: 245-52. 
}
\note{ 
This function is best called by kin.cohort than directly
}
\seealso{ 
\code{\link{kin.cohort}}, \code{\link{print.kin.cohort}}, \code{\link{plot.kin.cohort}} 
}
\examples{
\dontrun{
data(kin.data)
attach(kin.data)
res.mml<- kc.marginal(age, cancer, gen1, rel, knots=c(30,40,50,60,70,80), f=0.02)
res.mml
}
}
\keyword{survival}