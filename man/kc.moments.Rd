\name{kc.moments}
\alias{kc.moments}

\title{ Kin-cohort estimation of penetrance by the method of moments }
\description{
This function estimates cumulative risk and hazard at given ages for carriers and noncarriers of a mutation based on the probands genotypes. It uses the method of moments described by Wacholder et al (1998)
}
\usage{
kc.moments(t, delta, genes, r, knots, f, pw = rep(1,length(t)), 
           set = NULL, B = 1, logrank = TRUE, subset, trace=FALSE)
}

\arguments{
  \item{t}{ time variable. Usually age at diagnosis or at last follow-up }
  \item{delta}{ disease status (1: event, 0: no event }
  \item{genes}{ genotype of proband numeric. A factor is preferred, otherwise numeric code of genotypes (1: noncarrier, 2:carrier, [3: homozygous carrier]) }
  \item{r}{ relationship with proband 1:parent, 2:sibling 3:offspring 0:proband. Probands will be excluded from analysis and offspring will be recoded 1 internally. }
  \item{knots}{ time points (ages) for cumulative risk and hazard estimates }
  \item{f}{ mutation allele frequency in the population }
  \item{pw}{ prior weights, if needed }
  \item{set}{ family id (only needed for bootstrap) }
  \item{B}{ number of boostrap samples (only needed for bootstrap) }
  \item{logrank}{ if logrank test is desired }
  \item{subset}{ logical condition to subset data }
  \item{trace}{ Show iterations for bootstrap }
}

\value{
  object of classes "kin.cohort" and "wacholder". 
  \item{cumrisk}{matrix of dimension (number of knots x 3) with cumulative risk festimates or noncarriers, carriers and the cumulative risk ratio}
  \item{knots}{vector of knots}
  \item{km}{object class survfit (package survival) }
  \item{logrank}{p-value of the logrank test}  
  \item{events}{matrix with number of events and person years per each knot}
  \item{call}{copy of call}

  if bootstrap confidence intervals are requested (B>1) then the returned object is of classes "kin.cohort.boot" and "wacholder" with previous items packed in value estimate and each bootstrap sample packed in matrices. 
}
\references{ 
Wacholder S, Hartge P, Struewing JP, Pee D, McAdams M, Lawrence B, Tucker MA. 
The kin-cohort study for estimating penetrance. American Journal of Epidemiology. 1998; 148: 623-9.
}
\note{ 
This function is best called by kin.cohort than directly
}
\seealso{ 
\code{\link{kin.cohort}}, \code{\link{print.kin.cohort}}, \code{\link{plot.kin.cohort}} 
}
\examples{
data(kin.data)
attach(kin.data)
res.km<- kc.moments(age, cancer, gen1, rel, knots=c(30,40,50,60,70,80), f=0.02)
res.km
}
\keyword{survival}