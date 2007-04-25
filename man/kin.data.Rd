\name{kin.data}
\alias{kin.data}
\docType{data}
\title{ sample data for kin-cohort analysis}
\description{
  Simulated data of a study on the penetrance of breast cancer for carriers 2 mutations.
}
\usage{data(kin.data)}
\format{
  A data frame with 15341 observations on the following 5 variables.
  \describe{
    \item{\code{age}}{age at diagnosis or at last follow-up}
    \item{\code{cancer}}{disease status (1: breast cancer, 0: no breast cancer}
    \item{\code{gen1}}{gen1 genotypes of proband}
    \item{\code{gen2}}{gen2 genotypes of proband}
	 \item{\code{rel}}{relationship with proband 1:parent or offspring, 2:sibling}
    \item{\code{family}}{family id}
  }
}
\examples{
data(kin.data)
}
\keyword{datasets}
