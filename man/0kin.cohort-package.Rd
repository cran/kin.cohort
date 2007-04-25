\name{kin.cohort-package}
\alias{kin.cohort-package}
\docType{package}
\title{
Analysis of kin-cohort studies
}
\description{
Analysis of kin-cohort studies. kin.cohort provides estimates of age-specific 
cumulative risk of a disease for carriers (penetrance) and noncarriers of a mutation. 
The cohorts are  retrospectively built from relatives of probands for whom the genotype is known. 
Currently the method of moments and marginal maximum likelihood are implemented. Confidence intervals 
are calculated from bootstrap samples.
}
\details{
\tabular{ll}{
Package: \tab kin.cohort\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2007-04-25\cr
License: \tab GPL version 2\cr
}
The main function, kin.cohort, will receive as arguments data vectors, the method of estimation desired,
currently only two methods are implemented: the method of moments based on 
Kaplan-Meier estimates (Wacholder´s et al., 1989) and the Marginal Maximum Likelihood method (Chatterjee and Wacholder, 2001) . 
If confidence intervals are desired, the number of bootstrap samples can be specified. 
Print and plot methods are provided for inspecting the output.
}
\author{
Victor Moreno, Nilanjan Chatterjee, Bhramar Mukherjee.

Maintainer: Victor Moreno <v.moreno@iconcologia.net>

}
\references{
Chatterjee N and Wacholder S. A Marginal Likelihood Approach for Estimating Penetrance from Kin-Cohort Designs. Biometrics. 2001; 57: 245-52.

Wacholder S, Hartge P, Struewing JP, Pee D, McAdams M, Lawrence B, Tucker MA. The kin-cohort study for estimating penetrance. American Journal of Epidemiology. 1998; 148: 623-9.
}
\keyword{ package }

}
