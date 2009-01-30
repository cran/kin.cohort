\name{kin.cohort}
\alias{kin.cohort}

\title{ Analysis of kin-cohort data }
\description{
This function estimates cumulative risk at given ages for carriers and noncarriers
of a mutation based on the probands genotypes. It can use the Marginal Maximum
Likelihood estimation method (Chatterjee and Wacholder, 2001) or the method of
moments (Wacholder et al, 2001). Bootstrap confidence intervals can be requested.
}
\usage{
kin.cohort(\dots, method = c("marginal", "mml", "chatterjee", 
                             "moments",  "km",  "watcholder"))
}
\arguments{
  \item{\dots}{ see \code{\link{kc.marginal}} and \code{\link{kc.moments} for details} }
  \item{method}{ choose estimation method: Marginal Maximum Likelihood
  (selected by "marginal", "mml", "chatterjee") or method of moments
  (selected by "moments", "km", "watcholder")}
}
\details{
This function is a wrapper that will call \code{\link{kc.marginal}} or \code{\link{kc.moments}} depending on the argument method.
}

\references{
Wacholder S, Hartge P, Struewing JP, Pee D, McAdams M, Lawrence B, Tucker MA. The kin-cohort study for estimating penetrance. American Journal of Epidemiology. 1998; 148: 623-9.

Chatterjee N and Wacholder S. A Marginal Likelihood Approach for Estimating Penetrance from Kin-Cohort Designs. Biometrics. 2001; 57: 245-52.
}
\author{ Victor Moreno, Nilanjan Chatterjee, Bhramar Mukherjee. 

Maintainer: Victor Moreno <v.moreno@iconcologia.net>
}
\seealso{ 
\code{\link{kc.marginal}}, \code{\link{kc.moments}} 
}
\examples{
data(kin.data)
attach(kin.data)

res.k<-   kin.cohort(age, cancer, gen1, rel, knots=c(30,40,50,60,70,80), f=0.02, 
                     method="km")
res.k          
plot(res.k)
plot(res.k,what="crr")

set.seed(1)
res.k.b<- kin.cohort(age, cancer, gen1, rel, knots=c(30,40,50,60,70,80), f=0.02, 
                     set=family, method="km", B=10)
res.k.b
plot(res.k.b)
plot(res.k.b,what="crr")

res.m<-   kin.cohort(age, cancer, gen1, rel, knots=c(30,40,50,60,70,80), f=0.02, 
                    method="mml")
res.m
plot(res.m)
plot(res.m, what="hr")

res.m2<-  kin.cohort(age, cancer, data.frame(gen1,gen2), rel, 
                     knots=c(30,40,50,60,70,80), f=c(0.02,0.01), method="mml")
res.m2
plot(res.m2)
plot(res.m2, what="hr")

set.seed(1)
res.m.b<- kin.cohort(age, cancer, gen1, rel, knots=c(30,40,50,60,70,80), f=0.02, 
                     set=family, method="mml", B=10)
res.m.b
plot(res.m.b)
plot(res.m.b, what="hr")
}
\keyword{survival}