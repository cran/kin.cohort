\name{methods}
\alias{print.kin.cohort}
\alias{print.kin.cohort.boot}
\alias{plot.kin.cohort}
\alias{plot.kin.cohort.boot}

\title{ methods for print and plot }
\description{
Functions to print a formatted output and produce plots
}
\usage{
\method{print}{kin.cohort}(x, descriptive = TRUE, cumrisk = TRUE, hazard = FALSE, survival = FALSE, 
        logrank = TRUE, HR = TRUE, digits = 5, \dots)

\method{print}{kin.cohort.boot}(x, cumrisk = TRUE, hazard = TRUE, HR = TRUE, conf = 0.95, 
        digits = 5, show = TRUE, \dots)

\method{plot}{kin.cohort}(x, what = c("cr", "hr", "crr"), min.age = min(x$knots), 
      max.age = max(x$knots), max.y, type, add=FALSE, color, line,  \dots)

\method{plot}{kin.cohort.boot}(x, conf = 0.95, what = c("cr", "hr", "crr"), min.age = min(x$knots), 
      max.age = max(x$knots), age.start = 0, max.y, type, median = FALSE, add = FALSE, color, line, \dots)
}
\arguments{
  \item{x}{ object to be printed or plotted }
  \item{descriptive}{ print table with number of events and person-years}
  \item{cumrisk}{print cumulative risk}
  \item{hazard}{print hazard}
  \item{survival}{print survival}
  \item{HR}{print harard ratios}
  \item{logrank}{print logrank p value}
  \item{digits}{digits for rounding}
  \item{show}{do not print}
  \item{conf}{ coverage for confidence intervals }
  \item{what}{ type of plot desired: cumulative risk ("cr"), hazard ratio ("hr", for marginal method only), cumulative risk ratio ("crr", for moments method only) }
  \item{min.age}{ Minimal age for plots }
  \item{max.age}{ Maximal age for plots }
  \item{age.start}{ initial age value  }
  \item{max.y}{ Max value for y axis }
  \item{type}{ type of line in plots }
  \item{add}{ If TRUE, then lines are added to current plot. Useful to compare analyses.}
  \item{color}{ change line colors using a vector of values }
  \item{line}{ change line pattern using a vector of values }
  \item{median}{ plot median of bootstrap samples instead of point estimates  }
  \item{\dots}{ additional arguments for print or plot }
}
\details{
Specific output and plot types can be selected with arguments
}
\keyword{methods}
\keyword{print}