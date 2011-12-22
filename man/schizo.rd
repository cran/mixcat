\name{schizo}

\alias{schizo}

\docType{data}

\title{National Institute of Mental Health shizophrenia study}

\description{Schizophrenia data from a randomized controlled trial
with patients assigned to either drug or placebo group.
"Severity of Illness" was measured, at weeks 0,1,...6, on a four category ordered scale: 
1. normal or borderline mentally ill, 2. mildly or moderately ill, 3. markedly ill, and 
4. severely or among the most extremely ill. Most of the observations where made
on weeks 0,1,3, and 6.}

\usage{data(schizo)}

\format{
  A data frame with 1603 observations on 437 subjects. Four numerical vectors contain information on
  \describe{
    \item{\code{id}}{patient ID.}
    \item{\code{y}}{ordinal response on a 4 category scale.}
    \item{\code{trt}}{treatment indicator: 1 for drug, 0 for placebo.}
    \item{\code{wk}}{week.}
  }
}

\source{http://tigger.uic.edu/~hedeker/ml.html}

\references{Hedeker, D. and Gibbons, R. (2006). Longitudinal Data Analysis. Wiley, Palo Alto, CA.}

\keyword{datasets}
