\name{summary.npmreg}

\alias{summary.npmreg}
\alias{print.npmreg}

\title{Summarizing mixed multinomial regression model fits}

\description{\code{summary} and \code{print} methods for objects of  type \code{npmreg}.}

\usage{
summary.npmreg(object,digits = max(3, getOption("digits") - 3),...)
print.npmreg(x,digits = max(3, getOption("digits") - 3),...)}

\arguments{
  \item{object}{an object of class \code{npmreg}.}
  \item{x}{an object of class \code{npmreg}.}
  \item{digits}{the minimum number of significant digits to be printed in values.}
  \item{\dots}{further arguments, which will mostly be ignored.}
}

\details{The function \code{npmlt} returns an object of class "npmreg".
The function \code{summary} (i.e., \code{summary.npmreg}) can be used to obtain or print a summary of the results,
and the function \code{print} (i.e., \code{print.npmreg}) to print the results.}

\value{Summary or print.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com} and John Hinde}

\seealso{\code{\link{npmlt}}}

\examples{
data(schizo)
attach(schizo)
fit1<-npmlt(y~trt*sqrt(wk),formula.npo=~trt,random=~1,id=id,k=2)
print(fit1)
summary(fit1)
}

\keyword{models}
\keyword{regression}
