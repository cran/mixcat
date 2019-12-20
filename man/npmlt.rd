\name{npmlt}
\alias{npmlt}

\title{Mixed effects cumulative link and logistic regression models}

\description{Fits cumulative logit and baseline logit and link mixed effects regression models with non- parametric distribution for the random effects.}

\usage{
npmlt(formula, formula.npo=~1, random=~1, id, k=1, eps=0.0001,
      start.int=NULL, start.reg=NULL, start.mp=NULL,
      start.m=NULL, link="clogit",
      EB=FALSE, maxit=500, na.rm=TRUE, tol=0.0001)
}

\arguments{
\item{formula}{a formula defining the response and the fixed, proportional odds, effects part of the model,
               e.g. \code{y ~ x}.}

\item{formula.npo}{a formula defining non proportional odds variables of the model.
                     A response is not needed as it has been provided in \code{formula}.
                     Intercepts need not be provided as they are always non proportional.
                     Variables in \code{formula.npo} must be a subset of the variables
                     that appear in the right hand side of \code{formula}, e.g. \code{ ~ x}.}

\item{random}{a formula defining the random part of the model. For instance, \code{random = ~1}
                defines a random intercept model, while \code{random = ~1+x} defines a model
                with random intercept and random slope for the variable \code{x}. If
                argument \code{k=1}, the resulting model is a fixed effects model (see below).
                Variables in \code{random} must be a subset of the variables that appear
                in the right hand side of \code{formula}.}

\item{id}{a factor that defines the primary sampling units, e.g.
            groups, clusters, classes, or individuals in longitudinal studies. These sampling units
            have their own random coefficient, as defined in \code{random}.
            If argument \code{id} is missing it is taken to be \code{id=seq(N)},
            where N is the total number of observations,
            suitable for overdispersed independent multinomial data.}

\item{k}{the number of mass points and masses for the non-parametric (discrete) random effects distribution.
           If \code{k=1} the function fits a fixed effects models, regerdless of the \code{random} specification,
           as with \code{k=1} the random effects distribution is degenerate at zero.}

\item{eps}{positive convergence tolerance \eqn{epsilon}. Convergence is declared when the maximum
             of the absolute value of the score vector is less than \eqn{epsilon}.}

\item{start.int}{a vector of length (number of categories minus one) with the starting values the
                   fixed intercept(s).}

\item{start.reg}{a vector with the starting values for the regression coefficients. One
                   starting value for the proportional odds effects and (number of categories minus one)
                   starting values for the non proportional effects, in the same order as they appear
                   in \code{formula}.}

\item{start.mp}{starting values for the mass points of the random effects distribution in the form:
                  (\code{k} starting values for the intercepts, \code{k} starting values for the
                  first random slope,...).}

\item{start.m}{starting values for the masses of the random effects distribution: a vector of
                 length \code{k} with non-negative elements that sum to 1.}

\item{link}{for a cumulative logit model set \code{link="clogit"} (default).
              For a baseline logit model, set \code{link="blogit"}. Baseline category
              is the last category.}

\item{EB}{if \code{EB=TRUE} the empirical Bayes estimates of the random effects are calculated
            and stored in the component \code{eBayes}. Further, fitted values of the linear predictor
            (stored in the component \code{fitted}) and fitted probabilities (stored in object \code{prob}) are
            obtained at the empirical Bayes estimates of the random
            effects. Otherwise, if \code{EB=FALSE} (default), empirical Bayes estimates are not
            calculated and fitted values of the linear predictors and probabilities are calculated at the zero
            value of the random effects.}

\item{maxit}{integer giving the maximal number of iterations of the fitting algorithm until convergence.
               By default this number is set to 500.}

\item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}

\item{tol}{positive tolerance level used for calculating generalised inverses (g-inverses). Consider matrix
             \eqn{A = P D P^T}, where \eqn{D=Diag\{eigen_i\}} is diagonal with entries the eigen values of \eqn{A}.
             Its g-inverse is calculated as \eqn{A^{-} = P D^{-} P^T}, where \eqn{D^{-}} is diagonal with entries
             \eqn{1/eigen_i} if \eqn{eigen_i > tol}, and \eqn{0} otherwise.}
}

\details{Maximizing a likelihood over an unspecified random effects distribution results in a
         discrete mass point estimate of this distribution (Laird, 1978; Lindsay, 1983).
         Thus, the terms `non-parametric' (NP) and `discrete' random effects distribution are used here
         interchangeably. Function \code{npmlt} allows the user to choose the number \code{k} of mass
         points/masses of the discrete distribution, a choice that should be based on the log-likelihood.
         Note that the mean of the NP distribution is constrained to be zero and thus for \code{k=1}
         the fitted model is equivalent to a fixed effects model. For \code{k>1} and a random slope in the model,
         the mass points are bivariate with a component that corresponds to the intercept and another
         that corresponds to the slope.

         General treatments of non-parametric modeling can be found in Aitkin, M. (1999) and Aitkin et al. (2009).
         For more details on multinomial data see Hartzel et al (2001).

         The response variable \code{y} can be binary or multinomial. A binary response should
         take values 1 and 2, and the function \code{npmlt} will model the probability of 1. For an ordinal
         response, taking values \eqn{1,\dots,q}, a cumulative logit model can be fit. Ignoring the random effects,
         such a model, with formula \code{y~x}, takes the form
         \deqn{log \frac{P(Y \le r)}{1-P(Y \le r)}=\beta_r + \gamma x,}{log\{P(Y \le r)/[1-P(Y \le r)]\}=
         \beta_r + \gamma x,}
         where \eqn{\beta_r, r=1,\dots,q-1}, are the cut-points and \eqn{\gamma} is the slope.
         Further, if argument formula.npo is specified as \code{~x}, the model becomes
         \deqn{log \frac{P(Y \le r)}{1-P(Y \le r)}=\beta_r + \gamma_r x,}{log\{P(Y \le r)/[1-P(Y \le r)]\}=
         \beta_r + \gamma_r x.}
         Similarly, for a nominal response with q categories, a baseline logit model can be fit.
         The fixed effects part of the model, \code{y~x}, takes the form,
         \deqn{log \frac{P(Y=r)}{P(Y=q)} = \beta_r + \gamma x,}{log[P(Y=r)/P(Y=q)] = \beta_r + \gamma x,}
         where \eqn{r=1,\dots,q-1.} Again, formula.npo can be specified as \code{~x}, in which
         case slope \eqn{\gamma} will be replaced by category specific slopes, \eqn{\gamma_r}.

         The user is provided with the option of specifying starting values for some or all the model
         parameters. This option allows for starting the algorithm at different starting points, in order
         to ensure that it has convered to the point of maximum likelihood. Further, if the fitting algorithm
         fails, the user can start by fitting a less complex model and use the estimates of this model
         as starting values for the more complex one.

         With reference to the \code{tol} argument, the fitting
         algorithm calculates g-inverses of two matrices: 1. the information matrix of the model,
         and 2. the covariance matrix of multinomial proportions. The covariance matrix of a
         multinomial proportion \eqn{p} of length \eqn{q} is calculated as \eqn{Diag\{p*\} -p* p*^T},
         where \eqn{p*} is of length \eqn{q-1}. A g-inverse for this matrix is needed because elements
         of \eqn{p*} can become zero or one.}

\value{The function \code{npmlt} returns an object of class \sQuote{npmreg}, a list containing at least the following components:

\item{call}{the matched call.}
\item{formula}{the formula supplied.}
\item{formula.npo}{the formula for the non proportional odds supplied.}
\item{random}{the random effects formula supplied.}
\item{coefficients}{a named vector of regression coefficients.}
\item{mass.points}{a vector or a table that contains the mass point estimates.}
\item{masses}{the masses (probabilities) corresponding to the mass points.}
\item{vcvremat}{the estimated variance-covariance matrix of the random effects.}
\item{var.cor.mat}{the estimated variance-covariance matrix of the random effects, with the upper triangular covariances replaced by the corresponding correlations.}
\item{m2LogL}{minus twice the maximized log-likelihood of the chosen model.}
\item{SE.coefficients}{a named vector of standard errors of the estimated regression coefficients.}
\item{SE.mass.points}{a vector or a table that contains the the standard errors of the estimated mass points.}
\item{SE.masses}{the standard errors of the estimated masses.}
\item{VRESE}{the standard errors of the estimates of the variances of random effects.}
\item{CVmat}{the inverse of the observed information matrix of the model.}
\item{eBayes}{if \code{EB=TRUE} it contains the empirical Bayes estimates of the random effects. Otherwise it  contains vector(s) of zeros.}
\item{fitted}{the fitted values of the linear predictors computed at the empirical Bayes estimates of the random effects, if \code{EB=TRUE}. Otherwise, if \code{EB=FALSE} (default) these fitted values are computed at the zero value of the random effects.}
\item{prob}{the estimated probabilities of observing a response at one of the categories. These probabilities are computed at the empirical Bayes estimates of the random effects, if \code{EB=TRUE}. If \code{EB=FALSE} (default) these estimated probabilities are computed at the zero value of the random effects.}
\item{nrp}{number of random slopes specified.}
\item{iter}{the number of iterations of the fitting algorithm.}
\item{maxit}{the maximal allowed number of iterations of the fitting algorithm until convergence.}
\item{flagcvm}{last iteration at which eigenvalue(s) of covariance matrix of multinomial variable were less than \code{tol} argument.}
\item{flaginfo}{last iteration at which eigenvalue(s) of model information matrix were less than \code{tol} argument.}
}

\references{
Aitkin, M. (1999). A general maximum likelihood analysis of variance components in generalized linear models. Biometrics 55, 117-128.

Aitkin, M., Francis, B., Hinde, J., and Darnell, R. (2009). Statistical Modelling in R. Oxford Statistical Science Series, Oxford, UK.

Hedeker, D. and Gibbons, R. (2006). Longitudinal Data Analysis. Wiley, Palo Alto, CA.

Hartzel, J., Agresti, A., and Caffo, B. (2001). Multinomial logit random effects models. Statistical Modelling, 1(2), 81-102.

Laird, N. (1978). Nonparametric maximum likelihood estimation of a mixing distribution. Journal of the
American Statistical Association, 73, 805-811.

Lindsay, B. G. (1983). The geometry of mixture likelihoods, Part II: The exponential family. The Annals
of Statistics, 11, 783-792.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\seealso{\code{\link{summary.npmreg}}}

\examples{
data(schizo)
attach(schizo)
\donttest{
npmlt(y~trt*sqrt(wk),formula.npo=~trt,random=~1+trt,id=id,k=2,EB=FALSE)
}
}

\keyword{models}
\keyword{regression}
