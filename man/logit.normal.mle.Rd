\name{logit.normal.mle}
\alias{logit.normal.mle}

\title{ Logistic-Normal Maximum Likelihood Estimation for Binary Data }
\description{

\code{\link{logit.normal.mle}} allows for a more general model than the simple logistic-normal model since the scale parameter, \code{sigma}, may depend on cluster level characteristics.

}
\usage{
logit.normal.mle(meanmodel, logSigma, id, n=NULL, beta=NULL, alpha=NULL, model="marginal", lambda=0.0, r = 20, 
			maxits=50, tol = 1e-3, data = sys.frame(sys.parent()) )

}

\arguments{
  \item{meanmodel}{
a symbolic description of the marginal model to be fit that generally takes the form \code{y ~ x}.  Further details of model specification are provided below.
}
  \item{logSigma}{
a symbolic description of the model used to estimate the dependence of \code{sigma} on covariates \code{z}, where \code{z} is a subset of covariates \code{x}. In general, \code{random} has the form \code{ ~ z}.
}
  \item{id}{Cluster identification variable.}
  \item{n}{Number of responses per binomial trial.}
  \item{beta}{Initial estimate of mean parameter.}
  \item{alpha}{Initial estimate of log variance component parameter. }
  \item{model}{"conditional" for classic GLMM, "marginal" for the marginalized model. }
  \item{lambda}{A likelihood penalty parameter (>= 0) for \code{alpha} }
  \item{r}{Number of Gauss-Hermite quadrature points. The user may choose r=3, 5, 10, 20, or 50.  The default value is r=20.}
  \item{maxits}{ Maximum number of iterations for convergence.  The default is maxits=50.}
  \item{tol}{Convergence criterion that specifies the absolute change in parameter estimates. }
  \item{data}{an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(marginal)}, typically the environment from which \code{logit.normal.mle} is called.}
}
\details{
\code{logit.model.mle} assumes that longitudinal binary measurements of \code{y} and possibly time-dependent exogenous covariates \code{x} are collected at times t(i,1), t(i,2),..., t(i,n(i)) for i = 1,2,...,N subjects or clusters.  Each of the clusters or individuals (\code{id}) need not have measurements at every interval. 

The \code{meanmodel} is specified symbolically as a formula.  A typical marginal model has the form \code{y ~ x}.  The covariates \code{x} are a series of terms separated by \code{+} which specify the marginal linear predictor for the longitudinal (binary) response vector \code{y}.  The parameter(s) \code{beta} denote the contrast in log odds of success for subgroups defined by covariates \code{x(t)}  The user may provide initial estimate(s) for \code{beta}.     

The parameter \code{alpha} describes how the individual-level heterogeneity \code{sigma_i} in the log odds varies as a function of covariates \code{z_i}.  The random model statement (\code{logSigma}) has the form \code{ ~ z}.  The covariates \code{z} are a subset of \code{x}, possibly just an intercept.  The user may provide initial estimates(s) for \code{alpha}.

All formulas have an implied intercept term.  To remove this, include \code{-1} on the right-hand side of the formula statement.  See \code{\link{formula}} for more details. 

}
\value{
Returns an object of class \code{logit.normal.mle}. The function \code{print.logit.normal.mle} may used to obtain and print a summary of the results.  See below for an example.
}
\references{ 
Diggle, P.J. Heagerty, P. Liang, K.Y. and Zeger, S.L. (2002) "Analysis of Longitudinal Data", 2nd Edition, Oxford University Press.

Heagerty, P.J. (1999) "Marginally Specified Logistic-Normal Models for Longitudinal Binary Data", Biometrics 55, 688-98.

}
}
\author{ Patrick Heagerty \email{heagerty@u.washington.edu} }


\seealso{\code{\link{print.logit.normal.mle}} }
\examples{
#
## Madras data example:
#
data(madras)
attach(madras)
model1 <- logit.normal.mle(meanmodel = y ~ gender+month+age+monthXage, logSigma = ~ 1 + age, id=id, model="marginal", data=madras)
		print.logit.normal.mle(model1)
}
\keyword{}

