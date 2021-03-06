\name{bGAMMControl}
\alias{bGAMMControl}
\concept{bGAMMControl}
\title{Control Values for \code{bGAMM} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{bGAMM} function.
}

\usage{
bGAMMControl(nue=0.1,add.fix=NULL,start=NULL,q_start=NULL, 
             OPT=TRUE,nbasis=20,spline.degree=3,
             diff.ord=2,sel.method="aic",steps=500,
             method="EM",overdispersion=FALSE)  
} 
    
\arguments{
  \item{nue}{weakness of the learner. Choose 0 < nue =< 1. Default is 0.1.}  
  \item{add.fix}{a vector specifying smooth terms, which are excluded from selection.}
  \item{start}{a vector containing starting values for fixed and random effects of suitable length. Default is a vector full of zeros.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal.}
  \item{OPT}{logical scalar. When \code{TRUE} the estimates at the optimal number of boosting steps, chosen by information criteria, are derived. If \code{FALSE},
    the estimates at the maximal number of boosting steps are derived. Default is \code{TRUE}.}
  \item{nbasis}{the number of b-spline basis functions for the modeling of smooth terms. Default is 20.}
  \item{spline.degree}{the degree of the B-spline polynomials. Default is 3.}
  \item{diff.ord}{the order of the difference penalty; must be lower than the degree of the B-spline polynomials (see previous argument). Default is 2.}
  \item{sel.method}{two different information criteria, "aic" or "bic", can be chosen, on which the selection step is based on. Default is "aic".}
  \item{steps}{the number of boosting interations. Default is 500.}
  \item{method}{two methods for the computation of the random-effects variance-covariance parameter estimates can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{bobyqa} function for optimization.
  Default is \code{EM}.}
  \item{overdispersion}{logical scalar. If \code{FALSE}, no scale parameter is derived, if \code{TRUE}, in each boosting iteration a scale parameter is estimated by use of Pearson residuals. 
  This can be used to fit overdispersed Poisson models. Default is \code{FALSE}.}
}

\value{
  a list with components for each of the possible arguments.
}

\author{
Andreas Groll \email{andreas.groll@stat.uni-muenchen.de}
}

\seealso{
  \code{\link{bGAMM}}, \code{\link[minqa]{bobyqa}}
}

\examples{
# decrease the maximum number of boosting iterations 
# and use BIC for selection
bGAMMControl(steps = 100, sel.method = "BIC")
}
