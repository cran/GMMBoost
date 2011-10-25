\name{OrdinalBoostControl}
\alias{OrdinalBoostControl}
\concept{OrdinalBoostControl}
\title{Control Values for \code{OrdinalBoost} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{bGLMM} function.
}

\usage{
OrdinalBoostControl(nue=0.1,lin=NULL,katvar=NULL,start=NULL,q_start=NULL, OPT=TRUE, 
                    sel.method="aic",steps=100,method="EM",maxIter=500)
} 
    
\arguments{
  \item{nue}{weakness of the learner. Choose 0 < nue =< 1. Default is 0.1.}  
  \item{lin}{a vector specifying fixed effects, which are excluded from selection.}
  \item{katvar}{a vector specifying category-specific covariates, which are also excluded from selection.}
  \item{start}{a vector containing starting values for fixed and random effects of suitable length. Default is a vector full of zeros.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal.}
  \item{OPT}{logical scalar. When \code{TRUE} the estimates at the optimal number of boosting steps, chosen by information criteria, are derived. If \code{FALSE},
    the estimates at the maximal number of boosting steps are derived. Default is \code{TRUE}.}
  \item{sel.method}{two different information criteria, "aic" or "bic", can be chosen, on which the selection step is based on. Default is "aic".}
  \item{steps}{the number of boosting interations. Default is 100.}
  \item{method}{two methods for the computation of the random-effects variance-covariance parameter estimates can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{bobyqa} function for optimization.
  Default is \code{EM}.}
  \item{maxIter}{the number of interations for the final Fisher scoring reestimation procedure. Default is 500.}
}

\value{
a list with components for each of the possible arguments.
}

\author{
Andreas Groll \email{andreas.groll@stat.uni-muenchen.de}
}

\seealso{
  \code{\link{OrdinalBoost}}, \code{\link[minqa]{bobyqa}}
}

\examples{
# decrease the maximum number of boosting iterations 
# and use BIC for selection
OrdinalBoostControl(steps = 10, sel.method = "BIC")
}
