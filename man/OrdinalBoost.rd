\name{OrdinalBoost}
\alias{OrdinalBoost}
\concept{OrdinalBoost}
\title{Fit Generalized Mixed-Effects Models}
\description{
  Fit a generalized linear mixed model with ordinal response.
}

\usage{
OrdinalBoost(fix=formula, rnd=formula, data,model="sequential",control=list())
}     
\arguments{
  \item{fix}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right. For categorical covariables use \code{as.factor(.)} in the formula. 
    Note, that the corresponding dummies are treated as a group and are updated blockwise}  
  \item{rnd}{a two-sided linear formula object describing the
    random-effects part of the model, with the grouping factor on the left of a
    \code{~} operator and the random terms, separated by \code{+} operators, on
    the right.}
  \item{data}{the data frame containing the variables named in
    \code{formula}.}
  \item{model}{
    Two models for repeatedly assessed ordinal scores, based on the threshold concept, are available, the "sequential" and the "cumulative" model. Default is "sequential".}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{OrdinalBoostControl}. Defaults to an empty list.}
}
\value{Generic functions such as \code{print}, \code{predict} and \code{summary} have methods to show the results of the fit. The \code{predict} function shows the estimated probabilities for the different categories
 for each observation, either for the data set of the \code{OrdinalBoost} object or for \code{newdata}. Default is \code{newdata=Null}. 
 It uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor).


   \item{call}{a list containing an image of the \code{OrdinalBoost} call that produced the object.}  
  \item{coefficients}{a vector containing the estimated fixed effects}
  \item{ranef}{a vector containing the estimated random effects.}
  \item{StdDev}{a scalar or matrix containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively.}
  \item{fitted.values}{a vector of fitted values.}
  \item{HatMatrix}{hat matrix corresponding to the final fit.}
  \item{IC}{a matrix containing the evaluated information criterion for the different covariates (columns) and for each boosting iteration (rows).}
  \item{IC_sel}{a vector containing the evaluated information criterion for the selected covariate at different boosting iterations.}
  \item{components}{a vector containing the selected components at different boosting iterations.}
  \item{opt}{number of optimal boosting steps with respect to AIC or BIC, respectively, if \code{OPT=TRUE}. Otherwise, \code{opt} is equal to the number of iterations.
  Note, that the boosting algorithm is also stopped, if it has converged with respect to the parameter estimates \code{[coefficients,ranef]} or with respect to the \code{IC_sel}.}
  \item{Deltamatrix}{a matrix containing the estimates of fixed and random effects (columns) for each boosting iteration (rows).}
  \item{Q_long}{a list containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively, for each boosting iteration.}
  \item{fixerror}{a vector with standrad errors for the fixed effects.}
  \item{ranerror}{a vector with standrad errors for the random effects.}
}


\author{
Andreas Groll \email{andreas.groll@stat.uni-muenchen.de}
}

\references{
Tutz, G. and A. Groll (2012). Likelihood-based boosting in binary and ordinal random effects models. 
\emph{Journal of Computational and Graphical Statistics}. To appear.
}


\seealso{
  \code{\link{OrdinalBoostControl}} 
}
\examples{ 
\dontrun{
data(knee)

# fit a sequential model
# (here only one step is performed in order to
# save computational time)

glm1 <- OrdinalBoost(pain ~ time + th + age + sex, rnd = list(id=~1),
        data = knee, model = "sequential", control = list(steps=1))

# see also demo("OrdinalBoost-knee") for more extensive examples
}}              
\keyword{models}
\keyword{methods}
