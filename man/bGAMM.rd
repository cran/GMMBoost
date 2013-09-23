\name{bGAMM}
\alias{bGAMM}
\concept{bGAMM}
\title{Fit Generalized Semiparametric Mixed-Effects Models}
\description{
  Fit a semiparametric mixed model or a generalized semiparametric mixed model.
}

\usage{
bGAMM(fix=formula, add=formula, rnd=formula, 
      data, lambda, family = NULL, control = list())
}     
\arguments{
  \item{fix}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right. For categorical covariables use \code{as.factor(.)} in the formula. 
    Note, that the corresponding dummies are treated as a group and are updated blockwise}  
  \item{add}{a one-sided linear formula object describing the
    additive part of the model, with the additive terms on the right side of a
    \code{~} operator, separated by \code{+} operators. The smooth terms 
    are expanded in B-spline basis functions, with a difference penalty apllied on adjacent spline coefficients.}  
  \item{rnd}{a two-sided linear formula object describing the
    random-effects part of the model, with the grouping factor on the left of a
    \code{~} operator and the random terms, separated by \code{+} operators, on
    the right.}
  \item{data}{the data frame containing the variables named in
    \code{formula}.}
 \item{lambda}{the smoothing parameter that controls the smoothness of the additive terms. 
 The optimal smoothing parameter is a tuning parameter of the procedure that has to be determined, 
 e.g. by use of information criteria or cross validation.}
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}. If \code{family} is missing then a
    linear mixed model is fit; otherwise a generalized linear mixed
    model is fit.}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{bGAMMControl}. Defaults to an empty list.}
}
\value{Generic functions such as \code{print}, \code{predict}, \code{summary} and \code{plot} have methods to show the results of the fit. 
The \code{predict} function uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor).
The \code{plot} function shows the estimated smooth functions. Single functions can be specified by a suitable vector in the \code{which} argument. 
Default is \code{which=Null} and all smooth functions (up to a maximum of nine) are shown.


  \item{call}{a list containing an image of the \code{bGLMM} call that produced the object.}  
  \item{coefficients}{a vector containing the estimated fixed effects}
  \item{ranef}{a vector containing the estimated random effects.}
  \item{spline.weights}{a vector containing the estimated spline coefficients.}
  \item{StdDev}{a scalar or matrix containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively.}
  \item{fitted.values}{a vector of fitted values.}
  \item{phi}{estimated scale parameter, if \code{overdispersion=TRUE} is used. Otherwise, it is equal to one.}  
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
  \item{smootherror}{a matrix with pointwise standard errors for the smooth function estimates.}
}
         

\author{
Andreas Groll \email{andreas.groll@stat.uni-muenchen.de}
}

\references{
Groll, A. and G. Tutz (2012). Regularization for Generalized Additive Mixed Models by Likelihood-Based Boosting. 
\emph{Methods of Information in Medicine} 51(2), 168--177.
}


\seealso{
  \code{\link{bGAMMControl}} 
}
\examples{  
data("soccer")

gamm1 <- bGAMM(points ~ ball.possession + tackles,
         ~ transfer.spendings + transfer.receits 
         + unfair.score + ave.attend + sold.out,
         rnd = list(team=~1), data = soccer, lambda = 1e+5,
         family = poisson(link = log), control = list(steps=200, overdispersion=TRUE,
         start=c(1,rep(0,25))))

plot(gamm1)

# see also demo("bGAMM-soccer")
}
\keyword{models}
\keyword{methods}
