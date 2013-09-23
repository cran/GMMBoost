\name{bGLMM}
\alias{bGLMM}
\concept{bGLMM}
\title{Fit Generalized Mixed-Effects Models}
\description{
  Fit a linear mixed model or a generalized linear mixed model.
}

\usage{
bGLMM(fix=formula, rnd=formula, data, family = NULL, control = list())
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
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}. If \code{family} is missing then a
    linear mixed model is fit; otherwise a generalized linear mixed
    model is fit.}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{bGLMMControl}. Defaults to an empty list.}
}
\value{Generic functions such as \code{print}, \code{predict} and \code{summary} have methods to show the results of the fit.
The \code{predict} function uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor).


   \item{call}{a list containing an image of the \code{bGLMM} call that produced the object.}  
  \item{coefficients}{a vector containing the estimated fixed effects}
  \item{ranef}{a vector containing the estimated random effects.}
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
}


\author{
Andreas Groll \email{andreas.groll@stat.uni-muenchen.de}
}

\references{
Tutz, G. and A. Groll (2010). Generalized linear mixed models based on boosting. In
T. Kneib and G. Tutz (Eds.), \emph{Statistical Modelling and Regression Structures - Festschrift
in the Honour of Ludwig Fahrmeir.} Physica.
}


\seealso{
  \code{\link{bGLMMControl}} 
}
\examples{    

data("soccer")
## linear mixed models
lm1 <- bGLMM(points ~ transfer.spendings + I(transfer.spendings^2)
       + ave.unfair.score + transfer.receits + ball.possession
       + tackles + ave.attend + sold.out, rnd = list(team=~1), data = soccer)
      
lm2 <- bGLMM(points~transfer.spendings + I(transfer.spendings^2)
       + ave.unfair.score + transfer.receits + ball.possession
       + tackles + ave.attend + sold.out, rnd = list(team=~1 + ave.attend), 
       data = soccer, control = list(steps=10, lin=c("(Intercept)","ave.attend"), 
       method="REML", nue=1, sel.method="bic"))

## linear mixed models with categorical covariates
lm3 <- bGLMM(points ~ transfer.spendings + I(transfer.spendings^2)
       + as.factor(red.card) + as.factor(yellow.red.card) 
       + transfer.receits + ball.possession + tackles + ave.attend
       + sold.out, rnd = list(team=~1), data = soccer, control = list(steps=10))


## generalized linear mixed model
glm1 <- bGLMM(points~transfer.spendings  + I(transfer.spendings^2)
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = poisson(link = log), data = soccer, 
        control = list(start=c(5,rep(0,31))))

   }
\keyword{models}
\keyword{methods}
