bGLMMControl<-function(nue=0.1,lin="(Intercept)",start=NULL,q_start=NULL, OPT=TRUE,sel.method="aic",steps=500,method="EM",overdispersion=FALSE)
{                       
if (lin[1]!="(Intercept)")
lin<-c("(Intercept)",lin)

list(nue = nue, lin = lin, start = start, q_start = q_start, OPT = OPT,
        sel.method = sel.method, steps = steps, method = method,
        overdispersion = overdispersion)
}
