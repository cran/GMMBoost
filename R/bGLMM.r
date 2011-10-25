#############################################################################################################################################               
###################################################### main boosting function #################################################################
est.bGLMM<-function(fix,rnd,data,family=gaussian(link = "identity"),control=list())
{                      
y <- model.response(model.frame(fix, data))

X <- model.matrix(fix, data)

very.old.names<-attr(terms(fix),"term.labels")
factor.names<-character()
for (i in 1:length(very.old.names))
{
if (substr(very.old.names[i],1,9)=="as.factor")
factor.names<-c(factor.names,very.old.names[i])
}

factor.list<-list()
if(length(factor.names)>0)
{
spl<-strsplit(factor.names,"\\(")

categ.names<-character()
for(uz in 1:length(spl))
categ.names<-c(categ.names,spl[[uz]][2])

spl2<-strsplit(categ.names,"\\)")
categ.names2<-character()
for(uz in 1:length(spl2))
categ.names2<-c(categ.names2,spl2[[uz]])
factor.names<-categ.names2

for (i in 1:length(factor.names))
factor.list[[i]]<-levels(as.factor(data[,factor.names[i]]))
}    


old.names<-attr(X,"dimnames")[[2]]

rndformula <- as.character(rnd)

trmsrnd <- terms(rnd[[1]])
newrndfrml <- "~ -1"
newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")

if(length(attr(trmsrnd, "variables"))>1)
{
newrndfrml <- paste(newrndfrml,  
         paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                     paste(lbl, names(rnd)[1], sep=":")
                 }), collapse=" + "), sep="+") }


W_start <- model.matrix(formula(newrndfrml), data)

rnlabels<-terms(formula(newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
k<-table(data[,colnames(data)==(names(rnd)[1])])   
n<-length(k)
s<-dim(W_start)[2]/n


if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}


print(paste("Iteration ", 1,sep=""))
control<-do.call(bGLMMControl, control)

if(sum(substr(control$lin,1,9)=="as.factor")>0)
{
group0<-substr(control$lin,1,9)=="as.factor"
spl0<-strsplit(control$lin[group0],"\\)")

control$lin<-control$lin[!group0]
for(ur in 1:length(spl0))
control$lin<-c(control$lin,old.names[is.element(substr(old.names,1,nchar(spl0[[ur]])),spl0[[ur]])])
}



if(length(control$lin)>1)
{
lin.out<-!is.element(very.old.names,control$lin[2:length(control$lin)])
very.old.names<-very.old.names[lin.out]
}


group<-substr(very.old.names,1,9)=="as.factor"

if(sum(group)>0)
{
spl<-strsplit(very.old.names[group],"\\(")

categ.names<-character()
for(uz in 1:length(spl))
categ.names<-c(categ.names,spl[[uz]][2])

spl2<-strsplit(categ.names,"\\)")
categ.names2<-character()
for(uz in 1:length(spl2))
categ.names2<-c(categ.names2,spl2[[uz]])

block<-numeric()
posi<-1
for(ip in 1:length(group))
{
  if(!group[ip])
  {
  block<-c(block,1)
  }else{
  block<-c(block,length(levels(as.factor(data[,categ.names2[posi]])))-1)
  posi<-posi+1
  }
}}else{
block<-rep(1,length(group))  
}

if(length(control$start)==0)
control$start<-c(rep(0,(length(control$lin)+n*s)))

if(length(control$q_start)==0)
{
control$q_start<-rep(0.1,s)
if(s>1)
control$q_start<-diag(control$q_start)
}

beta_null<-control$start[1:length(control$lin)]
ranef_null<-control$start[(length(control$lin)+1):(length(control$lin)+n*s)]
q_start<-control$q_start

N<-length(y)
if (all (X==0))
{
X<-rep(1,N)
X<-as.matrix(X)
}
# add Intercept 
if (!all (X[,1]==1))
X<-cbind(1,X)

Z_fastalles<-X

no.sel<-is.element(attr(X,"dimnames")[[2]],control$lin)

U<-X[,!no.sel]
X<-as.matrix(X[,no.sel])

lin<-length(control$lin)

m<-length(block)
m_alt<-ncol(U)              

lin0<-sum(beta_null!=0)
if(s==1)
{
Q_start<-diag(q_start^2,s)
p_start<-c(rep(0,lin0),rep((q_start^2)^(-1),n*s))
P_start<-diag(p_start)
}else{
Q_start<-q_start
P_start<-matrix(0,lin0+n*s,lin0+n*s)
for(jf in 1:n)
P_start[(lin0+(jf-1)*s+1):(lin0+jf*s),(lin0+(jf-1)*s+1):(lin0+jf*s)]<-chol2inv(chol(q_start))
}





if(lin>1)
{
Eta_start<-X%*%beta_null+W%*%ranef_null
}else{
Eta_start<-rep(beta_null,N)+W%*%ranef_null
}

D_start<-as.vector(family$mu.eta(Eta_start))
Mu_start<-as.vector(family$linkinv(Eta_start))
Sigma_start<-as.vector(family$variance(Mu_start))
W0_inv<-D_start*1/Sigma_start*D_start


Z_start<-cbind(X[,beta_null!=0],W)
M0<-Z_start%*%chol2inv(chol((t(Z_start)%*%(Z_start*W0_inv))+P_start))%*%t(Z_start*W0_inv)

blocksum<-rep(0,m)
for (jg in 1:m)
blocksum[jg]<-sum(block[1:jg])


Mu<-as.vector(family$linkinv(Eta_start))

Z_alles<-cbind(X,U,W)

#########################  some definitions ################

Z_r<-list()


Z_r[[1]]<-cbind(X,as.matrix(U[,1:block[1]]))
for (r in 2:m)
Z_r[[r]]<-cbind(X,as.matrix(U[,(blocksum[r-1]+1):blocksum[r]]))

Z_r_gross<-list()
for (r in 1:m)
Z_r_gross[[r]]<-cbind(Z_r[[r]],W)


Sigma<-as.vector(family$variance(Mu))

if(control$overdispersion==T)
{
phi_start<-(sum((y-Mu)^2/Mu))/(N-sum(diag(M0)))
Sigma<-Sigma*phi_start
}

########################################################## some definitions ################################################

Delta<-matrix(0,control$steps,(lin+m_alt+s*n))


Standard<-matrix(0,control$steps,(lin+m_alt+s*n))

Standard_ma<-list()

for (ijk in 1:m)
Standard_ma[[ijk]]<-as.numeric()

Delta[1,1:lin]<-beta_null
Delta[1,(lin+m_alt+1):(lin+m_alt+s*n)]<-t(ranef_null)

Delta_start<-Delta[1,]

IC<-matrix(0,nrow=control$steps,ncol=m)
IC_stop<-rep(0,control$steps)
komp<-rep(0,control$steps)

E<-diag(N)

M1<-list()
diag_H1<-rep(0,m)
FinalHat<-list()

phi<-rep(0,control$steps)

Q<-array(0,c(s,s,control$steps+1))
Q[,,1]<-Q_start
######### derive D and score_r and F_r  ###########
D<-as.vector(family$mu.eta(Eta_start))


Delta_r<-list()
for(r in 1:m)
Delta_r[[r]]<-rep(0,lin+block[r]+s*n)


Eta_r<-matrix(0,N,m)
Mu_r<-matrix(0,N,m)

grp<-rep(1,m_alt)
for (r in 2:m_alt)
grp[((r-1)+1):(r)]<-r



if(s==1)
{
P_alles<-c(rep(0,lin+m_alt),rep((Q_start^(-1)),n*s))
P_alles<-diag(P_alles)
}else{
P_alles<-matrix(0,lin+m_alt+n*s,lin+m_alt+n*s)
for(jf in 1:n)
P_alles[(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s),(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s)]<-chol2inv(chol(Q_start))
}


score_alles<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_alles<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P_alles

Inv_r<-list()


l=1
for (r in 1:m)
{      
  if(r==1)
  {
  Inv_r[[1]]<-chol2inv(chol(F_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[1]]%*%score_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]
  }else{
  Inv_r[[r]]<-chol2inv(chol(F_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[r]]%*%score_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]
  }

############# derive predictor Eta_r ########
Eta_r[,r]<-Eta_start+Z_r_gross[[r]]%*%Delta_r[[r]]

############# new  Mu_r ##############
Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))

################### derive Hat-Matrix  ##############
M1[[r]]<-(Z_r_gross[[r]]*sqrt(Sigma*D*1/Sigma*D))%*%Inv_r[[r]]%*%t(Z_r_gross[[r]]*sqrt(D*1/Sigma*D*1/Sigma))
diag_H1[r]<-sum(diag(M0+M1[[r]]%*%(E-M0)))

######## likelihood of information criterion   ########
    if (family$family=="poisson")
    IC[1,r]<-sum(y*log(Mu_r[,r])-Mu_r[,r])#-0.5*pena

    if (family$family=="binomial")
    IC[1,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena

    if (family$family=="gaussian")
    IC[1,r]<-sum(y*Mu_r[,r]-0.5*(Mu_r[,r]^2))#-0.5*pena
    
   ############ AIC ##############
    if (control$sel.method=="aic")
    IC[1,r]<-(-IC[1,r])+diag_H1[r]

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC[1,r]<-(-IC[1,r])+0.5*diag_H1[r]*log(n)
}

   mi_ic<-match(min(IC[1,]),IC[1,])

#################### update ###########################
if(mi_ic==1)   
{
g<-c(rep(T,lin),is.element(grp,1:blocksum[1]),rep(T,s*n))
}else{
g<-c(rep(T,lin),is.element(grp,(blocksum[mi_ic-1]+1):blocksum[mi_ic]),rep(T,s*n))
}

Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]<-control$nue*Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]
Delta[1,g]<-Delta[1,g]+t(Delta_r[[mi_ic]])

Nue_ma<-c(rep(1,lin),rep(control$nue,block[mi_ic]),rep(1,s*n))

Uu<-(E-(cbind(Z_r[[mi_ic]],W)*sqrt(Sigma*D*1/Sigma*D))%*%(Inv_r[[mi_ic]]*Nue_ma)%*%t(cbind(Z_r[[mi_ic]],W)*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)
FinalHat[[1]]<-E-Uu

defre<-sum(diag(FinalHat[[1]]))

komp[1]<-mi_ic

if(mi_ic==1)
{
X1<-cbind(X,as.matrix(U[,1:block[1]]))
}else{
X1<-cbind(X,as.matrix(U[,(blocksum[mi_ic-1]+1):blocksum[mi_ic]]))
}


if(s==1)
{
PP<-rep((Q_start^(-1)),n*s)
PP<-diag(PP)
}else{
PP<-matrix(0,n*s,n*s)
for(jf in 1:n)
PP[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Q_start))
}


R1<-chol2inv(chol((t(X1)%*%(X1*D*1/Sigma*D))-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP))%*%t(t(X1)%*%(W*D*1/Sigma*D)))))
R2<-(t(X1*D*1/Sigma)-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP)))%*%t(W*D*1/Sigma))

Q_mi<-control$nue*R1[(lin+block[mi_ic]):(lin+block[mi_ic]),]%*%R2%*%(E-M0)

Standard_ma[[mi_ic]]<-Q_mi

if(mi_ic==1)
{
Standard[1,(lin+1):(lin+blocksum[1])]<-sqrt(diag(Standard_ma[[1]]%*%t(Standard_ma[[1]]*Sigma)))
}else{
Standard[1,(lin+blocksum[mi_ic-1]+1):(lin+blocksum[mi_ic])]<-sqrt(diag(Standard_ma[[mi_ic]]%*%t(Standard_ma[[mi_ic]]*Sigma)))
}   
Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-as.vector(family$variance(Mu))

D<-as.vector(family$mu.eta(Eta))


if (control$method=="EM")
{
aktuell_vec<-!is.element(Delta[1,],0)
X_aktuell<-Z_alles[,aktuell_vec] 

m_aktuell<-sum(block[is.element(1:m,komp)])

if(s==1)
{
P_gross<-c(rep(0,lin+m_aktuell),rep((Q_start^(-1)),n*s))
P_gross<-diag(P_gross)
}else{
P_gross<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
for(jf in 1:n)
P_gross[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q_start))
}


F_gross<-t(X_aktuell)%*%(X_aktuell*D*1/Sigma*D)+P_gross


############################# update Q ################
InvFisher<-chol2inv(chol(F_gross))

Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[1,(lin+m+1):(lin+m+s)]%*%t(Delta[1,(lin+m+1):(lin+m+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)])

#Q1<-InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,seq(from=lin+1,to=lin+1+(s-1)*n,by=n)]%*%t(Delta[1,seq(from=lin+1,to=lin+1+(s-1)*n,by=n)])
#for (i in 2:n)
#Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,seq(from=lin+i,to=lin+i+(s-1)*n,by=n)]%*%t(Delta[1,seq(from=lin+i,to=lin+i+(s-1)*n,by=n)])

Q1<-1/n*Q1
}else{

Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:(lin+m_alt)]     

aktuell_vec<-!is.element(Delta[1,1:(lin+m_alt)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]


if(s==1)
{
optim.obj<-bobyqa(sqrt(Q_start),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, model="no.model",  lower = 1e-14, upper = 20)
Q1<-as.matrix(optim.obj$par)^2
}else{
q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
upp<-rep(up1,length(q_start_vec))
low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, model="no.model", lower=low,upper=upp))
Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}
}

    if (family$family=="poisson")
    IC_stop[1]<-sum(y*log(Mu)-Mu)

    if (family$family=="binomial")
    IC_stop[1]<-sum(y*log(Mu)+(1-y)*log(1-Mu))

    if (family$family=="gaussian")
    IC_stop[1]<-sum(y*Mu-0.5*(Mu^2))
    
   ############ AIC ##############
    if (control$sel.method=="aic")
    IC_stop[1]<-(-IC_stop[1])+defre

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC_stop[1]<-(-IC_stop[1])+0.5*defre*log(n)

if(control$overdispersion==T)
{
phi[1]<-(sum((y-Mu)^2/Mu))/(N-defre)
Sigma<-Sigma*phi[1]
}



Q[,,2]<-Q1
###############################################################################################################################################
################################################################### Boost ###################################################################
if(control$steps!=1)
{
for (l in 2:control$steps)
{
print(paste("Iteration ", l,sep=""))



if(s==1)
{
P_alles<-c(rep(0,lin+m_alt),rep((Q1^(-1)),n*s))
P_alles<-diag(P_alles)
}else{
P_alles<-matrix(0,lin+m_alt+n*s,lin+m_alt+n*s)
for(jf in 1:n)
P_alles[(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s),(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s)]<-chol2inv(chol(Q1))
}


score_alles<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_alles<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P_alles

  for (r in 1:m)
  {      
  if(r==1)
  {
  Inv_r[[1]]<-chol2inv(chol(F_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[1]]%*%score_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]
  }else{
  Inv_r[[r]]<-chol2inv(chol(F_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[r]]%*%score_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]
  }

  Eta_r[,r]<-Eta+Z_r_gross[[r]]%*%Delta_r[[r]]

  Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))

  M1[[r]]<-(Z_r_gross[[r]]*sqrt(Sigma*D*1/Sigma*D))%*%Inv_r[[r]]%*%t(Z_r_gross[[r]]*sqrt(D*1/Sigma*D*1/Sigma))
  diag_H1[r]<-sum(diag(E-((E-M1[[r]])%*%Uu)))



   ######## likelihood of information criterion    ########
      if (family$family=="poisson")
      IC[l,r]<-sum(y*log(Mu_r[,r])-Mu_r[,r])#-0.5*pena

      if (family$family=="binomial")
      IC[l,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena
      
      if (family$family=="gaussian")
      IC[l,r]<-sum(y*Mu_r[,r]-0.5*(Mu_r[,r]^2))#-0.5*pena

   ############ AIC ##############
      if (control$sel.method=="aic")
      IC[l,r]<-(-IC[l,r])+diag_H1[r]

   ############ BIC ###############
      if (control$sel.method=="bic")
      IC[l,r]<-(-IC[l,r])+0.5*diag_H1[r]*log(n)
   }

   mi_ic<-match(min(IC[l,]),IC[l,])
   
   if(mi_ic==1)   
   {
   g<-c(rep(T,lin),is.element(grp,1:blocksum[1]),rep(T,s*n))
   }else{
   g<-c(rep(T,lin),is.element(grp,(blocksum[mi_ic-1]+1):blocksum[mi_ic]),rep(T,s*n))
   }
   
   Delta[l,]<-Delta[l-1,]

   Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]<-control$nue*Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]
   Delta[l,g]<-Delta[l,g]+t(Delta_r[[mi_ic]])                                                      

komp[l]<-mi_ic

Nue_ma<-c(rep(1,lin),rep(control$nue,block[mi_ic]),rep(1,s*n))

Uu<-(E-((cbind(Z_r[[mi_ic]],W)*sqrt(Sigma*D*1/Sigma*D))%*%(Inv_r[[mi_ic]]*Nue_ma)%*%t(cbind(Z_r[[mi_ic]],W)*sqrt(D*1/Sigma*D*1/Sigma))))%*%Uu
FinalHat[[l]]<-E-Uu
defre<-sum(diag(FinalHat[[l]]))
  
if(mi_ic==1)
{
X1<-cbind(X,as.matrix(U[,1:block[1]]))
}else{
X1<-cbind(X,as.matrix(U[,(blocksum[mi_ic-1]+1):blocksum[mi_ic]]))
}

if(s==1)
{
PP<-rep((Q1^(-1)),n*s)
PP<-diag(PP)
}else{
PP<-matrix(0,n*s,n*s)
for(jf in 1:n)
PP[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Q1))
}

R1<-chol2inv(chol((t(X1)%*%(X1*D*1/Sigma*D))-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP))%*%t(t(X1)%*%(W*D*1/Sigma*D)))))
R2<-(t(X1*D*1/Sigma)-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP)))%*%t(W*D*1/Sigma))

Q_mi<-control$nue*R1[(lin+block[mi_ic]):(lin+block[mi_ic]),]%*%R2%*%(E-FinalHat[[l-1]])
   
    if(is.matrix(Standard_ma[[mi_ic]]))
    {
    Standard_ma[[mi_ic]]<-Standard_ma[[mi_ic]]+Q_mi
    }else{
    Standard_ma[[mi_ic]]<-Q_mi
    }

    for(jh in as.numeric(levels(as.factor(komp[1:l])))[as.numeric(levels(as.factor(komp[1:l])))!=0])
    {
    if(jh==1)
    {
    Standard[l,(lin+1):(lin+blocksum[1])]<-sqrt(diag(Standard_ma[[1]]%*%t(Standard_ma[[1]]*Sigma)))
    }else{
    Standard[l,(lin+blocksum[jh-1]+1):(lin+blocksum[jh])]<-sqrt(diag(Standard_ma[[jh]]%*%t(Standard_ma[[jh]]*Sigma)))
    }}


          
Eta<-Z_alles%*%Delta[l,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-as.vector(family$variance(Mu))

D<-as.vector(family$mu.eta(Eta))



if (control$method=="EM")
{
aktuell_vec<-!is.element(Delta[l,],0)
X_aktuell<-Z_alles[,aktuell_vec] 

m_aktuell<-sum(block[is.element(1:m,komp)])


if(s==1)
{
P_gross<-c(rep(0,lin+m_aktuell),rep((Q1^(-1)),n*s))
P_gross<-diag(P_gross)
}else{
P_gross<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
for(jf in 1:n)
P_gross[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q1))
}


F_gross<-t(X_aktuell)%*%(X_aktuell*D*1/Sigma*D)+P_gross


############################# update Q ################
InvFisher<-chol2inv(chol(F_gross))

Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[l,(lin+m+1):(lin+m+s)]%*%t(Delta[l,(lin+m+1):(lin+m+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)])

#Q1<-InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,seq(from=lin+1,to=lin+1+(s-1)*n,by=n)]%*%t(Delta[1,seq(from=lin+1,to=lin+1+(s-1)*n,by=n)])
#for (i in 2:n)
#Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,seq(from=lin+i,to=lin+i+(s-1)*n,by=n)]%*%t(Delta[1,seq(from=lin+i,to=lin+i+(s-1)*n,by=n)])


Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:(lin+m_alt)]

aktuell_vec<-!is.element(Delta[l,1:(lin+m_alt)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{        
optim.obj<-try(bobyqa(sqrt(Q1),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, model="no.model",  lower = 1e-12, upper = 20))
Q1<-as.matrix(optim.obj$par)^2
}else{
Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, model="no.model", lower=low,upper=upp))

Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])


#### Check for positiv definitness ########
for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
       Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}
}





    if (family$family=="poisson")
    IC_stop[l]<-sum(y*log(Mu)-Mu)

    if (family$family=="binomial")
    IC_stop[l]<-sum(y*log(Mu)+(1-y)*log(1-Mu))

    if (family$family=="gaussian")
    IC_stop[l]<-sum(y*Mu-0.5*(Mu^2))
    
   ############ AIC ##############
    if (control$sel.method=="aic")
    IC_stop[l]<-(-IC_stop[l])+defre

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC_stop[l]<-(-IC_stop[l])+0.5*defre*log(n)


if(control$overdispersion==T)
{
phi[l]<-(sum((y-Mu)^2/Mu))/(N-defre)
Sigma<-Sigma*phi[l]
}

Q[,,l+1]<-Q1

kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l,]^2))


if(kritval<1e-8)
break


kritval2<-abs(IC_stop[l-1]-IC_stop[l])/abs(IC_stop[l])
      
if(kritval2<1e-8)
break


}}
if (l<control$steps)
{
for (wj in (l+1):control$steps)
{
Delta[wj,]<-Delta[l,]
Q[,,wj+1]<-Q1
IC[wj,]<-IC[l,]
IC_stop[wj]<-IC_stop[l]
}}


##########


opt<-match(min(IC_stop),IC_stop)


if (control$OPT==TRUE)
{
Delta_neu<-Delta[opt,]
}else{
Delta_neu<-Delta[l,]
}


Eta_opt<-Z_alles%*%Delta_neu

Mu_opt<-as.vector(family$linkinv(Eta_opt))

Sigma_opt<-as.vector(family$variance(Mu_opt))    

D_opt<-as.vector(family$mu.eta(Eta_opt))

aktuell_vec<-!is.element(Delta_neu,0)
X_aktuell<-Z_alles[,aktuell_vec] 


m_aktuell<-sum(block[is.element(1:m,komp[1:opt])])


if(s==1)
{
P_opt<-c(rep(0,lin+m_aktuell),rep((Q[,,opt+1]^(-1)),n*s))
P_opt<-diag(P_opt)
}else{
P_opt<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
for(jf in 1:n)
P_opt[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q[,,opt+1]))
}



F_opt<-t(X_aktuell)%*%(X_aktuell*D_opt*1/Sigma_opt*D_opt)+P_opt


Inv_F_opt<-chol2inv(chol(F_opt))
Lin_errors<-sqrt(diag(Inv_F_opt))


Standard_errors<-Standard[opt,]
g3<-c(rep(T,lin),rep(F,m_alt),rep(T,s*n))
g4<-c(rep(T,lin),rep(F,m_aktuell),rep(T,s*n))
Standard_errors[g3]<-Lin_errors[g4]

y_dach<-as.vector(family$linkinv(Eta_opt))

if(s==1)
{
Qfinal<-sqrt(Q[,,opt+1])
}else{
Qfinal<-Q[,,opt+1]
}

if(!is.matrix(Qfinal))
Qfinal<-as.matrix(Qfinal)
colnames(Qfinal)<-random.labels
rownames(Qfinal)<-random.labels

names(Delta_neu)[1]<-"(Intercept)"
names(Standard_errors)[1]<-"(Intercept)"
if(lin>1)
{
names(Delta_neu)[2:lin]<-colnames(X)[2:lin]
names(Standard_errors)[2:lin]<-colnames(X)[2:lin]
}
names(Delta_neu)[(lin+1):(lin+m_alt)]<-colnames(U)
names(Standard_errors)[(lin+1):(lin+m_alt)]<-colnames(U)

permut<-numeric()
for(ys in 1:length(old.names))
permut<-c(permut,match(old.names[ys],names(Delta_neu)))

Delta_neu[1:length(old.names)]<-Delta_neu[permut]
Standard_errors[1:length(old.names)]<-Standard_errors[permut]

names(Delta_neu)[1:length(old.names)]<-old.names
names(Standard_errors)[1:length(old.names)]<-old.names


Delta[,1:length(old.names)]<-Delta[,permut]

names(Delta_neu)[(lin+m_alt+1):(lin+m_alt+n*s)]<-colnames(W)
names(Standard_errors)[(lin+m_alt+1):(lin+m_alt+n*s)]<-colnames(W)
colnames(Delta)<-c(old.names,colnames(W))

colnames(IC)<-very.old.names

comp<-character()
for(oi in 1:length(komp))
comp<-c(comp,very.old.names[komp[oi]])

ret.obj=list()
ret.obj$IC<-IC
ret.obj$IC_sel<-IC_stop
ret.obj$opt<-opt
ret.obj$Deltamatrix<-Delta
ret.obj$ranef<-Delta_neu[(lin+m_alt+1):(lin+m_alt+n*s)]
ret.obj$coefficients<-Delta_neu[1:(lin+m_alt)]
ret.obj$fixerror<-Standard_errors[1:(lin+m_alt)]
ret.obj$ranerror<-Standard_errors[(lin+m_alt+1):(lin+m_alt+n*s)]
ret.obj$components<-comp
ret.obj$Q_long<-Q
ret.obj$Q<-Qfinal
ret.obj$y_hat<-Mu_opt
ret.obj$HatMatrix<-FinalHat[[opt]]
ret.obj$phi<-phi[opt]
ret.obj$family<-family
ret.obj$fix<-fix
ret.obj$newrndfrml<-newrndfrml
ret.obj$subject<-names(rnd)[1]
ret.obj$data<-data
ret.obj$factor.names<-factor.names
ret.obj$factor.list<-factor.list
return(ret.obj)
}



bGLMM <- function(fix=formula, rnd=formula, data,  family=NULL, control=list()) UseMethod("bGLMM")

bGLMM.formula <- function(fix,rnd,...)
{
est <- est.bGLMM(fix,rnd,...)
est$fitted.values <- est$y_hat
est$StdDev <- est$Q
est$call <- match.call()
class(est) <- "bGLMM"
est
}


print.bGLMM <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}




summary.bGLMM <- function(object, ...)
{
se <- object$fixerror
zval <- coefficients(object) / se
TAB <- cbind(Estimate = coefficients(object),
StdErr = se,
z.value = zval,
p.value = 2*pnorm(-abs(zval)))
res <- list(call=object$call,
coefficients=TAB,StdDev=object$StdDev)
class(res) <- "summary.bGLMM"
res
}


print.summary.bGLMM <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}


predict.bGLMM <- function(object,newdata=NULL,...)
{
if(is.null(newdata))
{
y<-fitted(object)
}else{                   
family<-object$family

if(length(object$factor.names>0))
{
for (i in 1:length(object$factor.names))
{
newdata[,object$factor.names[i]]<-factor(newdata[,object$factor.names[i]],levels=object$factor.list[[i]])
}}

X <- model.matrix(object$fix, newdata)

subj.new<-levels(as.factor(newdata[,object$subject]))
subj.old<-levels(as.factor(object$data[,object$subject]))
subj.test<-is.element(subj.new,subj.old)
subj.ok<-subj.new[subj.test]

krit.random<-!all(!is.element(subj.new,subj.old))
if(krit.random)
{
W_start <- model.matrix(formula(object$newrndfrml), newdata)
}else{
W_start <- NULL
}

rnlabels<-terms(formula(object$newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
s<-length(random.labels)
k<-table(newdata[,colnames(newdata)==(object$subject)])   
n<-length(k)

if(s>1)
{for (i in 2:s)
subj.test<-cbind(subj.test,subj.test)
subj.test<-as.vector(t(subj.test))
}


if(krit.random)
{
if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}
y<- as.vector(family$linkinv(X%*%object$coef))
rand.ok<-is.element(newdata[,object$subject],subj.ok)
W.neu<-W[,subj.test]
y[rand.ok]<- family$linkinv(cbind(X,W.neu)[rand.ok,]%*%c(object$coef,object$ranef[match(colnames(W.neu),names(object$ranef))]))
}else{
W<-NULL
y<- as.vector(family$linkinv(X%*%object$coef))
}
}
y
}
