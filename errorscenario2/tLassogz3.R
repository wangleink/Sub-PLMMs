# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))

rm(list=ls())
library(MASS)
library(mvtnorm)
library(caret)
library(gbm)
library(randomForest)
library(nleqslv)
library(glmnet)
library(foreach)
library(doParallel)

source('Fun_2.R')

ncore=10

method="lasso"

tree.num=500
K=5
N=1*1e5
p=3;q=3;d=200

r1=600
r2=c(600,800,1000,1200)


beta0=rep(1,p)
alpha0=rep(0.5,q)
total0=c(beta0)+c(2*diag(1,q,p)%*%alpha0)
indirect0=total0-beta0

cc=0

case=3
nsim=500


#####Generate Data
#y,a,x
set.seed(cc)
Data<-GData(case,N,beta0,alpha0,p,q,d)


Y<-Data$Y
X<-Data$X
M<-Data$M
Z<-Data$Z
XMZ<-cbind(X,M,Z)
XZ<-cbind(X,Z)
gfun<-Data$gfun
EX0<-Data$EX0
EM0<-Data$EM0
L0<-Data$L0
hfun<-Data$h0
XM<-cbind(X,M)
EXM0<-cbind(EX0,EM0)



Direct.Beta.Unif=Direct.Se.Unif=Direct.Var.Unif=Direct.Cp.Unif=array(NA,c(length(r2),p,nsim))
Direct.Beta.Opt=Direct.Se.Opt=Direct.Var.Opt=Direct.Cp.Opt=array(NA,c(length(r2),p,nsim))
Direct.Beta.Oracle=Direct.Se.Oracle=Direct.Var.Oracle=Direct.Cp.Oracle=array(NA,c(length(r2),p,nsim))
Direct.Beta.linear=Direct.Se.linear=Direct.Var.linear=Direct.Cp.linear=array(NA,c(length(r2),p,nsim))

Total.Beta.Unif=Total.Se.Unif=Total.Var.Unif=Total.Cp.Unif=array(NA,c(length(r2),p,nsim))
Total.Beta.Opt=Total.Se.Opt=Total.Var.Opt=Total.Cp.Opt=array(NA,c(length(r2),p,nsim))
Total.Beta.Oracle=Total.Se.Oracle=Total.Var.Oracle=Total.Cp.Oracle=array(NA,c(length(r2),p,nsim))
Total.Beta.linear=Total.Se.linear=Total.Var.linear=Total.Cp.linear=array(NA,c(length(r2),p,nsim))

Indirect.Beta.Unif=Indirect.Se.Unif=Indirect.Var.Unif=Indirect.Cp.Unif=array(NA,c(length(r2),p,nsim))
Indirect.Beta.Opt=Indirect.Se.Opt=Indirect.Var.Opt=Indirect.Cp.Opt=array(NA,c(length(r2),p,nsim))
Indirect.Beta.Oracle=Indirect.Se.Oracle=Indirect.Var.Oracle=Indirect.Cp.Oracle=array(NA,c(length(r2),p,nsim))
Indirect.Beta.linear=Indirect.Se.linear=Indirect.Var.linear=Indirect.Cp.linear=array(NA,c(length(r2),p,nsim))




cl <- makeCluster(ncore)
registerDoParallel(cl)
packages <- c("glmnet","mvtnorm","MASS","foreach","doParallel","caret","nleqslv","glmnet","randomForest","ncvreg","geepack","gee")
clusterExport(cl, varlist="packages")
clusterEvalQ(cl, lapply(packages, require, character.only=TRUE))

ALLresult <- foreach(run=1:nsim, .combine=abind::abind, .multicombine=TRUE, .init=array(0,c(nrow=240,ncol=p,ndim=0))) %dopar% {
  
  set.seed(run+cc*nsim)
  cat(run)
  
  ##r0 estimation
  idx.init=c()
  idx.init=sample(1:N,r1,T,prob=rep(1/N,N))
  
  Y.init=Y[idx.init]
  X.init=X[idx.init,]
  M.init=M[idx.init,]
  Z.init=Z[idx.init,]
  XM.init=XM[idx.init,]
  
  L0.init=L0[idx.init]
  EX0.init=EX0[idx.init,]
  EM0.init=EM0[idx.init,]
  EXM0.init=EXM0[idx.init,]
  
  
  ##Pred estimation
  Pred<-TandP(Y.init,X.init,M.init,Z.init,method,Z)
  YZ.pred<-Pred[[1]]
  XZ.pred<-Pred[[2]]
  MZ.pred<-Pred[[3]]
  XMZ.pred<-cbind(XZ.pred,MZ.pred)
  
  #Beta.init estimation
  InitResult=DML(Y.init,X.init,M.init,Z.init,K=K,method)
  
  beta.xm.init<-InitResult[[1]]
  beta.total.init<-InitResult[[2]]
  
  ####### Probabilities of direct
  ###A
  Phi.direct.init.ginv=ginv(t(XM.init-XMZ.pred[idx.init])%*%(XM.init-XMZ.pred[idx.init]))
  D=cbind(diag(1,p,p),matrix(0,nrow = p,ncol = q))
  
  pi2.direct=sqrt((Y-(XM-XMZ.pred)%*%beta.xm.init-YZ.pred)^2*rowSums(t(D%*%Phi.direct.init.ginv%*%t(XM-XMZ.pred))^2))
  pi2.direct<-pmax( pi2.direct,1e-6)
  pi.direct.Opt<-pi2.direct/sum(pi2.direct)
  
  ###U
  pi.direct.Unif<-rep(1/N,N)
  
  
  ####### Probabilities of total
  
  ###A
  Phi.total.init.ginv=ginv(t(X.init-XZ.pred[idx.init])%*%(X.init-XZ.pred[idx.init]))
  
  pi2.total=sqrt((Y-(X-XZ.pred)%*%beta.total.init-YZ.pred)^2*rowSums(t(Phi.total.init.ginv%*%t(X-XZ.pred))^2))
  pi2.total<-pmax( pi2.total,1e-6)
  pi.total.Opt<-pi2.total/sum(pi2.total)
  
  ###U
  pi.total.Unif<-rep(1/N,N)
  
  
  ####### Probabilities of indirect
  
  ###A
  pi2.indirect<-rep(NA,N)
  for (int in 1:N) {
    a_i<-((Y[int]-(XM[int,]-XMZ.pred[int,])%*%beta.xm.init-YZ.pred[int])^2 )*(norm(Phi.total.init.ginv%*%(X[int,]-XZ.pred[int,])-D%*%Phi.direct.init.ginv%*%(XM[int,]-XMZ.pred[int,]),"2")^2)
    b_i<-(((X[int,]-XZ.pred[int,])%*%beta.total.init-(XM[int,]-XMZ.pred[int,])%*%beta.xm.init)^2)*(norm(Phi.total.init.ginv%*%(X[int,]-XZ.pred[int,]),"2")^2)
    pi2.indirect[int]=sqrt(a_i+b_i)
  }
  pi2.indirect<-pmax(pi2.indirect,1e-6)
  pi.indirect.Opt<-pi2.indirect/sum(pi2.indirect)
  

  ###U
  pi.indirect.Unif<-rep(1/N,N)
  
  
  
  
  
  #######################################################Oracle and linear
  
  #Beta.init estimation when functions is known
  
  ###direct
  Y.oracle.init<-Y.init-gfun[idx.init]
  tmp.xm.oracle=coef(lm(Y.oracle.init~XM.init-1))

  ##opt
  Phi.direct.init.ginv=ginv(t(XM.init)%*%(XM.init))
  D=cbind(diag(1,p,p),matrix(0,nrow = p,ncol = q))
  
  Y.oracle<-Y-gfun
  pi2.direct.oracle=sqrt(( Y.oracle-(XM)%*%tmp.xm.oracle)^2*rowSums(t(D%*%Phi.direct.init.ginv%*%t(XM))^2))
  pi2.direct.oracle<-pmax( pi2.direct.oracle,1e-6)
  pi.direct.Opt.oracle<-pi2.direct.oracle/sum(pi2.direct.oracle)
  

  ###total
  Y.total.oracle.init<-Y.init-gfun[idx.init]-(hfun[idx.init,]%*%alpha0)
  tmp.total.oracle=coef(lm(Y.total.oracle.init~X.init-1))
  
  Y.total.oracle<-Y-gfun-(hfun%*%alpha0)
  ##opt
  Phi.total.init.ginv=ginv(t(X.init)%*%(X.init))
  pi2.total.oracle=sqrt((Y.total.oracle-(X)%*%tmp.total.oracle)^2*rowSums(t(Phi.total.init.ginv%*%t(X))^2))
  pi2.total.oracle<-pmax( pi2.total.oracle,1e-6)
  pi.total.Opt.oracle<-pi2.total.oracle/sum(pi2.total.oracle)
  
  ##indirect
  #opt
  pi2.indirect.oracle<-rep(NA,N)
  for (int in 1:N) {
    a_i<-((Y.oracle[int]-(XM[int,])%*%tmp.xm.oracle)^2)*(norm(Phi.total.init.ginv%*%(X[int,])-D%*%Phi.direct.init.ginv%*%(XM[int,]),"2")^2)
    b_i<-(((X[int,])%*%tmp.total.oracle-(XM[int,])%*%tmp.xm.oracle)^2)*(norm(Phi.total.init.ginv%*%(X[int,]),"2")^2)
    pi2.indirect.oracle[int]=sqrt(a_i+b_i)
  }
  pi2.indirect.oracle<-pmax(pi2.indirect.oracle,1e-6)
  pi.indirect.Opt.oracle<-pi2.indirect.oracle/sum(pi2.indirect.oracle)
  
  
  
  
  
  ################################################ linear #########################################
  
  ###direct
  XMZ<-cbind(X,M,Z)
  XMZ.init=cbind(X.init,M.init,Z.init)
  tmp.xm.linear=coef(lm(Y.init~XMZ.init-1))
  
  ##opt
  Phi.direct.init.ginv=ginv(t(XMZ.init)%*%(XMZ.init))
  #D.linear=cbind(diag(1,p,p),matrix(0,nrow = p,ncol = q+d))
  #pi2.direct.linear=sqrt((Y-(XMZ)%*%tmp.xm.linear)^2*rowSums(t(D.linear%*%Phi.direct.init.ginv%*%t(XMZ))^2))
  pi2.direct.linear=sqrt((Y-(XMZ)%*%tmp.xm.linear)^2*rowSums(t(Phi.direct.init.ginv%*%t(XMZ))^2))
  
  pi2.direct.linear<-pmax(pi2.direct.linear,1e-6)
  pi.direct.Opt.linear<-pi2.direct.linear/sum(pi2.direct.linear)
  
  
  ###total
  XZ.init=cbind(X.init,Z.init)
  XZ<-cbind(X,Z)
  tmp.total.linear=coef(lm(Y.init~XZ.init-1))
  
  ##opt
  Phi.total.init.ginv=ginv(t(XZ.init)%*%(XZ.init))
  #D.linear=cbind(diag(1,p,p),matrix(0,nrow = p,ncol =d))
  #pi2.total.linear=sqrt((Y.init-(XZ)%*%tmp.total.linear)^2*rowSums(t(D.linear%*%Phi.total.init.ginv%*%t(XZ))^2))
  pi2.total.linear=sqrt((Y-(XZ)%*%tmp.total.linear)^2*rowSums(t(Phi.total.init.ginv%*%t(XZ))^2))
  pi2.total.linear<-pmax(pi2.total.linear,1e-6)
  pi.total.Opt.linear<-pi2.total.linear/sum(pi2.total.linear)
  
  ##indirect
  pi2.indirect.linear<-rep(NA,N)
  D.linear.direct=cbind(diag(1,p,p),matrix(0,nrow = p,ncol = q+d))
  D.linear.total=cbind(diag(1,p,p),matrix(0,nrow = p,ncol = d))
  for (int in 1:N) {
    a_i<-((Y[int]-(XMZ[int,])%*%tmp.xm.linear)^2)*(norm(D.linear.total%*%Phi.total.init.ginv%*%(XZ[int,])-D.linear.direct%*%Phi.direct.init.ginv%*%(XMZ[int,]),"2")^2)
    b_i<-(((XZ[int,])%*%tmp.total.linear-(XMZ[int,])%*%tmp.xm.linear)^2)*(norm(D.linear.total%*%Phi.total.init.ginv%*%(XZ[int,]),"2")^2)
    pi2.indirect.linear[int]=sqrt(a_i+b_i)
  }
  pi2.indirect.linear<-pmax(pi2.indirect.linear,1e-6)
  pi.indirect.Opt.linear<-pi2.indirect.linear/sum(pi2.indirect.linear)

  
  # pi.direct.Opt.linear<-rep(1/N,N)
  # pi.indirect.Opt.linear<-rep(1/N,N)
  # pi.total.Opt.linear<-rep(1/N,N)
  
  #########################################################
 
  beta.direct.Unif=se.direct.Unif=var.direct.Unif=cp.direct.Unif=c()
  beta.direct.Opt=se.direct.Opt=var.direct.Opt=cp.direct.Opt=c()
  beta.direct.Opt.oracle=se.direct.Opt.oracle=var.direct.Opt.oracle=cp.direct.Opt.oracle=c()
  beta.direct.Opt.linear=se.direct.Opt.linear=var.direct.Opt.linear=cp.direct.Opt.linear=c()
  
  beta.total.Unif=se.total.Unif=var.total.Unif=cp.total.Unif=c()
  beta.total.Opt=se.total.Opt=var.total.Opt=cp.total.Opt=c()
  beta.total.Opt.oracle=se.total.Opt.oracle=var.total.Opt.oracle=cp.total.Opt.oracle=c()
  beta.total.Opt.linear=se.total.Opt.linear=var.total.Opt.linear=cp.total.Opt.linear=c()
  
  beta.indirect.Unif=se.indirect.Unif=var.indirect.Unif=cp.indirect.Unif=c()
  beta.indirect.Opt=se.indirect.Opt=var.indirect.Opt=cp.indirect.Opt=c()
  beta.indirect.Opt.oracle=se.indirect.Opt.oracle=var.indirect.Opt.oracle=cp.indirect.Opt.oracle=c()
  beta.indirect.Opt.linear=se.indirect.Opt.linear=var.indirect.Opt.linear=cp.indirect.Opt.linear=c()
  
  ###################################### Subsampling
  for(rr in r2){
    
  
    ######direct
    
    ##Opt
    direct.Opt_result=sampling_fun_direct(rr,pi.direct.Opt, XMZ.pred, YZ.pred)
    beta.direct.Opt=rbind(beta.direct.Opt,direct.Opt_result$est)
    se.direct.Opt=rbind(se.direct.Opt,direct.Opt_result$se)
    var.direct.Opt=rbind(var.direct.Opt,direct.Opt_result$var)
    cp.direct.Opt=rbind(cp.direct.Opt,direct.Opt_result$cp)
    
    ##Unif
    direct.Unif_result=sampling_fun_direct(rr,pi.direct.Unif, XMZ.pred, YZ.pred)
    beta.direct.Unif=rbind(beta.direct.Unif,direct.Unif_result$est)
    se.direct.Unif=rbind(se.direct.Unif,direct.Unif_result$se)
    var.direct.Unif=rbind(var.direct.Unif,direct.Unif_result$var)
    cp.direct.Unif=rbind(cp.direct.Unif,direct.Unif_result$cp)
    
    
    #Opt.oracle
    idx.sel=PI.sel=c()
    idx.sel=sample(1:N,rr,T,prob=pi.direct.Opt.oracle);PI.sel=pi.direct.Opt.oracle[idx.sel]
    Y.sel=Y[idx.sel];X.sel=X[idx.sel,];XM.sel=XM[idx.sel,];gfun.sel=gfun[idx.sel]
    Y.oracle.sel<-Y.sel-gfun.sel
    
    equ.oracle.direct=function(beta){
      h=apply(c((Y.oracle.sel-XM.sel%*%beta)/PI.sel)*XM.sel,2,mean)
      return(h)
    }
    tmp.xm.oracle=nleqslv(rep(0,p+q),equ.oracle.direct,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
  
    grad.sel=t(XM.sel)%*%diag(c(1/PI.sel))%*%(XM.sel)
    score.sel=c(Y.oracle.sel-(XM.sel)%*%tmp.xm.oracle)*(XM.sel)
    score.sel2=t(score.sel)%*%diag(c(1/(PI.sel)^2+rr/(PI.sel)))%*%score.sel
    tmp.se=sqrt(diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel)))[1:p]
    tmp.var=diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel))[1:p]
    tmp.cp=abs(tmp.xm.oracle[1:p]-beta0)<tmp.se*qnorm(1-0.05/2)
    
    beta.direct.Opt.oracle=rbind(beta.direct.Opt.oracle, tmp.xm.oracle[1:p])
    se.direct.Opt.oracle=rbind(se.direct.Opt.oracle, tmp.se)
    var.direct.Opt.oracle=rbind(var.direct.Opt.oracle, tmp.var)
    cp.direct.Opt.oracle=rbind(cp.direct.Opt.oracle,tmp.cp)
    
    
    
    #Opt.linear
    idx.sel=PI.sel=c()
    idx.sel=sample(1:N,rr,T,prob=pi.direct.Opt.linear);PI.sel=pi.direct.Opt.linear[idx.sel]
    Y.sel=Y[idx.sel];X.sel=X[idx.sel,];XM.sel=XM[idx.sel,];XMZ.sel=XMZ[idx.sel,]
    
    equ.linear.direct=function(beta){
      h=apply(c((Y.sel-XMZ.sel%*%beta)/PI.sel)*XMZ.sel,2,mean)
      return(h)
    }
    tmp.xm.linear=nleqslv(rep(0,p+q+d),equ.linear.direct,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    grad.sel=t(XMZ.sel)%*%diag(c(1/PI.sel))%*%(XMZ.sel)
    score.sel=c(Y.sel-(XMZ.sel)%*%tmp.xm.linear)*(XMZ.sel)
    score.sel2=t(score.sel)%*%diag(c(1/(PI.sel)^2+rr/(PI.sel)))%*%score.sel
    tmp.se=sqrt(diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel)))[1:p]
    tmp.var=diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel))[1:p]
    tmp.cp=abs(tmp.xm.linear[1:p]-beta0)<tmp.se*qnorm(1-0.05/2)
    
    beta.direct.Opt.linear=rbind(beta.direct.Opt.linear, tmp.xm.linear[1:p])
    se.direct.Opt.linear=rbind(se.direct.Opt.linear, tmp.se)
    var.direct.Opt.linear=rbind(var.direct.Opt.linear, tmp.var)
    cp.direct.Opt.linear=rbind(cp.direct.Opt.linear,tmp.cp)
    
    
    
    
    
    
    ######total
    
    #Opt
    total.Opt_result=sampling_fun_total(rr,pi.total.Opt, XZ.pred, YZ.pred)
    beta.total.Opt=rbind(beta.total.Opt,total.Opt_result$est)
    se.total.Opt=rbind(se.total.Opt,total.Opt_result$se)
    var.total.Opt=rbind(var.total.Opt,total.Opt_result$var)
    cp.total.Opt=rbind(cp.total.Opt,total.Opt_result$cp)
    
    #Unif
    total.Unif_result=sampling_fun_total(rr,pi.total.Unif, XZ.pred, YZ.pred)
    beta.total.Unif=rbind(beta.total.Unif,total.Unif_result$est)
    se.total.Unif=rbind(se.total.Unif,total.Unif_result$se)
    var.total.Unif=rbind(var.total.Unif,total.Unif_result$var)
    cp.total.Unif=rbind(cp.total.Unif,total.Unif_result$cp)
    
    #Opt.oracle
    idx.sel=PI.sel=c()
    idx.sel=sample(1:N,rr,T,prob=pi.total.Opt.oracle);PI.sel=pi.total.Opt.oracle[idx.sel]
    Y.sel=Y[idx.sel];X.sel=X[idx.sel,];XM.sel=XM[idx.sel,];gfun.sel=gfun[idx.sel];hfun.sel=hfun[idx.sel,]
    Y.oracle.sel<-Y.sel-gfun.sel-hfun.sel%*%alpha0
    
    equ.oracle.total=function(beta){
      h=apply(c((Y.oracle.sel-X.sel%*%beta)/PI.sel)*X.sel,2,mean)
      return(h)
    }
    tmp.total.oracle=nleqslv(rep(0,p),equ.oracle.total,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    
    grad.sel=t(X.sel)%*%diag(c(1/PI.sel))%*%(X.sel)
    score.sel=c(Y.oracle.sel-(X.sel)%*%tmp.total.oracle)*(X.sel)
    score.sel2=t(score.sel)%*%diag(c(1/(PI.sel)^2+rr/(PI.sel)))%*%score.sel
    tmp.se=sqrt(diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel)))
    tmp.var=diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel))
    tmp.cp=abs(tmp.total.oracle-total0)<tmp.se*qnorm(1-0.05/2)

    beta.total.Opt.oracle=rbind(beta.total.Opt.oracle,tmp.total.oracle)
    se.total.Opt.oracle=rbind(se.total.Opt.oracle, tmp.se)
    var.total.Opt.oracle=rbind(var.total.Opt.oracle,tmp.var)
    cp.total.Opt.oracle=rbind(cp.total.Opt.oracle,tmp.cp)
    
    
    ##Opt.linear
    idx.sel=PI.sel=c()
    idx.sel=sample(1:N,rr,T,prob=pi.total.Opt.linear);PI.sel=pi.total.Opt.linear[idx.sel]
    Y.sel=Y[idx.sel];XZ.sel=XZ[idx.sel,]
    equ.linear.total=function(beta){
      h=apply(c((Y.sel-XZ.sel%*%beta)/PI.sel)*XZ.sel,2,mean)
      return(h)
    }
    tmp.total.linear=nleqslv(rep(0,p+d),equ.linear.total,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    grad.sel=t(XZ.sel)%*%diag(c(1/PI.sel))%*%(XZ.sel)
    score.sel=c(Y.sel-(XZ.sel)%*%tmp.total.linear)*(XZ.sel)
    score.sel2=t(score.sel)%*%diag(c(1/(PI.sel)^2+rr/(PI.sel)))%*%score.sel
    tmp.se=sqrt(diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel)))[1:p]
    tmp.var=diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel))[1:p]
    tmp.cp=abs(tmp.total.linear[1:p]-total0)<tmp.se*qnorm(1-0.05/2)
    
    beta.total.Opt.linear=rbind(beta.total.Opt.linear,tmp.total.linear[1:p])
    se.total.Opt.linear=rbind(se.total.Opt.linear, tmp.se)
    var.total.Opt.linear=rbind(var.total.Opt.linear,tmp.var)
    cp.total.Opt.linear=rbind(cp.total.Opt.linear,tmp.cp)
    
    
    
  
    ######Indirect
    
    #Opt
    indirect.Opt_result=sampling_fun_indirect(rr,pi.indirect.Opt, XZ.pred,XMZ.pred, YZ.pred)
    beta.indirect.Opt=rbind(beta.indirect.Opt,indirect.Opt_result$est)
    se.indirect.Opt=rbind(se.indirect.Opt,indirect.Opt_result$se)
    var.indirect.Opt=rbind(var.indirect.Opt,indirect.Opt_result$var)
    cp.indirect.Opt=rbind(cp.indirect.Opt,indirect.Opt_result$cp)
    
    #Unif
    indirect.Unif_result=sampling_fun_indirect(rr,pi.indirect.Unif, XZ.pred, XMZ.pred,YZ.pred)
    beta.indirect.Unif=rbind(beta.indirect.Unif,indirect.Unif_result$est)
    se.indirect.Unif=rbind(se.indirect.Unif,indirect.Unif_result$se)
    var.indirect.Unif=rbind(var.indirect.Unif,indirect.Unif_result$var)
    cp.indirect.Unif=rbind(cp.indirect.Unif,indirect.Unif_result$cp)
    
    
    ###Opt Oracle
    idx.sel=PI.sel=c()
    idx.sel=sample(1:N,rr,T,prob=pi.indirect.Opt.oracle);PI.sel=pi.indirect.Opt.oracle[idx.sel]
    Y.sel=Y[idx.sel];X.sel=X[idx.sel,];XM.sel=XM[idx.sel,];gfun.sel=gfun[idx.sel];hfun.sel=hfun[idx.sel,]
    Y.total.oracle.sel<-Y.sel-gfun.sel-hfun.sel%*%alpha0
    Y.direct.oracle.sel<-Y.sel-gfun.sel
    
    equ.oracle.total=function(beta){
      h=apply(c((Y.total.oracle.sel-X.sel%*%beta)/PI.sel)*X.sel,2,mean)
      return(h)
    }
    tmp.total.oracle=nleqslv(rep(0,p),equ.oracle.total,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    equ.oracle.direct=function(beta){
      h=apply(c((Y.direct.oracle.sel-XM.sel%*%beta)/PI.sel)*XM.sel,2,mean)
      return(h)
    }
    tmp.xm.oracle=nleqslv(rep(0,p+q),equ.oracle.direct,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    grad.sel.total=t(X.sel)%*%diag(c(1/PI.sel))%*%(X.sel)/(N*rr)
    grad.sel.direct=t(XM.sel)%*%diag(c(1/PI.sel))%*%(XM.sel)/(N*rr)
    
    score.sel.sigma1=score.sel.sigma2=matrix(NA,nrow =rr ,ncol = p)
    for (int in 1:rr) {
      score.sel.sigma1[int,]=c(Y.direct.oracle.sel-(XM.sel)%*%tmp.xm.oracle)[int]*
        (ginv(grad.sel.total)%*%(X.sel)[int,]-D%*%ginv(grad.sel.direct)%*%(XM.sel)[int,])
      score.sel.sigma2[int,]=c((X.sel)%*%tmp.total.oracle-(XM.sel)%*%tmp.xm.oracle)[int]*
        (ginv(grad.sel.total)%*%(X.sel)[int,])
    }
    
    score.sel2.sigma1=score.sel2.sigma2=matrix(0,nrow =p ,ncol = p)
    for (int in 1:rr) {
      score.sel2.sigma1=score.sel2.sigma1+score.sel.sigma1[int,]%*%t(score.sel.sigma1[int,])*(1/(PI.sel[int])^2+rr/PI.sel[int])/(N^2*rr^2)
      score.sel2.sigma2=score.sel2.sigma2+score.sel.sigma2[int,]%*%t(score.sel.sigma2[int,])*(1/(PI.sel[int])^2+rr/PI.sel[int])/(N^2*rr^2)
    }
    Sigma.sel<-score.sel2.sigma1+score.sel2.sigma2
    
    
    tmp.se=sqrt(diag(Sigma.sel))
    tmp.var=diag(Sigma.sel)
    tmp.indirect.oracle=tmp.total.oracle-c(tmp.xm.oracle[1:p])
    tmp.cp=abs(tmp.indirect.oracle-indirect0)<tmp.se*qnorm(1-0.05/2)
    
    beta.indirect.Opt.oracle=rbind(beta.indirect.Opt.oracle, tmp.indirect.oracle)
    se.indirect.Opt.oracle=rbind(se.indirect.Opt.oracle, tmp.se)
    var.indirect.Opt.oracle=rbind(var.indirect.Opt.oracle,tmp.var)
    cp.indirect.Opt.oracle=rbind(cp.indirect.Opt.oracle,tmp.cp)
    
    
    ###################linear
    idx.sel=PI.sel=c()
    idx.sel=sample(1:N,rr,T,prob=pi.indirect.Opt.linear);PI.sel=pi.indirect.Opt.linear[idx.sel]
    Y.sel=Y[idx.sel];X.sel=X[idx.sel,];XMZ.sel=XMZ[idx.sel,];XZ.sel=XZ[idx.sel,]
    
    equ.linear.direct=function(beta){
      h=apply(c((Y.sel-XMZ.sel%*%beta)/PI.sel)*XMZ.sel,2,mean)
      return(h)
    }
    tmp.xm.linear=nleqslv(rep(0,p+q+d),equ.linear.direct,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    equ.linear.total=function(beta){
      h=apply(c((Y.sel-XZ.sel%*%beta)/PI.sel)*XZ.sel,2,mean)
      return(h)
    }
    tmp.total.linear=nleqslv(rep(0,p+d),equ.linear.total,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
    
    grad.sel.total=t(XZ.sel)%*%diag(c(1/PI.sel))%*%(XZ.sel)/(N*rr)
    grad.sel.direct=t(XMZ.sel)%*%diag(c(1/PI.sel))%*%(XMZ.sel)/(N*rr)
    
    score.sel.sigma1=score.sel.sigma2=matrix(NA,nrow =rr ,ncol = p)
    for (int in 1:rr) {
      score.sel.sigma1[int,]=c(Y.sel-(XMZ.sel)%*%tmp.xm.linear)[int]*
        (D.linear.total%*%ginv(grad.sel.total)%*%(XZ.sel)[int,]-D.linear.direct%*%ginv(grad.sel.direct)%*%(XMZ.sel)[int,])
      score.sel.sigma2[int,]=c((XZ.sel)%*%tmp.total.linear-(XMZ.sel)%*%tmp.xm.linear)[int]*
        (D.linear.total%*%ginv(grad.sel.total)%*%(XZ.sel)[int,])
    }
    
    score.sel2.sigma1=score.sel2.sigma2=matrix(0,nrow =p ,ncol = p)
    for (int in 1:rr) {
      score.sel2.sigma1=score.sel2.sigma1+score.sel.sigma1[int,]%*%t(score.sel.sigma1[int,])*(1/(PI.sel[int])^2+rr/PI.sel[int])/(N^2*rr^2)
      score.sel2.sigma2=score.sel2.sigma2+score.sel.sigma2[int,]%*%t(score.sel.sigma2[int,])*(1/(PI.sel[int])^2+rr/PI.sel[int])/(N^2*rr^2)
    }
    Sigma.sel<-score.sel2.sigma1+score.sel2.sigma2
    
    
    tmp.se=sqrt(diag(Sigma.sel))
    tmp.var=diag(Sigma.sel)
    tmp.indirect.linear=tmp.total.linear[1:p]-c(tmp.xm.linear[1:p])
    tmp.cp=abs(tmp.indirect.linear-indirect0)<tmp.se*qnorm(1-0.05/2)
    
    beta.indirect.Opt.linear=rbind(beta.indirect.Opt.linear, tmp.indirect.linear)
    se.indirect.Opt.linear=rbind(se.indirect.Opt.linear, tmp.se)
    var.indirect.Opt.linear=rbind(var.indirect.Opt.linear,tmp.var)
    cp.indirect.Opt.linear=rbind(cp.indirect.Opt.linear,tmp.cp)
    
    
    
  
  }#end R2
  
  AA<-matrix(NA,nrow =240 ,ncol = p)
  AA[1:4,]<-beta.direct.Opt;AA[5:8,]<-se.direct.Opt;AA[9:12,]<-cp.direct.Opt;AA[13:16,]<-var.direct.Opt;AA[17:20,]<-0
  AA[21:24,]<-beta.direct.Unif;AA[25:28,]<-se.direct.Unif;AA[29:32,]<-cp.direct.Unif;AA[33:36,]<-var.direct.Unif;AA[37:40,]<-0
  AA[41:44,]<-beta.direct.Opt.oracle;AA[45:48,]<-se.direct.Opt.oracle;AA[49:52,]<-cp.direct.Opt.oracle;AA[53:56,]<-var.direct.Opt.oracle;AA[57:60,]<-0
  AA[61:64,]<-beta.direct.Opt.linear;AA[65:68,]<-se.direct.Opt.linear;AA[69:72,]<-cp.direct.Opt.linear;AA[73:76,]<-var.direct.Opt.linear;AA[77:80,]<-0
  AA[81:84,]<-beta.total.Opt;AA[85:88,]<-se.total.Opt;AA[89:92,]<-cp.total.Opt;AA[93:96,]<-var.total.Opt;AA[97:100,]<-0
  AA[101:104,]<-beta.total.Unif;AA[105:108,]<-se.total.Unif;AA[109:112,]<-cp.total.Unif;AA[113:116,]<-var.total.Unif;AA[117:120,]<-0
  AA[121:124,]<-beta.total.Opt.oracle;AA[125:128,]<-se.total.Opt.oracle;AA[129:132,]<-cp.total.Opt.oracle;AA[133:136,]<-var.total.Opt.oracle;AA[137:140,]<-0
  AA[141:144,]<-beta.total.Opt.linear;AA[145:148,]<-se.total.Opt.linear;AA[149:152,]<-cp.total.Opt.linear;AA[153:156,]<-var.total.Opt.linear;AA[157:160,]<-0
  AA[161:164,]<-beta.indirect.Opt; AA[165:168,]<-se.indirect.Opt; AA[169:172,]<-cp.indirect.Opt; AA[173:176,]<-var.indirect.Opt;AA[177:180,]<-0
  AA[181:184,]<-beta.indirect.Unif; AA[185:188,]<-se.indirect.Unif; AA[189:192,]<-cp.indirect.Unif; AA[193:196,]<-var.indirect.Unif ;AA[197:200,]<-0
  AA[201:204,]<-beta.indirect.Opt.oracle; AA[205:208,]<-se.indirect.Opt.oracle; AA[209:212,]<-cp.indirect.Opt.oracle; AA[213:216,]<-var.indirect.Opt.oracle ;AA[217:220,]<-0
  AA[221:224,]<-beta.indirect.Opt.linear; AA[225:228,]<-se.indirect.Opt.linear; AA[229:232,]<-cp.indirect.Opt.linear; AA[233:236,]<-var.indirect.Opt.linear ;AA[237:240,]<-0
   
  return(AA)

}#end core

stopCluster(cl)


for(int in 1:nsim)
{
  

  Direct.Beta.Opt[,,int]=ALLresult[1:4,,int]; Direct.Se.Opt[,,int]=ALLresult[5:8,,int]
  Direct.Var.Opt[,,int]=ALLresult[13:16,,int]; Direct.Cp.Opt[,,int]=ALLresult[9:12,,int]
  Direct.Beta.Unif[,,int]=ALLresult[21:24,,int]; Direct.Se.Unif[,,int]=ALLresult[25:28,,int]
  Direct.Var.Unif[,,int]=ALLresult[33:36,,int]; Direct.Cp.Unif[,,int]=ALLresult[29:32,,int]
  Direct.Beta.Oracle[,,int]=ALLresult[41:44,,int]; Direct.Se.Oracle[,,int]=ALLresult[45:48,,int]
  Direct.Var.Oracle[,,int]=ALLresult[53:56,,int]; Direct.Cp.Oracle[,,int]=ALLresult[49:52,,int]
  Direct.Beta.linear[,,int]=ALLresult[61:64,,int]; Direct.Se.linear[,,int]=ALLresult[65:68,,int]
  Direct.Var.linear[,,int]=ALLresult[73:76,,int]; Direct.Cp.linear[,,int]=ALLresult[69:72,,int]
  
  
  Total.Beta.Opt[,,int]=ALLresult[81:84,,int]; Total.Se.Opt[,,int]=ALLresult[85:88,,int]
  Total.Var.Opt[,,int]=ALLresult[93:96,,int]; Total.Cp.Opt[,,int]=ALLresult[89:92,,int]
  Total.Beta.Unif[,,int]=ALLresult[101:104,,int]; Total.Se.Unif[,,int]=ALLresult[105:108,,int]
  Total.Var.Unif[,,int]=ALLresult[113:116,,int]; Total.Cp.Unif[,,int]=ALLresult[109:112,,int]
  Total.Beta.Oracle[,,int]=ALLresult[121:124,,int]; Total.Se.Oracle[,,int]=ALLresult[125:128,,int]
  Total.Var.Oracle[,,int]=ALLresult[133:136,,int]; Total.Cp.Oracle[,,int]=ALLresult[129:132,,int]
  Total.Beta.linear[,,int]=ALLresult[141:144,,int]; Total.Se.linear[,,int]=ALLresult[145:148,,int]
  Total.Var.linear[,,int]=ALLresult[153:156,,int]; Total.Cp.linear[,,int]=ALLresult[149:152,,int]
  
  
  Indirect.Beta.Opt[,,int]=ALLresult[161:164,,int]; Indirect.Se.Opt[,,int]=ALLresult[165:168,,int]
  Indirect.Var.Opt[,,int]=ALLresult[173:176,,int]; Indirect.Cp.Opt[,,int]=ALLresult[169:172,,int]
  Indirect.Beta.Unif[,,int]=ALLresult[181:184,,int]; Indirect.Se.Unif[,,int]=ALLresult[185:188,,int]
  Indirect.Var.Unif[,,int]=ALLresult[193:196,,int]; Indirect.Cp.Unif[,,int]=ALLresult[189:192,,int]
  Indirect.Beta.Oracle[,,int]=ALLresult[201:204,,int]; Indirect.Se.Oracle[,,int]=ALLresult[205:208,,int]
  Indirect.Var.Oracle[,,int]=ALLresult[213:216,,int]; Indirect.Cp.Oracle[,,int]=ALLresult[209:212,,int]
  Indirect.Beta.linear[,,int]=ALLresult[221:224,,int]; Indirect.Se.linear[,,int]=ALLresult[225:228,,int]
  Indirect.Var.linear[,,int]=ALLresult[233:236,,int]; Indirect.Cp.linear[,,int]=ALLresult[229:232,,int]
  
}



run=nsim+1


BETA0=array(rep(matrix(rep(beta0,length(r2)),length(r2),length(beta0),byrow=T),(run-1)),
            c(length(r2),length(beta0),(run-1)))

#####direct

bias.direct.Opt=apply(apply(Direct.Beta.Opt[,,(1:(run-1))]-BETA0,c(1,2),mean),1,mean)
abias.direct.Opt=apply(apply(abs(Direct.Beta.Opt[,,(1:(run-1))]-BETA0),c(1,2),mean),1,mean)
mse.direct.Opt=apply(apply( (Direct.Beta.Opt[,,(1:(run-1))]-BETA0)^2,c(1,2),mean),1,sum)
emse.direct.Opt=apply(apply(Direct.Var.Opt[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.direct.Opt=apply(apply(Direct.Beta.Opt[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.direct.Opt=apply(apply(Direct.Se.Opt[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.direct.Opt=2*1.96*ese.direct.Opt
covp.direct.Opt=apply(apply(Direct.Cp.Opt[,,(1:(run-1))],c(1,2),mean),1,mean)


bias.direct.Unif=apply(apply(Direct.Beta.Unif[,,(1:(run-1))]-BETA0,c(1,2),mean),1,mean)
abias.direct.Unif=apply(apply(abs(Direct.Beta.Unif[,,(1:(run-1))]-BETA0),c(1,2),mean),1,mean)
mse.direct.Unif=apply(apply( (Direct.Beta.Unif[,,(1:(run-1))]-BETA0)^2,c(1,2),mean),1,sum)
emse.direct.Unif=apply(apply(Direct.Var.Unif[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.direct.Unif=apply(apply(Direct.Beta.Unif[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.direct.Unif=apply(apply(Direct.Se.Unif[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.direct.Unif=2*1.96*ese.direct.Unif
covp.direct.Unif=apply(apply(Direct.Cp.Unif[,,(1:(run-1))],c(1,2),mean),1,mean)

bias.direct.Opt.oracle=apply(apply(Direct.Beta.Oracle[,,(1:(run-1))]-BETA0,c(1,2),mean),1,mean)
abias.direct.Opt.oracle=apply(apply(abs(Direct.Beta.Oracle[,,(1:(run-1))]-BETA0),c(1,2),mean),1,mean)
mse.direct.Opt.oracle=apply(apply( (Direct.Beta.Oracle[,,(1:(run-1))]-BETA0)^2,c(1,2),mean),1,sum)
emse.direct.Opt.oracle=apply(apply(Direct.Var.Oracle[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.direct.Opt.oracle=apply(apply(Direct.Beta.Oracle[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.direct.Opt.oracle=apply(apply(Direct.Se.Oracle[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.direct.Opt.oracle=2*1.96*ese.direct.Opt.oracle
covp.direct.Opt.oracle=apply(apply(Direct.Cp.Oracle[,,(1:(run-1))],c(1,2),mean),1,mean)


bias.direct.Opt.linear=apply(apply(Direct.Beta.linear[,,(1:(run-1))]-BETA0,c(1,2),mean),1,mean)
abias.direct.Opt.linear=apply(apply(abs(Direct.Beta.linear[,,(1:(run-1))]-BETA0),c(1,2),mean),1,mean)
mse.direct.Opt.linear=apply(apply( (Direct.Beta.linear[,,(1:(run-1))]-BETA0)^2,c(1,2),mean),1,sum)
emse.direct.Opt.linear=apply(apply(Direct.Var.linear[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.direct.Opt.linear=apply(apply(Direct.Beta.linear[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.direct.Opt.linear=apply(apply(Direct.Se.linear[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.direct.Opt.linear=2*1.96*ese.direct.Opt.linear
covp.direct.Opt.linear=apply(apply(Direct.Cp.linear[,,(1:(run-1))],c(1,2),mean),1,mean)



###total


GAMMA0=array(rep(matrix(rep(total0,length(r2)),length(r2),length(beta0),byrow=T),(run-1)),
            c(length(r2),length(beta0),(run-1)))

bias.total.Opt=apply(apply(Total.Beta.Opt[,,(1:(run-1))]-GAMMA0,c(1,2),mean),1,mean)
abias.total.Opt=apply(apply(abs(Total.Beta.Opt[,,(1:(run-1))]-GAMMA0),c(1,2),mean),1,mean)
mse.total.Opt=apply(apply( (Total.Beta.Opt[,,(1:(run-1))]-GAMMA0)^2,c(1,2),mean),1,sum)
emse.total.Opt=apply(apply(Total.Var.Opt[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.total.Opt=apply(apply(Total.Beta.Opt[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.total.Opt=apply(apply(Total.Se.Opt[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.total.Opt=2*1.96*ese.total.Opt
covp.total.Opt=apply(apply(Total.Cp.Opt[,,(1:(run-1))],c(1,2),mean),1,mean)

bias.total.Unif=apply(apply(Total.Beta.Unif[,,(1:(run-1))]-GAMMA0,c(1,2),mean),1,mean)
abias.total.Unif=apply(apply(abs(Total.Beta.Unif[,,(1:(run-1))]-GAMMA0),c(1,2),mean),1,mean)
mse.total.Unif=apply(apply( (Total.Beta.Unif[,,(1:(run-1))]-GAMMA0)^2,c(1,2),mean),1,sum)
emse.total.Unif=apply(apply(Total.Var.Unif[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.total.Unif=apply(apply(Total.Beta.Unif[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.total.Unif=apply(apply(Total.Se.Unif[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.total.Unif=2*1.96*ese.total.Unif
covp.total.Unif=apply(apply(Total.Cp.Unif[,,(1:(run-1))],c(1,2),mean),1,mean)


bias.total.Opt.oracle=apply(apply(Total.Beta.Oracle[,,(1:(run-1))]-GAMMA0,c(1,2),mean),1,mean)
abias.total.Opt.oracle=apply(apply(abs(Total.Beta.Oracle[,,(1:(run-1))]-GAMMA0),c(1,2),mean),1,mean)
mse.total.Opt.oracle=apply(apply( (Total.Beta.Oracle[,,(1:(run-1))]-GAMMA0)^2,c(1,2),mean),1,sum)
emse.total.Opt.oracle=apply(apply(Total.Var.Oracle[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.total.Opt.oracle=apply(apply(Total.Beta.Oracle[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.total.Opt.oracle=apply(apply(Total.Se.Oracle[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.total.Opt.oracle=2*1.96*ese.total.Opt.oracle
covp.total.Opt.oracle=apply(apply(Total.Cp.Oracle[,,(1:(run-1))],c(1,2),mean),1,mean)

bias.total.Opt.linear=apply(apply(Total.Beta.linear[,,(1:(run-1))]-GAMMA0,c(1,2),mean),1,mean)
abias.total.Opt.linear=apply(apply(abs(Total.Beta.linear[,,(1:(run-1))]-GAMMA0),c(1,2),mean),1,mean)
mse.total.Opt.linear=apply(apply( (Total.Beta.linear[,,(1:(run-1))]-GAMMA0)^2,c(1,2),mean),1,sum)
emse.total.Opt.linear=apply(apply(Total.Var.linear[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.total.Opt.linear=apply(apply(Total.Beta.linear[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.total.Opt.linear=apply(apply(Total.Se.linear[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.total.Opt.linear=2*1.96*ese.total.Opt.linear
covp.total.Opt.linear=apply(apply(Total.Cp.linear[,,(1:(run-1))],c(1,2),mean),1,mean)



Beta.Indirect0=array(rep(matrix(rep(indirect0,length(r2)),length(r2),length(beta0),byrow=T),(run-1)),
             c(length(r2),length(beta0),(run-1)))

bias.indirect.Opt=apply(apply(Indirect.Beta.Opt[,,(1:(run-1))]-Beta.Indirect0,c(1,2),mean),1,mean)
abias.indirect.Opt=apply(apply(abs(Indirect.Beta.Opt[,,(1:(run-1))]-Beta.Indirect0),c(1,2),mean),1,mean)
mse.indirect.Opt=apply(apply( (Indirect.Beta.Opt[,,(1:(run-1))]-Beta.Indirect0)^2,c(1,2),mean),1,sum)
emse.indirect.Opt=apply(apply(Indirect.Var.Opt[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.indirect.Opt=apply(apply(Indirect.Beta.Opt[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.indirect.Opt=apply(apply(Indirect.Se.Opt[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.indirect.Opt=2*1.96*ese.indirect.Opt
covp.indirect.Opt=apply(apply(Indirect.Cp.Opt[,,(1:(run-1))],c(1,2),mean),1,mean)


bias.indirect.Unif=apply(apply(Indirect.Beta.Unif[,,(1:(run-1))]-Beta.Indirect0,c(1,2),mean),1,mean)
abias.indirect.Unif=apply(apply(abs(Indirect.Beta.Unif[,,(1:(run-1))]-Beta.Indirect0),c(1,2),mean),1,mean)
mse.indirect.Unif=apply(apply( (Indirect.Beta.Unif[,,(1:(run-1))]-Beta.Indirect0)^2,c(1,2),mean),1,sum)
emse.indirect.Unif=apply(apply(Indirect.Var.Unif[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.indirect.Unif=apply(apply(Indirect.Beta.Unif[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.indirect.Unif=apply(apply(Indirect.Se.Unif[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.indirect.Unif=2*1.96*ese.indirect.Unif
covp.indirect.Unif=apply(apply(Indirect.Cp.Unif[,,(1:(run-1))],c(1,2),mean),1,mean)

bias.indirect.Opt.oracle=apply(apply(Indirect.Beta.Oracle[,,(1:(run-1))]-Beta.Indirect0,c(1,2),mean),1,mean)
abias.indirect.Opt.oracle=apply(apply(abs(Indirect.Beta.Oracle[,,(1:(run-1))]-Beta.Indirect0),c(1,2),mean),1,mean)
mse.indirect.Opt.oracle=apply(apply( (Indirect.Beta.Oracle[,,(1:(run-1))]-Beta.Indirect0)^2,c(1,2),mean),1,sum)
emse.indirect.Opt.oracle=apply(apply(Indirect.Var.Oracle[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.indirect.Opt.oracle=apply(apply(Indirect.Beta.Oracle[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.indirect.Opt.oracle=apply(apply(Indirect.Se.Oracle[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.indirect.Opt.oracle=2*1.96*ese.indirect.Opt.oracle
covp.indirect.Opt.oracle=apply(apply(Indirect.Cp.Oracle[,,(1:(run-1))],c(1,2),mean),1,mean)


bias.indirect.Opt.linear=apply(apply(Indirect.Beta.linear[,,(1:(run-1))]-Beta.Indirect0,c(1,2),mean),1,mean)
abias.indirect.Opt.linear=apply(apply(abs(Indirect.Beta.linear[,,(1:(run-1))]-Beta.Indirect0),c(1,2),mean),1,mean)
mse.indirect.Opt.linear=apply(apply( (Indirect.Beta.linear[,,(1:(run-1))]-Beta.Indirect0)^2,c(1,2),mean),1,sum)
emse.indirect.Opt.linear=apply(apply(Indirect.Var.linear[,,(1:(run-1))],c(1,2),mean),1,sum)
sd.indirect.Opt.linear=apply(apply(Indirect.Beta.linear[,,(1:(run-1))],c(1,2),sd),1,mean)
ese.indirect.Opt.linear=apply(apply(Indirect.Se.linear[,,(1:(run-1))],c(1,2),mean),1,mean)
AL.indirect.Opt.linear=2*1.96*ese.indirect.Opt.linear
covp.indirect.Opt.linear=apply(apply(Indirect.Cp.linear[,,(1:(run-1))],c(1,2),mean),1,mean)







Direct.Opt=rbind(bias.direct.Opt,mse.direct.Opt,emse.direct.Opt,sd.direct.Opt,ese.direct.Opt,covp.direct.Opt)
Direct.Unif=rbind(bias.direct.Unif,mse.direct.Unif,emse.direct.Unif,sd.direct.Unif,ese.direct.Unif,covp.direct.Unif)
Direct.Opt.oracle=rbind(bias.direct.Opt.oracle,mse.direct.Opt.oracle,emse.direct.Opt.oracle,sd.direct.Opt.oracle,ese.direct.Opt.oracle,covp.direct.Opt.oracle)
Direct.Opt.linear=rbind(bias.direct.Opt.linear,mse.direct.Opt.linear,emse.direct.Opt.linear,sd.direct.Opt.linear,ese.direct.Opt.linear,covp.direct.Opt.linear)



Total.Opt=rbind(bias.total.Opt,mse.total.Opt,emse.total.Opt,sd.total.Opt,ese.total.Opt,covp.total.Opt)
Total.Unif=rbind(bias.total.Unif,mse.total.Unif,emse.total.Unif,sd.total.Unif,ese.total.Unif,covp.total.Unif)
Total.Opt.oracle=rbind(bias.total.Opt.oracle ,mse.total.Opt.oracle,emse.total.Opt.oracle,sd.total.Opt.oracle,ese.total.Opt.oracle,covp.total.Opt.oracle)
Total.Opt.linear=rbind(bias.total.Opt.linear,mse.total.Opt.linear,emse.total.Opt.linear,sd.total.Opt.linear,ese.total.Opt.linear,covp.total.Opt.linear)



Indirect.Opt=rbind(bias.indirect.Opt,mse.indirect.Opt,emse.indirect.Opt,sd.indirect.Opt,ese.indirect.Opt,covp.indirect.Opt)
Indirect.Unif=rbind(bias.indirect.Unif,mse.indirect.Unif,emse.indirect.Unif,sd.indirect.Unif,ese.indirect.Unif,covp.indirect.Unif)
Indirect.Opt.oracle=rbind(bias.indirect.Opt.oracle,mse.indirect.Opt.oracle,emse.indirect.Opt.oracle,sd.indirect.Opt.oracle,ese.indirect.Opt.oracle,covp.indirect.Opt.oracle)
Indirect.Opt.linear=rbind(bias.indirect.Opt.linear,mse.indirect.Opt.linear,emse.indirect.Opt.linear,sd.indirect.Opt.linear,ese.indirect.Opt.linear,covp.indirect.Opt.linear)



# result=as.matrix(rbind(cbind(Direct.Opt.oracle[,1],Direct.Opt[,1],Direct.Unif[,1],Direct.Opt.linear[,1],Total.Opt.oracle[,1],Total.Opt[,1],Total.Unif[,1],Total.Opt.linear[,1],Indirect.Opt.oracle[,1],Indirect.Opt[,1],Indirect.Unif[,1],Indirect.Opt.linear[,1]),
#                        cbind(Direct.Opt.oracle[,2],Direct.Opt[,2],Direct.Unif[,2],Direct.Opt.linear[,2],Total.Opt.oracle[,2],Total.Opt[,1],Total.Unif[,2],Total.Opt.linear[,2],Indirect.Opt.oracle[,2],Indirect.Opt[,2],Indirect.Unif[,2],Indirect.Opt.linear[,2]),
#                        cbind(Direct.Opt.oracle[,3],Direct.Opt[,3],Direct.Unif[,3],Direct.Opt.linear[,3],Total.Opt.oracle[,3],Total.Opt[,1],Total.Unif[,3],Total.Opt.linear[,3],Indirect.Opt.oracle[,3],Indirect.Opt[,3],Indirect.Unif[,3],Indirect.Opt.linear[,3]),
#                        cbind(Direct.Opt.oracle[,4],Direct.Opt[,4],Direct.Unif[,4],Direct.Opt.linear[,4],Total.Opt.oracle[,4],Total.Opt[,1],Total.Unif[,4],Total.Opt.linear[,4],Indirect.Opt.oracle[,4],Indirect.Opt[,4],Indirect.Unif[,4],Indirect.Opt.linear[,4])))
# 
# colnames(result)=c("Direcrt:Opt.oracle","Opt","Unif","linear","Total:Opt.oracle","Opt","Unif","linear","Indirect:Oracle","Opt","Unif","linear")
# rownames(result)=c("r400,bias","mse","emse","sd","se","CP",
#                    "r500,bias","mse","emse","sd","se","CP",
#                    "r600,bias","mse","emse","sd","se","CP",
#                    "r700,bias","mse","emse","sd","se","CP")
# 
# 
# result


result1=as.matrix(rbind(cbind(Direct.Opt.oracle[,1],Direct.Opt.linear[,1],Direct.Opt[,1],Direct.Unif[,1]),
                       cbind(Direct.Opt.oracle[,2],Direct.Opt.linear[,2],Direct.Opt[,2],Direct.Unif[,2]),
                       cbind(Direct.Opt.oracle[,3],Direct.Opt.linear[,3],Direct.Opt[,3],Direct.Unif[,3]),
                       cbind(Direct.Opt.oracle[,4],Direct.Opt.linear[,4],Direct.Opt[,4],Direct.Unif[,4])))

colnames(result1)=c("Direcrt:Opt.oracle","linear","Opt","Unif")
rownames(result1)=c("r400,bias","mse","emse","sd","se","CP",
                   "r500,bias","mse","emse","sd","se","CP",
                   "r600,bias","mse","emse","sd","se","CP",
                   "r700,bias","mse","emse","sd","se","CP")

result2=as.matrix(rbind(cbind(Total.Opt.oracle[,1],Total.Opt.linear[,1],Total.Opt[,1],Total.Unif[,1]),
                        cbind(Total.Opt.oracle[,2],Total.Opt.linear[,2],Total.Opt[,2],Total.Unif[,2]),
                        cbind(Total.Opt.oracle[,3],Total.Opt.linear[,3],Total.Opt[,3],Total.Unif[,3]),
                        cbind(Total.Opt.oracle[,4],Total.Opt.linear[,4],Total.Opt[,4],Total.Unif[,4])))

colnames(result2)=c("Total:Opt.oracle","linear","Opt","Unif")
rownames(result2)=c("r400,bias","mse","emse","sd","se","CP",
                    "r500,bias","mse","emse","sd","se","CP",
                    "r600,bias","mse","emse","sd","se","CP",
                    "r700,bias","mse","emse","sd","se","CP")

result3=as.matrix(rbind(cbind(Indirect.Opt.oracle[,1],Indirect.Opt.linear[,1],Indirect.Opt[,1],Indirect.Unif[,1]),
                        cbind(Indirect.Opt.oracle[,2],Indirect.Opt.linear[,2],Indirect.Opt[,2],Indirect.Unif[,2]),
                        cbind(Indirect.Opt.oracle[,3],Indirect.Opt.linear[,3],Indirect.Opt[,3],Indirect.Unif[,3]),
                        cbind(Indirect.Opt.oracle[,4],Indirect.Opt.linear[,4],Indirect.Opt[,4],Indirect.Unif[,4])))

colnames(result3)=c("Indirect:Oracle","linear","Opt","Unif")
rownames(result3)=c("r400,bias","mse","emse","sd","se","CP",
                    "r500,bias","mse","emse","sd","se","CP",
                    "r600,bias","mse","emse","sd","se","CP",
                    "r700,bias","mse","emse","sd","se","CP")


result1
result2
result3



save.image("Lassogz3terror.RData")



