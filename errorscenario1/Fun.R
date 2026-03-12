GData<-function(case,N,beta0,alpha0,p,q,d){

  X<-matrix(NA,nrow =N,ncol = p )
  M<-matrix(NA,nrow = N,ncol = q)
  Z<-matrix(NA,nrow = N,ncol = d)
  Y<-rep(NA,N)

  EX0=matrix(NA,nrow =N,ncol = p )
  EM0=matrix(NA,nrow =N,ncol = q )
  L0=rep(NA,N)

  sigmaz<-matrix(NA,d,d)
  for (j in 1:d) {
    for (jj in 1:d) {
      sigmaz[j,jj]=0.5^{abs(j-jj)}
    }
  }

  Z=rmvnorm(N,mean=rep(0,d),sigma=sigmaz)

  sigma1<-matrix(NA,p,p)
  for (j in 1:p) {
    for (jj in 1:p) {
      sigma1[j,jj]=0.5^{abs(j-jj)}
    }
  }


  # errorx<-rmvnorm(N,mean=rep(0,p),sigma=sigma1/sqrt(2))
  # X[,1]=sin(Z[,1]/2)+errorx[,1]
  # X[,2]=(Z[,2])/2+errorx[,2]
  # X[,3]=0.1*exp((Z[,3])/2)+errorx[,3]
  # 
  # EX0[,1]=sin(Z[,1]/2)
  # EX0[,2]=(Z[,2])/2
  # EX0[,3]=0.1*exp((Z[,3])/2)
  
  errorx<-rmvnorm(N,mean=rep(0,p),sigma=sigma1/2)
  
  X[,1]=0.1*Z[,1]*Z[,2]+errorx[,1]
  X[,2]=0.1*exp((Z[,3]))+errorx[,2]
  X[,3]=log(1+abs(Z[,4]+Z[,5]))+errorx[,3]

   EX0[,1]=0.1*Z[,1]*Z[,2]
   EX0[,2]=0.1*exp((Z[,3]))
   EX0[,3]=log(1+abs(Z[,4]+Z[,5]))
  
  sigma2<-matrix(NA,q,q)
  for (j in 1:q) {
    for (jj in 1:q) {
      sigma2[j,jj]=0.5^{abs(j-jj)}
    }
  }

  errm= rmvnorm(N,mean=rep(0,q),sigma=sigma2/sqrt(2))
  Gamma=2*diag(1,q,p)
  xeta0=rbind(diag(1,q),matrix(0,nrow = d-q,ncol = q))
  M=X%*%Gamma+cos(Z%*%xeta0)+errm

  EM0=EX0%*%Gamma+cos(Z%*%xeta0)
  h0=cos(Z%*%xeta0)

  s2=10
  gamma0=rep(0,d)
  for(j in 1:s2){gamma0[j]=0.4*(1+j/(2*s2))}

  switch(case,
         {
           gfun=Z%*%gamma0
         },{
           gfun=2*exp(Z%*%gamma0)/(1+exp(Z%*%gamma0))
         },{
           #gfun=0.1*Z[,1]*Z[,2]+0.1*Z[,3]*Z[,4]+0.1*Z[,5]^2-0.5*(sin(Z[,6]))^2+0.5*cos(Z[,7])+1/(1+Z[,8]^2)-1/(1+exp(Z[,9]))+0.25*(Z[,10]>0)
           gfun=sin(Z%*%gamma0)+0.5*log(abs(Z%*%gamma0)/(1+abs(Z%*%gamma0)))
         }
         )

  Y=X%*%beta0+M%*%alpha0+gfun+rnorm(N,0,4)

  L0=EX0%*%beta0+EM0%*%alpha0+gfun


  return(list('Y'=Y,'X'=X,'M'=M,'Z'=Z,'gfun'=gfun,'EX0'=EX0,'EM0'=EM0,'L0'=L0,'h0'=h0))
}






L=function(R,C,givenEstimator,Ctest){
  
  Data=cbind(R,C)
  colnames(Data)[1]='R'
  
  tt=length(table(R))
  if(tt==2){
    Data$R=gsub('0','N',Data$R)
    Data$R=gsub('1','P',Data$R)
    Data$R=as.factor(Data$R)
  }
  
  fitControl=trainControl(method="cv",number=5,verboseIter=FALSE)
  
  if((tt==2)&(givenEstimator %in% c('svmLinear2'))){
    Fit=train(R~.,data=Data,method=givenEstimator,trControl=fitControl,
              verbose=FALSE,probability=TRUE)
  }else if(givenEstimator=='nnet'){
    Fit=train(R~.,data=Data,method="nnet",
              #trControl=fitControl,
              linout = TRUE,
              preProcess = c('center', 'scale'),
              maxit = 500,
              tuneGrid = expand.grid(size = 2, decay = 0),
              trControl = trainControl(method = "none", seeds = 1),
              verbose=FALSE) 
  }else{
    Fit=train(R~.,data=Data,method=givenEstimator,trControl=fitControl,
              verbose=FALSE) 
  }
  
  if(tt==2){
    Rtest=predict(Fit,Ctest,type='prob')[,2]
  }else{
    Rtest=predict(Fit,Ctest)
  }
  
  return(Rtest)
}




##############Given train Data and Predict test data
TandP<-function(Y.T,X.T,M.T,Z.T,method,Z.P){
  
  xz.pred<-matrix(NA,nrow=dim(Z.P)[1] ,ncol = p)
  mz.pred<-matrix(NA,nrow=dim(Z.P)[1] ,ncol = p)
  yz.pred<-rep(NA,nrow=dim(Z.P)[1])
  
  if(method=="lasso"){
    for(j in 1:p){
      modelX.Z<-cv.glmnet(Z.T,X.T[,j],family="gaussian",type.measure="mse",nfolds=5)
      xz.pred[,j]=predict(modelX.Z,Z.P,s="lambda.min")
    }
    for(j in 1:q){
      modelM.Z<-cv.glmnet(Z.T,M.T[,j],family="gaussian",type.measure="mse",nfolds=5)
      mz.pred[,j]=predict(modelM.Z,Z.P,s="lambda.min")
    }
    modelY.Z=cv.glmnet(Z.T,Y.T,family="gaussian",type.measure="mse",nfolds=5)
      yz.pred=predict(modelY.Z,Z.P,s="lambda.min")
  }
  if(method=="gbm"){
    for(j in 1:p){
      xz.pred[,j]=L(X.T[,j],as.data.frame(Z.T),givenEstimator=method,Z.P)
    } 
    for(j in 1:q){
      mz.pred[,j]=L(M.T[,j],as.data.frame(Z.T),givenEstimator=method,Z.P)
    }
    yz.pred=L(Y.T,as.data.frame(Z.T),givenEstimator=method,Z.P)
  }
  if(method=="rf"){
    for(j in 1:p){
      data.T=data.frame(y=X.T[,j],Z.T)
      modelX.Z=randomForest(y~.,data=data.T,ntree=tree.num)
      xz.pred[,j]=predict(modelX.Z,data.frame(Z.P))
    }
    for(j in 1:q){
      data.T=data.frame(y=M.T[,j],Z.T)
      modelM.Z=randomForest(y~.,data=data.T,ntree=tree.num)
      mz.pred[,j]=predict(modelM.Z,data.frame(Z.P))
    }
    data.T=data.frame(y=Y.T,Z.T)
    modelY.Z=randomForest(y~.,data=data.T,ntree=tree.num)
    yz.pred=predict(modelY.Z,data.frame(Z.P))
  }
  
  return(list('Y.pre'= yz.pred,'X.pre'=  xz.pred,'M.pre'= mz.pred))
}





Split<-function(K,input){
  INT<-c(1:input)
  Fold<-createFolds(INT,K)
  I1 = rep(0,input)
  for(i in 1:K){
    I1[Fold[[i]]]<-i
  }
  return(I1)
}

DML=function(Y.T,X.T,M.T,Z.T,K,method){
  
  n=length(Y.T)
  I2=Split(K,n)
  
  xz.pred<-matrix(NA,nrow=n,ncol = p)
  mz.pred<-matrix(NA,nrow=n,ncol = p)
  yz.pred<-rep(NA,n)
  
  
  for(kk in 1:K){
    idj<-which(I2==kk); idNj<-which(I2!=kk)
     
    X.train<-X.T[idNj,];M.train<-M.T[idNj,];Y.train<-Y[idNj]; Z.train<-Z.T[idNj,];
    Z.test<-Z.T[idj,]
    
    if(method=="lasso"){
      for(j in 1:p){
        modelX.Z<-cv.glmnet(Z.train,X.train[,j],family="gaussian",type.measure="mse",nfolds=5)
        xz.pred[idj,j]=predict(modelX.Z,Z.test,s="lambda.min")
      }
      for(j in 1:q){
        modelM.Z<-cv.glmnet(Z.train,M.train[,j],family="gaussian",type.measure="mse",nfolds=5)
        mz.pred[idj,j]=predict(modelM.Z,Z.test,s="lambda.min")
      }
      
      modelY.Z=cv.glmnet(Z.train,Y.train,family="gaussian",type.measure="mse",nfolds=5)
      yz.pred[idj]=predict(modelY.Z,Z.test,s="lambda.min")
      
    }else if(method=="gbm"){
      for(j in 1:p){
        xz.pred[idj,j]=L(X.train[,j],as.data.frame(Z.train),givenEstimator=method,Z.test)
      } 
      for(j in 1:q){
        mz.pred[idj,j]=L(M.train[,j],as.data.frame(Z.train),givenEstimator=method,Z.test)
      }
      yz.pred=L(Y.train,as.data.frame(Z.train),givenEstimator=method,Z.test)
      
    }else if(method=="rf"){
      for(j in 1:p){
        data.T=data.frame(y=X.train[,j],Z.train)
        modelX.Z=randomForest(y~.,data=data.T,ntree=tree.num)
        xz.pred[idj,j]=predict(modelX.Z,data.frame(Z.test))
      }
      for(j in 1:q){
        data.T=data.frame(y=M.train[,j],Z.train)
        modelM.Z=randomForest(y~.,data=data.T,ntree=tree.num)
        mz.pred[idj,j]=predict(modelM.Z,data.frame(Z.test))
      }
      data.T=data.frame(y=Y.train,Z.train)
      modelY.Z=randomForest(y~.,data=data.T,ntree=tree.num)
      yz.pred=predict(modelY.Z,data.frame(Z.test))
    }#end ML algorithm

    }#end data-splitting K
  
   XM.T<-cbind(X.T,M.T)
   xmz.pred<-cbind(xz.pred,mz.pred)
  
   score_direct=function(beta){
    h=apply(c((Y.T-(XM.T-xmz.pred)%*%beta-yz.pred))*(XM.T-xmz.pred),2,mean)
    return(h)}
    
   tmp.direct.xm=nleqslv(rep(0,p+q),score_direct,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
  
   score_total=function(beta){
     h=apply(c((Y.T-(X.T-xz.pred)%*%beta-yz.pred))*(X.T-xz.pred),2,mean)
    return(h)}
   
   tmp.total=nleqslv(rep(0,p),score_total,control=list(xtol=1e-8,ftol=1e-8,btol=1e-8))$x
  
   return(list('beta.xm'=tmp.direct.xm,'beta.total'=tmp.total))
}



sampling_fun_direct=function(r.sel,pi.sel,pred.XM,pred.Y){
  idx.sel=PI.sel=c()
  idx.sel=sample(1:N,r.sel,T,prob=pi.sel);PI.sel=pi.sel[idx.sel]
  Y.sel=Y[idx.sel]
  XM.sel=XM[idx.sel,]
  Pred.XM.sel=pred.XM[idx.sel,]
  Pred.Y.sel=pred.Y[idx.sel] 
  
  score.sel=function(beta){
    h=apply(c((Y.sel-(XM.sel-Pred.XM.sel)%*%beta-Pred.Y.sel)/PI.sel)*(XM.sel-Pred.XM.sel),2,mean)
    return(h)
  }
  
  beta.sel=nleqslv(rep(0,p+q),score.sel,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  
  
  grad.sel=t(XM.sel-Pred.XM.sel)%*%diag(c(1/PI.sel))%*%(XM.sel-Pred.XM.sel)
  score.sel=c(Y.sel-(XM.sel-Pred.XM.sel)%*%beta.sel-Pred.Y.sel)*(XM.sel-Pred.XM.sel)
  score.sel2=t(score.sel)%*%diag(c(1/(PI.sel)^2+r.sel/(PI.sel)))%*%score.sel
  
  
  tmp.se=sqrt(diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel)))[1:p]
  tmp.var=diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel))[1:p]
  tmp.cp=abs(beta.sel[1:p]-beta0)<tmp.se*qnorm(1-0.05/2)
  return(list(est=beta.sel[1:p],se=tmp.se,var=tmp.var,cp=tmp.cp))
}




sampling_fun_total=function(r.sel,pi.sel,pred.X,pred.Y){
  idx.sel=PI.sel=c()
  idx.sel=sample(1:N,r.sel,T,prob=pi.sel);PI.sel=pi.sel[idx.sel]
  Y.sel=Y[idx.sel]
  X.sel=X[idx.sel,]
  Pred.X.sel=pred.X[idx.sel,]
  Pred.Y.sel=pred.Y[idx.sel] 
  
  score.sel=function(beta){
    h=apply(c((Y.sel-(X.sel-Pred.X.sel)%*%beta-Pred.Y.sel)/PI.sel)*(X.sel-Pred.X.sel),2,mean)
    return(h)
  }
  
  beta.sel=nleqslv(rep(0,p),score.sel,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  
  
  grad.sel=t(X.sel-Pred.X.sel)%*%diag(c(1/PI.sel))%*%(X.sel-Pred.X.sel)
  score.sel=c(Y.sel-(X.sel-Pred.X.sel)%*%beta.sel-Pred.Y.sel)*(X.sel-Pred.X.sel)
  score.sel2=t(score.sel)%*%diag(c(1/(PI.sel)^2+r.sel/(PI.sel)))%*%score.sel
  
  
  tmp.se=sqrt(diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel)))
  tmp.var=diag(ginv(grad.sel)%*%score.sel2%*%ginv(grad.sel))
  tmp.cp=abs(beta.sel-total0)<tmp.se*qnorm(1-0.05/2)
  return(list(est=beta.sel,se=tmp.se,var=tmp.var,cp=tmp.cp))
}


sampling_fun_indirect=function(r.sel,pi.sel,pred.X,pred.XM,pred.Y){
  idx.sel=PI.sel=c()
  idx.sel=sample(1:N,r.sel,T,prob=pi.sel);PI.sel=pi.sel[idx.sel]
  Y.sel=Y[idx.sel]
  X.sel=X[idx.sel,]
  XM.sel=XM[idx.sel,]
  
  Pred.X.sel=pred.X[idx.sel,]
  Pred.XM.sel=pred.XM[idx.sel,]
  Pred.Y.sel=pred.Y[idx.sel] 
  
  score.sel.total=function(beta){
    h=apply(c((Y.sel-(X.sel-Pred.X.sel)%*%beta-Pred.Y.sel)/PI.sel)*(X.sel-Pred.X.sel),2,mean)
    return(h)
  }
  beta.sel.total=nleqslv(rep(0,p),score.sel.total,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  
  score.sel.direct=function(beta){
    h=apply(c((Y.sel-(XM.sel-Pred.XM.sel)%*%beta-Pred.Y.sel)/PI.sel)*(XM.sel-Pred.XM.sel),2,mean)
    return(h)
  }
  
  beta.sel.xm=nleqslv(rep(0,p+q),score.sel.direct,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  
  
  grad.sel.total=t(X.sel-Pred.X.sel)%*%diag(c(1/PI.sel))%*%(X.sel-Pred.X.sel)/(N*r.sel)
  grad.sel.direct=t(XM.sel-Pred.XM.sel)%*%diag(c(1/PI.sel))%*%(XM.sel-Pred.XM.sel)/(N*r.sel)
  
  
  score.sel.sigma1=score.sel.sigma2=matrix(NA,nrow =r.sel ,ncol = p)
  for (int in 1:r.sel) {
    score.sel.sigma1[int,]=c(Y.sel-(XM.sel-Pred.XM.sel)%*%beta.sel.xm
                  -Pred.Y.sel)[int]*(ginv(grad.sel.total)%*%(X.sel-Pred.X.sel)[int,]-D%*%ginv(grad.sel.direct)%*%(XM.sel-Pred.XM.sel)[int,])
    score.sel.sigma2[int,]=c((X.sel-Pred.X.sel)%*%beta.sel.total-(XM.sel-Pred.XM.sel)%*%beta.sel.xm
                           )[int]*(ginv(grad.sel.total)%*%(X.sel-Pred.X.sel)[int,])
  }
  
  score.sel2.sigma1=score.sel2.sigma2=matrix(0,nrow =p ,ncol = p)
  for (int in 1:r.sel) {
    score.sel2.sigma1=score.sel2.sigma1+score.sel.sigma1[int,]%*%t(score.sel.sigma1[int,])*(1/(PI.sel[int])^2+r.sel/PI.sel[int])/(N^2*r.sel^2)
    score.sel2.sigma2=score.sel2.sigma2+score.sel.sigma2[int,]%*%t(score.sel.sigma2[int,])*(1/(PI.sel[int])^2+r.sel/PI.sel[int])/(N^2*r.sel^2)
  }
  Sigma.sel<-score.sel2.sigma1+score.sel2.sigma2
  
  
  tmp.se=sqrt(diag(Sigma.sel))
  tmp.var=diag(Sigma.sel)
  
  beta.sel=beta.sel.total-c(beta.sel.xm[1:p])
  tmp.cp=abs(beta.sel-indirect0)<tmp.se*qnorm(1-0.05/2)
  
  return(list(est=beta.sel,se=tmp.se,var=tmp.var,cp=tmp.cp))
}








