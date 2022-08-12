context("testing cpop")

##dyn.load("coeff.updateR.so")

########
# CPOP algorithm for finding the best segmentation of data for a change-in-slope model
#
# best is defined in terms of minimising
#
#  sum_{i=1}^n1/(sigma_i^2)  (y_i-f_i)^2+m*beta  -(*)
#
# where y_i are data, f_i is the fitted mean at point i; m is the number of changepoints
# and we consider all continuous piecewise-linear functions for f. Changepoints thus
# correspond to changes in slope.
#
# Input is y -- vector of data; and beta a positive constant that penalises additional
# changes; sigsquared -- an estimate of residual variance
#
# useCprune determine whether (faster) C code is used within the algorithm;
# printinteration determines whether updates of progress are printed.
#
# Output: list of minimum value of (*) together with the inferred changepoints. This list will include 0 and n as first/last entries.
#
########

####CHANGE-- INPUT x AS LIST OF LOCATIONS
####sigsquared an vary with x

CPOP.uneven.var<-function(y,x,beta,sigsquared=1,useCprune=FALSE,printiteration=FALSE){
  
  if(useCprune) dyn.load("prune2R.so")
  
  n<-length(y)
  if(length(sigsquared)!=n) sigsquared=rep(sigsquared[1],n)
  
  S<-0
  for(i in 1:n){
    S[i+1]<-S[i]+y[i]/sigsquared[i]
  }
 
  SS<-0
  for(i in 1:n){
    SS[i+1]<-SS[i]+y[i]^2/sigsquared[i]
  }  ##preprocessing
  
  ####ADDITIONAL PRE-PREOCESSING FOR UNEVEN COMPONENTS
  SX<-0
  for(i in 1:n){
    SX[i+1]<-SX[i]+x[i]/sigsquared[i]
  }
  SX2<-0
  for(i in 1:n){
    SX2[i+1]<-SX2[i]+x[i]^2/sigsquared[i]
  }
  SXY<-0
  for(i in 1:n){
    SXY[i+1]<-SXY[i]+x[i]*y[i]/sigsquared[i]
  }
  SP<-0
  for(i in 1:n){
    SP[i+1]<-SP[i]+1/sigsquared[i]
  }
  
  x0=c(2*x[1]-x[2],x)
  
  coeffs<-matrix(0,ncol=5,nrow=1) #first two columns are current time point and most recent changepoint, final three are coefficients for cost
  coeffs[1,5]<--beta
  coeffs[1,1:2]<-c(0,0)
  CPvec<-c("0") #vector storing changepoint values, not used in code but required as an output 
  lencoo<-c()
  lencur<-c()
 ######
  ## code minimise 1/(2sigma^2)*RSS rather than RSS/sigma^2
  ## hence need to have sigma^2 
 #### 
 sigsquared=sigsquared/2

  for(taustar in 1:n){
    new.CPvec<-paste(CPvec,taustar,sep=",")
    ##update coefficients --THIS HAS BEEN CHANGED FROM CPOP CODE
    ##CURRENTLY CHANGE ONLY IN R CODE VERSION
   #if(useC==FALSE){
     new.coeffs=coeff.update.uneven.var(coeffs,S,SXY,SS,SX,SX2,SP,x0,taustar,beta)
    #}    else if(useC==TRUE){
    #  new.coeffs=coeff.update.c(coeffs,S,SJ,SS,taustar,sigsquared,beta)
    #} else{stop("useC must be a TRUE or FALSE value")
    #}
   # if(taustar==2){browser()}
    if(taustar!=n){ #skip pruning on last step
    ###################################################pruning bit##########  
    if(length(new.coeffs[,1])>1){
      ##added###
      #keep1=prune1(new.coeffs,taustar) ##first pruning
      keep1=1:length(new.coeffs[,1])
      new.coeffs.p=new.coeffs[keep1,]
      new.CPvec=new.CPvec[keep1]
      ###########
      if(sum(keep1)>1){
       if(useCprune==F){
         keep2=prune2b(new.coeffs.p) }##find set of functions to keep
       else if(useCprune==T){
         keep2=prune2.c(new.coeffs.p) } 
       else{stop("useCprune must be a TRUE or FALSE value")
       }
       new.coeffs.p=new.coeffs.p[keep2,]
        new.CPvec=new.CPvec[keep2]
      }
    }else{
      new.coeffs.p=new.coeffs
    }
    ####PELT PRUNE############################
    if(taustar>2){
      keeppelt=peltprune(new.coeffs,beta)
      coeffs<-coeffs[keeppelt,]
      CPvec<-CPvec[keeppelt]
    }
    ##########################################
    }
    else{new.coeffs.p<-new.coeffs}
    CPvec<-c(CPvec,new.CPvec) #prunes both CPvec vector and coeffs matrix
    coeffs<-rbind(coeffs,new.coeffs.p)
    lencoo[taustar]<-length(coeffs[,1])
    lencur[taustar]<-length(new.coeffs.p)/5
    #####################################################
    if(printiteration==TRUE){
    if(taustar%%100==0) cat("Iteration ",taustar,"Functions-stored",lencoo[taustar],lencur[taustar],"\n")}
    else if(printiteration!=FALSE){stop("printiteration must be a TRUE or FALSE value")}
  }
  
  coeffscurr<-coeffs[coeffs[,1]==n,] #matrix of coeffs for end time t=n
  if(!is.matrix(coeffscurr)){coeffscurr<-t(as.matrix(coeffscurr))} #makes sure coeffscurr is in the right format
  ttemp<-coeffscurr[,5]-(coeffscurr[,4]^2)/(4*coeffscurr[,3])
  mttemp<-min(ttemp)
  num<-which(ttemp==mttemp)
  
  CPveccurr<-CPvec[coeffs[,1]==n]
  CPS<-eval(parse(text=paste("c(",CPveccurr[num],")")))
#####
  ##code has an additional additive factor of (n/2) * log(2*pi*sigma^2) in cost
  ## hence we remove this so mttemp is the minimum of the correct cost
#####
  #mttemp=mttemp - (n/2)*log(2*pi*sigsquared)
  return(list(min.cost=mttemp,changepoints=CPS)) #return min cost and changepoints
}


########################################################################################
###################### coeff update#####################################################
### NEW VERSION -- UNEVEN LOCATIONS AND VARYING MEAN
###avoids loop
########################################################################################

coeff.update.uneven.var=function(coeffs,S,SXY,SS,SX,SX2,SP,x0,taustar,beta){
  
  coeff.new<-coeffs
  coeff.new[,2]=coeffs[,1] 
  coeff.new[,1]<-taustar
  
  sstar<-coeff.new[,2]
  Xs<-x0[sstar+1]
  Xt<-x0[taustar+1]
  seglen=Xt-Xs
  
  n.obs=taustar-sstar
  A<-(SX2[taustar+1]-SX2[sstar+1]-2*Xs*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xs^2)/(seglen^2)
  B<- 2*( (Xt+Xs)*(SX[taustar+1]-SX[sstar+1])-(SP[taustar+1]-SP[sstar+1])*Xt*Xs-(SX2[taustar+1]-SX2[sstar+1]))/(seglen^2)
  C<-(-2)/(seglen)*(SXY[taustar+1]-SXY[sstar+1]-Xs*(S[taustar+1]-S[sstar+1]))
  D<- (SS[taustar+1]-SS[sstar+1])
  E<-(-2)/(seglen)*(Xt*(S[taustar+1]-S[sstar+1])-(SXY[taustar+1]-SXY[sstar+1]))
  FF<-(SX2[taustar+1]-SX2[sstar+1]-2*Xt*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xt^2)/(seglen^2)
  
  m=length(sstar)
  ind1=(1:m)[FF==0 & coeffs[,3]==0 & B==0]
  ind2=(1:m)[FF==0 & coeffs[,3]==0 & B!=0]
  ind3=(1:m)[!(FF==0 & coeffs[,3]==0)]
  if(length(ind1)>0){
    coeff.new[ind1,5]<-coeffs[ind1,5]+D[ind1]+beta
    coeff.new[ind1,4]<-C[ind1]
    coeff.new[ind1,3]<-A[ind1] 
  }
  
  if(length(ind2)>0){
    coeff.new[ind2,5]<-coeffs[ind2,5]+(-E[ind2]-coeffs[ind2,4])/B[ind2]+beta
    coeff.new[ind2,4]<-0
    coeff.new[ind2,3]<-0
  }
  
  if(length(ind3)>0){
    coeff.new[ind3,5]<-coeffs[ind3,5]+D[ind3]-(coeffs[ind3,4]+E[ind3])^2/(4*(coeffs[ind3,3]+FF[ind3]))+beta
    coeff.new[ind3,4]<-C[ind3]-(coeffs[ind3,4]+E[ind3])*B[ind3]/(2*(coeffs[ind3,3]+FF[ind3]))
    coeff.new[ind3,3]<-A[ind3]-(B[ind3]^2)/(4*(coeffs[ind3,3]+FF[ind3]))
  }
  
  return(coeff.new)
}





##########################################################################################################
##second pruning
## again x is matrix of quadratics
##version to avoid nested loops
##########################################################################################################

prune2b=function(x){
  Sets<-list()
  n=length(x[,1])
  vec=(1:n)
  
  tcurr= -Inf
  
  whichfun<-which(x[,3]==min(x[,3])) #which element of vec gives min value at -Infinity--smallest theta^2 coeff; then largest theta coeff; then smallest constant
  whichfun<-whichfun[which(x[whichfun,4]==max(x[whichfun,4]))]
  whichfun<-whichfun[which(x[whichfun,5]==min(x[whichfun,5]))]
  
  Sets[[whichfun]]<-c(tcurr)
  diffcoeffs=matrix(NA,nrow=n,ncol=3)
  intercepts=rep(NA,n)
  disc=rep(NA,n)
  while(length(vec)>1){ #while functions being considered is bigger than 1
    intercepts[1:n]<-NA
    diffcoeffs[1:(length(vec)),]<-t(t(x[vec,3:5])-x[whichfun,3:5]) #difference between coeffs at i and current function
    disc[1:(length(vec))]<-diffcoeffs[1:(length(vec)),2]^2-4*diffcoeffs[1:(length(vec)),1]*diffcoeffs[1:(length(vec)),3] #discriminent of difference quad
    
    ind1=(1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]==0] ##disc>0 for quadratic to cross.
    ind2=(1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]!=0] ##disc>0 for quadratic to cross.
    
    if(length(ind1)>0){
      r1= - diffcoeffs[ind1,3]/diffcoeffs[ind1,2]
      if(sum(r1>tcurr)>0){
        intercepts[ind1[r1>tcurr]]= r1[r1>tcurr]
      }
    }
    if(length(ind2)>0){
      r1=(-diffcoeffs[ind2,2]-sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
      r2=(-diffcoeffs[ind2,2]+sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
      ##only want roots if > tcurr
      if(sum(r1>tcurr)>0){
        intercepts[ind2[r1>tcurr]]=r1[r1>tcurr]
      }
      if(sum(r1<=tcurr & r2>tcurr)>0){
        intercepts[ind2[r1<=tcurr & r2>tcurr]]=r2[r1<=tcurr & r2>tcurr]  
      }
    }
    
    loggy<-!is.na(intercepts)
    loggy[vec==whichfun]<-T
    if(!sum(!is.na(intercepts))==0){ #if at least one intercept value is not na
      tcurr<-min(intercepts,na.rm=T) #change tcurr to first intercept
      whichfunnew<-vec[which(intercepts==tcurr)[1]] #whichfunnew is set as value which first intercept occurs     
      Sets[[whichfun]]<-c(Sets[[whichfun]],tcurr) #add intercept to current function opt interval (to close it)
      if(whichfunnew>length(Sets)){Sets[[whichfunnew]]<-c(tcurr)}else{
        Sets[[whichfunnew]]<-c(Sets[[whichfunnew]],tcurr)} #add intercept to new fucntion interval (opening it)
      whichfun<-whichfunnew #change current function to new function
      
    }
    vec<-vec[loggy[(1:length(vec))]]
    
  }
  Sets[[whichfun]]<-c(Sets[[whichfun]],Inf)
  
  
  output1 <- do.call(rbind,lapply(Sets,length))
  
  return(which(output1[,1]!=0))
}



###########################################################################################################
###################################PELT pruning function###################################################
###########################################################################################################

peltprune=function(x,beta){
  minx<-x[,5]-x[,4]^2/(4*x[,3])
  return(which(minx<=(min(minx)+2*beta)))  
}

###########################################################################################################
##################coverts null to na#######################################################################
###########################################################################################################

null2na<-function(vec){
  if(is.null(vec)){vec<-NA}
  return(vec)
}


################end#######################


#######FUNCTION FOR FITTING
### Input y-data; x- locations; out -- output from CPOP; and sigsquared -- vector of variances
###
CPOP.fit=function(y,x,out.changepoints,sigsquared){
  n=length(y)
  if(length(sigsquared)!=n) sigsquared=rep(sigsquared[1],n)
  p=length(out.changepoints)
  W=diag(sigsquared^-1)
  X=matrix(NA,nrow=n,ncol=p)
  X[,1]=1
  X[,2]=x-x[1]
  if(p>2){
    for(i in 2:(p-1)){
      X[,i+1]=c(rep(0,out.changepoints[i]-1),x[(out.changepoints[i]):n]-x[out.changepoints[i]])
    }
  }
  XTX=t(X)%*%W%*%X
  beta=as.vector(solve(XTX)%*%t(X)%*%W%*%y)
  fit=X%*%beta
  residuals=y-fit
  return(list(fit=fit,residuals=residuals,X=X,pars=beta))
}

##### function to simulate mean
## x data locations
## mu(x)=sum_{i=1}^K a_k max(x-tau_k,0)
##where changepoints is the vector tau_{1:K}
## and slope is the vetor of a_{1:K} -- the change in slope
change.in.slope.mean=function(x,changepoints,change.slope){
  K=length(changepoints)
  mu=rep(0,length(x))
  for(k in 1:K) mu=mu+change.slope[k]*pmax(x-changepoints[k],0)
  return(mu)
}


test_that("test that cpop predicts the correct RSS and changepoints using default parameter values",
{
   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y<-mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=1)
   fit=CPOP.fit(y,x,out$changepoints,1)
   RSS<-sum(fit$residuals^2)
   cpop.res<-cpop(y,x)
   cpop.RSS<-sum(fitted(cpop.res)$RSS)
   expect_equal(cpop.RSS,RSS)
   expect_equal(changepoints(cpop.res)$location,x[out$changepoints[2:4]])
})






test_that("test 1 - test that cpop predicts the correct RSS and changepoints using default parameter values",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y<-mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=1)
   fit=CPOP.fit(y,x,out$changepoints,1)
   RSS<-sum(fit$residuals^2)    
   cpop.res<-cpop(y,x)
   cpop.RSS<-sum(fitted(cpop.res)$RSS)
   expect_equal(cpop.RSS,RSS)
   expect_equal(changepoints(cpop.res)$location,x[out$changepoints[2:4]])
})

test_that("test 2 - test that cpop predicts the correct RSS and changepoints using non default beta value",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2,sigsquared=1)
   fit=CPOP.fit(y,x,out$changepoints,1)
   RSS<-sum(fit$residuals^2)
   cpop.res<-cpop(y,x,beta=2)
   cpop.RSS<-sum(fitted(cpop.res)$RSS)
   expect_equal(cpop.RSS,RSS)
   expect_equal(changepoints(cpop.res)$location,x[out$changepoints[2:24]])
})





test_that("test 3 - test that fitted predicts the correct RSS",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=1)
   out$changepoints
   #[1]   0  22  52 95 200
   out$min.cost
   #[1] 199.0554
   fit=CPOP.fit(y,x,out$changepoints,1)
   #CHECK
   # RSS<-sum(fit$residuals^2)+2*log(length(x))*(length(out$changepoints)-2)
   RSS<-sum(fit$residuals^2)
   cpop.res<-cpop(y,x,sd=1)
   expect_equal(sum(fitted(cpop.res)$RSS),RSS)
})


test_that("test 4 - test that fitted predicts the correct RSS",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=1)
   out$changepoints
   #[1]   0  22  52 95 200
   out$min.cost
   #[1] 199.0554
   fit=CPOP.fit(y,x,out$changepoints,1)
   #CHECK
   # RSS<-sum(fit$residuals^2)+2*log(length(x))*(length(out$changepoints)-2)
   RSS<-sum(fit$residuals^2)
   cpop.res<-cpop(y,x,sd=1)
   expect_equal(sum(fitted(cpop.res)$RSS),RSS)
})


test_that("test 5 - test that results from fitted can be used to calculate the cost correctly",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=1)
   out$changepoints
   #[1]   0  22  52 95 200
   out$min.cost
   #[1] 199.0554
   fit=CPOP.fit(y,x,out$changepoints,1)
   #CHECK
   cost<-sum(fit$residuals^2)+2*log(length(x))*(length(out$changepoints)-2)
   cpop.res<-cpop(y,x,sd=1)
   expect_equal(cost(cpop.res),cost)
})



test_that("test 6 - test default value of sd",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=1)
   out$changepoints
   #[1]   0  22  52 95 200
   out$min.cost
   #[1] 199.0554
   fit=CPOP.fit(y,x,out$changepoints,1)
   #CHECK
   cost<-sum(fit$residuals^2)+2*log(length(x))*(length(out$changepoints)-2)
   cpop.res<-cpop(y,x)
   expect_equal(cost(cpop.res),cost)
})


test_that("test 7 - test non default values of sd",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=2)
   out$changepoints
   #[1]   0  22  52 95 200
   out$min.cost
   #[1] 199.0554
   fit=CPOP.fit(y,x,out$changepoints,2)
   #CHECK
   cost<-sum(fit$residuals^2/2)+2*log(length(x))*(length(out$changepoints)-2)
   cpop.res<-cpop(y,x,sd=sqrt(2))
   expect_equal(cost(cpop.res),cost)
   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   mu=change.in.slope.mean(x,changepoints,change.slope)
   y=mu+rnorm(200)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=4)
   out$changepoints
   #[1]   0  22  52 95 200
   out$min.cost
   #[1] 199.0554
   fit=CPOP.fit(y,x,out$changepoints,4)
   #CHECK
   cost<-sum(fit$residuals^2/4)+2*log(length(x))*(length(out$changepoints)-2)
   cpop.res<-cpop(y,x,sd=2)
   expect_equal(cost(cpop.res),cost)
})


test_that("test 8 - test for the effects of setting minseglen greater than shortest distance between changepoints",
{

   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   y<-simchangeslope(x,changepoints,change.slope,1)
   out=CPOP.uneven.var(y,x,beta=2*log(length(x)),sigsquared=4)
   out$changepoints
   out$min.cost
   fit=CPOP.fit(y,x,out$changepoints,4)
   cost<-sum(fit$residuals^2/4)+2*log(length(x))*(length(out$changepoints)-2)
   cpop.res<-cpop(y,x,sd=2)
   expect_equal(cost(cpop.res),cost)
   cpop.minseglen.res<-cpop(y,x,sd=2,minseglen=22)
   expect_equal(cost(cpop.res),cost(cpop.minseglen.res))
   cpop.minseglen.res<-cpop(y,x,sd=2,minseglen=23)
   expect_false(isTRUE(all.equal(cost(cpop.res),cost(cpop.minseglen.res))))
   cpop.minseglen.res<-cpop(y,x,sd=2,minseglen=26)
   expect_equal(changepoints(cpop.minseglen.res)$location[1],110)
})


test_that("test 9 - test use of non default value for grid",
{

   set.seed(1)
   x<-1:200
   changepoints<-c(0,25,50,100)
   change.slope<-c(0.2,-0.3,0.2,-0.1)
   y<-simchangeslope(x,changepoints,change.slope,1)
   cpop.res<-cpop(y,x,sd=2)
   cpop.grid.res<-cpop(y,x,grid=c(0.5,1.5,99.0),sd=2)
   expect_false(isTRUE(all.equal(cost(cpop.res),cost(cpop.grid.res))))
   expect_equal(changepoints(cpop.grid.res)$location[1],99)
})



test_that("test 10 - test use of non unit locations data",
{

   set.seed(0)
   x <- seq(0,1,0.01)
   n <- length(x)
   sigma <- rep(0.1,n)
   mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
   y <- rnorm(n,mu,sigma)
   # use the locations in x
   out=CPOP.uneven.var(y,x,beta=2*log(length(y)),sigsquared=0.01)
   cpop.res <- cpop(y,x,beta=2*log(length(y)),sd=0.1)
   fit=CPOP.fit(y,x,out$changepoints,0.1)
   RSS<-sum(fit$residuals^2)
   expect_equal(sum(fitted(cpop.res)$RSS),RSS)
})





test_that("test 11 - default x values ",
{
 # generate some test data
 set.seed(0)
 x <- seq(0,1,0.01)
 n <- length(x)
 sigma <- rep(0.1,n)
 mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
 y <- rnorm(n,mu,sigma)
 
 # use the locations in x
 res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
 fitted.base <- fitted(res)
 cpts.base <- changepoints(res)
 estimates.base <- estimate(res,x)

 # without locations (note explicit paramater names)
 res <- cpop(y,beta=2*log(length(y)),sd=sigma)
 expect_equal(estimate(res,1:length(y)-1)$y_hat,estimates.base$y_hat)
 expect_equal(fitted(res)[,7],fitted.base[,7])
 
})


test_that("test 12 - shifted x values ",
{
 # generate some test data
 set.seed(0)
 x <- seq(0,1,0.01)
 n <- length(x)
 sigma <- rep(0.1,n)
 mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
 y <- rnorm(n,mu,sigma)
 
 # use the locations in x
 res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
 fitted.base <- fitted(res)
 cpts.base <- changepoints(res)
 estimates.base <- estimate(res,x)

 x <- x + 1
 res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
 expect_equal(estimate(res,x)$y_hat,estimates.base$y_hat)
 expect_equal(fitted(res)[,7],fitted.base[,7])

 x <- x - 1
 for(i in 1:10)
 {
   x <- x + rnorm(1,0,10)
   res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
   expect_equal(estimate(res,x)$y_hat,estimates.base$y_hat)
   expect_equal(fitted(res)[,7],fitted.base[,7])
 }

})
















