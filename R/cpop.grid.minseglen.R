

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

cpop.grid.minseglen<-function(y,x=NULL,grid=NULL,minseg=0,beta,sigsquared=1,printiteration=FALSE,PRUNE.APPROX=FALSE)
{


CPOP.grid.minseg<-function(y,x=NULL,grid=NULL,minseg=0,beta,sigsquared=1,printiteration=FALSE,PRUNE.APPROX=FALSE){

  n<-length(y)
  if(is.null(x)) x <-  1:n
  if(is.null(grid)) grid <- x  ###SET grid to x if not defined
  if(max(grid)<max(x)) grid <- c(grid,max(x))
  if(length(sigsquared)!=n) sigsquared <- rep(sigsquared[1],n)
  
  if(minseg<0) minseg <- 0
  
  ngrid <- length(grid) ##ngrid is the number of grid points -- below n becomes ngrid throughout
  
  ##CHANGE HERE -- x0 is grid with an addition of a 0th grid point. Defined in terms of grid not x
  x0 <- c(2*x[1]-x[2],grid)
  
  if(minseg*2> max(x0)-min(x0)) stop("Minseg is greater than half the grid width => no changepoints possible")
  
  ##THESE UPDATES ARE CHANGED.
  S<-0
  SS<-0
  SX<-0
  SX2<-0
  SXY<-0
  SP<-0
  for(i in 1:ngrid){
    ##set index to be which observations are in the interval
    index <- (1:n)[x>x0[i]&x<=x0[i+1]]
    if(length(index)>0){##change in summaries -- sum over observations
      S[i+1]<-S[i]+sum(y[index]/sigsquared[index])
      SS[i+1]<-SS[i]+sum(y[index]^2/sigsquared[index])
      SX[i+1]<-SX[i]+sum(x[index]/sigsquared[index])
      SX2[i+1]<-SX2[i]+sum(x[index]^2/sigsquared[index])
      SXY[i+1]<-SXY[i]+sum(x[index]*y[index]/sigsquared[index])
      SP[i+1]<-SP[i]+sum(1/sigsquared[index])
    }else{ ##no observations so no change in summaries
      S[i+1]<-S[i]
      SS[i+1]<-SS[i]
      SX[i+1]<-SX[i]
      SX2[i+1]<-SX2[i]
      SXY[i+1]<-SXY[i]
      SP[i+1]<-SP[i]
    }
  }
  #######
  
  coeffs<-matrix(0,ncol=5,nrow=1) #first two columns are current time point and most recent changepoint, final three are coefficients for cost
  coeffs[1,5]<--beta
  coeffs[1,1:2]<-c(0,0)
  CPvec<-c("0") #vector storing changepoint values, not used in code but required as an output 
  
  for(taustar in 1:ngrid){ ##n->ngrid is the number of iterations
    if(grid[taustar]-x0[1]>=minseg){ ##ONLY START ADDING CHANGES IF FIRST SEG > minseg
      ##check on which values of most recent change are possible
      index <- (1:length(CPvec))[grid[taustar]-x0[coeffs[,1]+1]>=minseg]     
      new.CPvec<-paste(CPvec[index],taustar,sep=",")
      ##update coefficients 
      new.coeffs <- coeff.update.uneven.var(coeffs[index,,drop=FALSE],S,SXY,SS,SX,SX2,SP,x0,taustar,beta)
      new.coeffs.p <- new.coeffs  
      if(taustar!=ngrid){##n->ngrid #skip pruning on last step
    ###################################################pruning bit##########  
        if(length(new.coeffs[,1])>1){
      ##added###
      ###########
        # keep2 <- prune2b(new.coeffs.p) ##find set of functions to keep
	keep2 <- prune2.c(new.coeffs.p)
        new.coeffs.p <- new.coeffs.p[keep2,]
        new.CPvec <- new.CPvec[keep2]
        }
    ####PELT PRUNE############################
        if(PRUNE.APPROX){
          j <- which.max(x0[-1]*(x0[taustar+2]-x0[-1]>=minseg))
          ###IDEA IS WE CAN PELT PRUNE BASED ON COEFFS FOR CHANGE AT j =>
          index <- (1:length(CPvec))[grid[j]-x0[coeffs[,1]+1]>=minseg]
          if(length(index)>1){
            coeffs.pelt <- coeff.update.uneven.var(coeffs[index,,drop=FALSE],S,SXY,SS,SX,SX2,SP,x0,j,beta)
            keeppelt<-peltprune(coeffs.pelt,1.5*beta) ##this holds for new.coeffs by choice above.
            if(length(keeppelt)<length(index)){
              coeffs<-coeffs[-index[-keeppelt],]
              CPvec<-CPvec[-index[-keeppelt]]
            }
          }
         }
        }
    ##########################################
    CPvec<-c(CPvec,new.CPvec) 
    coeffs<-rbind(coeffs,new.coeffs.p)
    }
    #####################################################
    if(printiteration==TRUE){
      if(taustar%%100==0) cat("Iteration ",taustar,"\n")}
      else if(printiteration!=FALSE){stop("printiteration must be a TRUE or FALSE value")}
  }
  coeffscurr<-coeffs[coeffs[,1]==ngrid,] #matrix of coeffs for end time t=n ##n->ngrid
  if(!is.matrix(coeffscurr)){coeffscurr<-t(as.matrix(coeffscurr))} #makes sure coeffscurr is in the right format
##HACK here to make ttemp correct if coeffscurr[,3]==0
  ttemp<-coeffscurr[,5]
  index <- (1: (length(coeffscurr)/5))[coeffscurr[,3]>0] 
  if(length(index)>0){ 
    ttemp[index]<-coeffscurr[index,5]-(coeffscurr[index,4]^2)/(4*coeffscurr[index,3])}
  mttemp<-min(ttemp)
  num<-which(ttemp==mttemp)
  
  CPveccurr<-CPvec[coeffs[,1]==ngrid] ##n->ngrid
  CPS<-eval(parse(text=paste("c(",CPveccurr[num],")")))

  return(list(min.cost=mttemp,changepoints.index=CPS,changepoints=grid[CPS])) #return min cost and changepoints
}


########################################################################################
###################### coeff update#####################################################
### NEW VERSION -- UNEVEN LOCATIONS AND VARYING MEAN
###avoids loop
########################################################################################

coeff.update.uneven.var <- function(coeffs,S,SXY,SS,SX,SX2,SP,x0,taustar,beta){
  
  coeff.new<-coeffs
  coeff.new[,2] <- coeffs[,1] 
  coeff.new[,1]<-taustar
  
  sstar<-coeff.new[,2]
  Xs<-x0[sstar+1]
  Xt<-x0[taustar+1]
  seglen <- Xt-Xs
  
  A<-(SX2[taustar+1]-SX2[sstar+1]-2*Xs*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xs^2)/(seglen^2)
  B<- 2*( (Xt+Xs)*(SX[taustar+1]-SX[sstar+1])-(SP[taustar+1]-SP[sstar+1])*Xt*Xs-(SX2[taustar+1]-SX2[sstar+1]))/(seglen^2)
  C<-(-2)/(seglen)*(SXY[taustar+1]-SXY[sstar+1]-Xs*(S[taustar+1]-S[sstar+1]))
  D<- (SS[taustar+1]-SS[sstar+1])
  E<-(-2)/(seglen)*(Xt*(S[taustar+1]-S[sstar+1])-(SXY[taustar+1]-SXY[sstar+1]))
  FF<-(SX2[taustar+1]-SX2[sstar+1]-2*Xt*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xt^2)/(seglen^2)
  
  m <- length(sstar)
  ind1 <- (1:m)[FF==0 & coeffs[,3]==0 & B==0]
  ind2 <- (1:m)[FF==0 & coeffs[,3]==0 & B!=0]
  ind3 <- (1:m)[!(FF==0 & coeffs[,3]==0)]
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



############################################################################################################
##first prune of functions
## xx is matrix of functions
############################################################################################################
prune1 <- function(xx,taustar){
  min.vals<-xx[,5]-(xx[,4]^2)/(2*xx[,3])
  m <- min(min.vals[xx[,2]==taustar-1])
  keep<-c(rep(T,length(min.vals[xx[,2]!=taustar-1])),(min.vals[xx[,2]==taustar-1])<=m)
  return(keep) ##keep only those with a smaller minimum.
}


##########################################################################################################
##second pruning
## again x is matrix of quadratics
##version to avoid nested loops
##########################################################################################################

prune2b <- function(x){
  Sets<-list()
  n <- length(x[,1])
  vec <- (1:n)
  
  tcurr <-  -Inf
  
  whichfun<-which(x[,3]==min(x[,3])) #which element of vec gives min value at -Infinity--smallest theta^2 coeff; then largest theta coeff; then smallest constant
  whichfun<-whichfun[which(x[whichfun,4]==max(x[whichfun,4]))]
  whichfun<-whichfun[which(x[whichfun,5]==min(x[whichfun,5]))]
  whichfun <- whichfun[1] #####Introduced this to avoid multiple minimima
  
  Sets[[whichfun]]<-c(tcurr)
  diffcoeffs <- matrix(NA,nrow=n,ncol=3)
  intercepts <- rep(NA,n)
  disc <- rep(NA,n)
  while(length(vec)>1){ #while functions being considered is bigger than 1
    intercepts[1:n]<-NA
    diffcoeffs[1:(length(vec)),]<-t(t(x[vec,3:5])-x[whichfun,3:5]) #difference between coeffs at i and current function
    disc[1:(length(vec))]<-diffcoeffs[1:(length(vec)),2]^2-4*diffcoeffs[1:(length(vec)),1]*diffcoeffs[1:(length(vec)),3] #discriminent of difference quad
    
    ind1 <- (1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]==0] ##disc>0 for quadratic to cross.
    ind2 <- (1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]!=0] ##disc>0 for quadratic to cross.
    
    if(length(ind1)>0){
      r1 <-  - diffcoeffs[ind1,3]/diffcoeffs[ind1,2]
      if(sum(r1>tcurr)>0){
        intercepts[ind1[r1>tcurr]] <- r1[r1>tcurr]
      }
    }
    if(length(ind2)>0){
      r1 <- (-diffcoeffs[ind2,2]-sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
      r2 <- (-diffcoeffs[ind2,2]+sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
      ##only want roots if > tcurr
      if(sum(r1>tcurr)>0){
        intercepts[ind2[r1>tcurr]] <- r1[r1>tcurr]
      }
      if(sum(r1<=tcurr & r2>tcurr)>0){
        intercepts[ind2[r1<=tcurr & r2>tcurr]] <- r2[r1<=tcurr & r2>tcurr]  
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

#############################################################################
#########second pruning in C##################################################
####################################################################################


prune2.c<-function(x){
  nrows<-dim(x)[[1]]
  # coutput<- .C("prune2R",as.double(as.vector((t(x)))),as.integer(nrows),Sets=integer(nrows))
  coutput<-prune(as.vector((t(x))),nrows)
  #result<-which(coutput$Sets==1)
  result<-which(coutput==1)
  return(result)
}

0

###########################################################################################################
###################################PELT pruning function###################################################
###########################################################################################################

peltprune <- function(x,beta){
  if(length(x)==5) return(1)
  minx<-x[,5]
  index <- (1:dim(x)[1])[x[,3]>0]
  if(length(index)>0) minx[index]<-x[index,5]-x[index,4]^2/(4*x[index,3]) ####NEED THIS TO BE x[,5] if x[,4]=x[,3]=0
  return(which(minx<=(min(minx)+2*beta)))  
}

###########################################################################################################
##################coverts null to na#######################################################################
###########################################################################################################


out<-CPOP.grid.minseg(y,x,grid,minseg,beta,sigsquared,FALSE,FALSE)
out$changepoints<-c(x[1],out$changepoints)
return(out)


################end#######################

}