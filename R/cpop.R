
.cpop.class<-setClass("cpop.class",representation(min.cost="numeric",
                                                  changepoints="numeric",
						  y="numeric",
						  x="numeric",
						  y_hat="numeric",
						  residuals="numeric",
						  parameters="numeric",
						  beta="numeric",
						  sigsquared="numeric"))

cpop.class<-function(y,x,beta,sigsquared,min.cost,changepoints,y_hat,residuals,parameters)
{
    .cpop.class(y=y,
                x=x,
		min.cost=min.cost,
		changepoints=changepoints,
		y_hat=as.numeric(y_hat),
		residuals=as.numeric(residuals),
		parameters=parameters,
		beta=beta,
		sigsquared=sigsquared)
}


#' A function for generating simulated multivariate data
#'
#' @name simulate
#'
#' @description Generates simulated data for use with \code{\link{cpop}}.
#' 
#' @param x A numeric vector of containing the locations of the data.
#' @param changepoints A numeric vector of changepoint locations.
#' @param change.slope A numeric vector indicating the change in slope at each changepoint. The initial slope is assumed to be 0.
#' @param sigma The residual standard deviation. Can be a single numerical value or a vector of values for the case of varying residual standard deviation.Default value is 1.
#' @return A vector of simulated y values.
#'
#' @examples
#' library(cpop)
#' set.seed(1)
#' changepoints=c(0,25,50,100)
#' change.slope=c(0.2,-0.3,0.2,-0.1)
#' x=1:200
#' sig=1+x/200
#' y<-simulate(x,changepoints,change.slope,sig)
#' res<-cpop(y,x,2*log(length(x)),sig^2)
#' summary(res)
#' plot(res)
#'
#' @export
simulate<-function(x,changepoints,change.slope,sigma=1)
{
  K=length(changepoints)
  mu=rep(0,length(x))
  for(k in 1:K)
  {
     mu=mu+change.slope[k]*pmax(x-changepoints[k],0)
  }
  y=mu+rnorm(length(x),0.0,sigma)
  return(y)
}


#' Visualisation of Changepoint Locations and Data
#'
#' @name plot
#'
#' @description Plot methods for S4 objects returned by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param x An instance of an cpop S4 class produced by \code{\link{cpop}}.
#' 
#' @return A ggplot object.
#'
#' @rdname plot-methods
#'
#' @aliases plot,cpop.class-method
#' 
#' @export 
setMethod("plot",signature=list("cpop.class"),function(x)
{
  obj <- x
  df <- data.frame("x"=obj@x,"y"=obj@y)
  cpts<-obj@changepoints
  # appease ggplot2
  y <- NULL
  p <- ggplot(data=df, aes(x=x, y=y))
  p <- p + geom_point(alpha=0.3)
  if(length(cpts) > 2)
  {
   cpts<-cpts[2:(length(cpts)-1)]
   for(cpt in cpts)
   {
     p <- p + geom_vline(xintercept = obj@x[cpt],color="red")
   }
  }
  # appease ggplot2
  xs <- ys <- xends <- yends <- NULL  
  cpts<-obj@changepoints
  cpts[1]<-1
  df<-data.frame("xs"=obj@x[cpts][1:(length(obj@x[cpts])-1)],
	         "ys"=obj@y_hat[cpts][1:(length(obj@y_hat[cpts])-1)],
		 "xends"=obj@x[cpts][2:length(obj@x[cpts])],
		 "yends"=obj@y_hat[cpts][2:length(obj@y_hat[cpts])])
  p <- p + geom_segment(data=df,aes(x=xs,y=ys,xend=xends,yend=yends))
  p <- p + theme_bw()
  return(p)
})



#' Summary of cpop Analysis.
#'
#' @name summary
#'
#' @description Summary method for results produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @rdname summary-methods
#'
#' @aliases summary,cpop.class-method
#'
#' @export
setMethod("summary",signature=list("cpop.class"),function(object)
{
  cat('\n',"cpop analysis with n = ",length(object@x)," and penalty (beta)  = ",object@beta,'\n\n',sep="")
  if(length(object@changepoints) == 2)
  {
    cat("No changepoints detected",'\n\n',sep="")
  }
  else
  {
    msg<-paste(length(object@changepoints)-2," ",sep="")
    if(length(object@changepoints) == 3)
     {
       msg<-paste(msg," changepoint",sep="")
     }
     else
     {
       msg<-paste(msg," changepoints",sep="")
     }
     msg<-paste(msg," detected at :",'\n',sep="")
     cat(msg)
     msg<-"index : "
     for(cpt in object@changepoints[2:(length(object@changepoints)-1)])
     {
       msg<-paste(msg,cpt,sep=" ")
     }
     msg<-paste(msg,'\n',sep="")
     cat(msg)
     msg<-"location : "
     for(location in object@x[object@changepoints[2:(length(object@changepoints)-1)]])
     {
       msg<-paste(msg,location,sep=" ")
     }
     msg<-paste(msg,'\n',sep="")
     cat(msg)
  }
  df <- fitted(object)
  cat("fitted values : ",'\n',sep="")
  print(df)
  cat('\n',"overall RSS = ",sum(df$RSS),'\n',sep="")
  cat("cost = ",sum(object@residuals^2/object@sigsquared)+2*log(length(object@x))*(length(object@changepoints)-2),'\n',sep="")
  invisible()
})


#' Displays an S4 object produced by cpop.
#'
#' @name show
#'
#' @description Displays an S4 object produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @rdname show-methods
#'
#' @aliases show,cpop.class-method
#' 
#'
#' @export
if(!isGeneric("show")) {setGeneric("show",function(object) {standardGeneric("show")})}
setMethod("show",signature=list("cpop.class"),function(object)
{
    summary(object)
    invisible()
})


#' Extract Model Fitted Values
#'
#' @name fitted
#'
#' @description Extracts the fitted valus produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @return A data frame containing the endpoint coordinates for each line segment fitted between
#' the detected changepoints. The data frame also contains the gradient and intercept values
#' for each segment and the corresponding residual sum of squares (RSS)..
#'
#' @rdname fitted-methods
#'
#' @aliases fitted,cpop.class-method
#'
#' @export
if(!isGeneric("fitted")) {setGeneric("fitted",function(object) {standardGeneric("fitted")})}
setMethod("fitted",signature=list("cpop.class"),
          function(object)
          {
	  cpts <- object@changepoints
          cpts[1] <- 1
  	  df<-data.frame("x0"=object@x[cpts][1:(length(object@x[cpts])-1)],
			 "y0"=object@y_hat[cpts][1:(length(object@y_hat[cpts])-1)],
		         "x1"=object@x[cpts][2:length(object@x[cpts])],
		         "y1"=object@y_hat[cpts][2:length(object@y_hat[cpts])])
          df <- cbind(df,data.frame("gradient"=(df$y1 - df$y0)/(df$x1 - df$x0)))
	  df <- cbind(df,data.frame("intercept"=df$y1 - df$gradient * df$x1))
	  cpts <- object@changepoints
          rss<-as.numeric(Map(function(i,j) sum(object@residuals[i:j]*object@residuals[i:j]),cpts[1:(length(cpts)-1)]+1,cpts[2:length(cpts)]))
	  df <- cbind(df,data.frame("RSS"=rss))
	  return(df)
          })	      


#' Changepoint Locations
#'
#' @name changepoints
#'
#' @description Creates a data frame containing the locations of the changepoints in terms of the index of the data and the value of the location at that index.
#'
#' @docType methods
#'
#' @rdname changepoints-methods
#'
if(!isGeneric("changepoints")) {setGeneric("changepoints",function(object) {standardGeneric("changepoints")})}

#' @name changepoints
#' @param object  An instance of an cpop S4 class produced by \code{\link{cpop}}.
#' 
#' @return A data frame.
#' 
#' @rdname changepoints-methods
#'
#' @aliases changepoints,cpop.class-method
#'
#' @examples
#'
#' library(cpop)
#' # generate some test data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sigma <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sigma)
#'
#' # use the locations in x
#' res <- cpop(y,x,2*log(length(y)),0.1)
#' changepoints(res)
#'
#' @export
setMethod("changepoints",signature=list("cpop.class"),
          function(object)
          {
	      if(length(object@changepoints) > 2)
	      {
	      	      df <- data.frame("index"=object@changepoints[2:(length(object@changepoints)-1)],"location"=object@x[object@changepoints[2:(length(object@changepoints)-1)]])	
	      }
	      else
	      {
	      	      df <- data.frame("index"=integer(0),"locaion"=numeric(0))	
	      }
	      return(df)
          })	      


#' cpop
#'
#'  Algorithm for finding the best segmentation of data for a change-in-slope model.
#' 
#' @param y A vector of length n containing the data.
#' @param x A vector of length n containing the locations of y. Default value is NULL, in which case the locations \code{x = 1:length(y)} are assumed.
#' @param beta A positive real value for the penalty incurred for adding a changepoint (prevents over-fitting).
#' @param sigsquared Estimate of residual variance. Can be a single numerical value or a vector of values for the case of varying residual variance. Default value is 1. 
#'
#' @return An instance of an S4 class of type cpop.class.
#'
#' @references \insertRef{doi:10.1080/10618600.2018.1512868}{cpop}
#'
#' @examples
#'
#' library(cpop)
#' # generate some test data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sigma <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sigma)
#'
#' # use the locations in x
#' res <- cpop(y,x,2*log(length(y)),0.1)
#' plot(res)
#'
#' # without locations (note explicit paramater names)
#' res <- cpop(y,beta=2*log(length(y)),sigsquared=0.1)
#' plot(res)
#'
#' # stretch the end of the data
#' x[75:101] <- x[75:101] + seq(from=0,by=0.2,length.out=27)
#' res <- cpop(y,x,2*log(length(y)),0.1)
#' plot(res)
#'  
#' @export
cpop<-function(y,x=NULL,beta=0,sigsquared=1)
{
    if(is.null(x))
    {
       x<-1:length(y)
    }
    if(length(sigsquared)!=length(y))
    {
      sigsquared=rep(sigsquared[1],length(y))
    }
    res<-CPOP.uneven_impl(y,x,beta,sigsquared,TRUE,FALSE)
    fit<-cpop.fit(y,x,res$changepoints,sigsquared)
    return(cpop.class(y,x,beta,sigsquared,res$min.cost,res$changepoints,fit$fit,fit$residuals,fit$pars))
}

CPOP.uneven_impl<-function(y,x,beta,sigsquared=1,useCprune=FALSE,printiteration=FALSE){
  
  # if(useCprune) dyn.load("prune2R.so")
  
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
  #A<-(SX2[taustar+1]-SX2[sstar+1]-2*Xs*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xs^2)/(seglen^2)
  #B<- 2*( (Xt+Xs)*(SX[taustar+1]-SX[sstar+1])-(SP[taustar+1]-SP[sstar+1])*Xt*Xs-(SX2[taustar+1]-SX2[sstar+1]))/(seglen^2)
  #C<-(-2)/(seglen)*(SXY[taustar+1]-SXY[sstar+1]-Xs*(S[taustar+1]-S[sstar+1]))
  #D<- (SS[taustar+1]-SS[sstar+1])
  #E<-(-2)/(seglen)*(Xt*(S[taustar+1]-S[sstar+1])-(SXY[taustar+1]-SXY[sstar+1]))
  #FF<-(SX2[taustar+1]-SX2[sstar+1]-2*Xt*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xt^2)/(seglen^2)

  coeffs.cpp<-coeffs_update_cpp(SXY,SX2,S,SS,Xs,SX,SP,seglen,taustar,sstar,Xt)
  A<-coeffs.cpp[[1]]
  B<-coeffs.cpp[[2]]
  C<-coeffs.cpp[[3]]
  D<-coeffs.cpp[[4]]
  E<-coeffs.cpp[[5]]
  FF<-coeffs.cpp[[6]]


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




########################################################################################
###################### coeff update#####################################################
###avoids loop
########################################################################################

coeff.update2=function(coeffs,S,SJ,SS,taustar,sigsquared,beta){
  
  coeff.new<-coeffs
  coeff.new[,2]=coeffs[,1] 
  coeff.new[,1]<-taustar
  
  sstar<-coeff.new[,2]
  seglen<-taustar-sstar
  A<-(seglen+1)*(2*seglen+1)/(12*seglen*sigsquared)
  B<- (seglen^2-1)/(6*seglen*sigsquared)
  C<-(-1)/(seglen*sigsquared)*(SJ[taustar+1]-SJ[sstar+1]-sstar*(S[taustar+1]-S[sstar+1]))
  D<-seglen/2*log(2*pi*sigsquared)+1/(2*sigsquared)*(SS[taustar+1]-SS[sstar+1])
  E<-(-1)*C-1/(sigsquared)*(S[taustar+1]-S[sstar+1])
  FF<-(seglen-1)*(2*seglen-1)/(12*seglen*sigsquared)
  
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
    coeff.new[ind2,5]<-(-E[ind2]-coeffs[ind2,4])/B[ind2]+beta
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

########################################################################################
#####coeff.update in C##################
########################################################################################
coeff.update.c<-function(coeffs,S,SJ,SS,taustar,sigsquared,beta){
  nrows<-dim(coeffs)
  coeffs<-as.double(as.vector((t(coeffs))))
  # coutput<- .C("coeffupdate",as.double(t(coeffs)),as.double(S),as.double(SJ),as.double(SS),as.integer(taustar),as.double(sigsquared),as.double(beta),as.integer(nrows[1]),coeffnewlong=double(nrows[1]*nrows[2]))
  coutput<-coeffupdate(t(coeffs),S,SJ,SS,taustar,sigsquared,beta,nrows[1],nrows[2])
  # result<-matrix(coutput$coeffnewlong,nrow=nrows[1],byrow=T)
  result<-matrix(coutput,nrow=nrows[1],byrow=T)
  return(result)
}

############################################################################################################
##first prune of functions
## x is matrix of functions
############################################################################################################
prune1=function(x,taustar){
  min.vals<-x[,5]-(x[,4]^2)/(2*x[,3])
  m=min(min.vals[x[,2]==taustar-1])
  keep<-c(rep(T,length(min.vals[x[,2]!=taustar-1])),(min.vals[x[,2]==taustar-1])<=m)
  return(keep) ##keep only those with a smaller minimum.
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

#############################################################################
#########second pruning in C##################################################
####################################################################################

prune2.c<-function(x){
  nrows<-dim(x)[[1]]
  # coutput<- .C("prune2R",as.double(as.vector((t(x)))),as.integer(nrows),Sets=integer(nrows))
  coutput<-prune(as.vector((t(x))),nrows)
  # result<-which(coutput$Sets==1)
  result<-which(coutput==1)
  return(result)
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


cpop.fit=function(y,x,out.changepoints,sigsquared)
{
  n=length(y)
  if(length(sigsquared)!=n) sigsquared=rep(sigsquared[1],n)
  p=length(out.changepoints)
  W=diag(sigsquared^-1)
  X=matrix(NA,nrow=n,ncol=p)
  X[,1]=1
  X[,2]=x-x[1]
  if(p>2)
  {
    for(i in 2:(p-1))
    {
      X[,i+1]=c(rep(0,out.changepoints[i]-1),x[(out.changepoints[i]):n]-x[out.changepoints[i]])
    }
  }
  XTX=t(X)%*%W%*%X
  beta=as.vector(solve(XTX)%*%t(X)%*%W%*%y)
  fit=X%*%beta
  residuals=y-fit
  return(list(fit=fit,residuals=residuals,X=X,pars=beta))
}