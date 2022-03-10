
.cpop.class<-setClass("cpop.class",representation(changepoints="numeric",
						  y="numeric",
						  x="numeric",
						  beta="numeric",
						  sd="numeric"))

cpop.class<-function(y,x,beta,sd,changepoints)
{
    .cpop.class(y=y,
                x=x,
		changepoints=changepoints,
		beta=beta,
		sd=sd)
}



#' A function for calculating the cost of a model fitted by cpop
#'
#' @name cost
#'
#' @description Calculates the penalised cost of a model fitted by cpop using the residual sum of squares and the penalty values.
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @return Numerical value of the penalised cost associated with the segmentations determined by using \code{\link{cpop}}   
#'

#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # determine changepoints
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # calculate the penalised cost 
#' cost(res)
#'
#' @rdname cost-methods
#'
#' @aliases cost,cpop.class-method
#'
#' @export
setGeneric("cost",function(object) {standardGeneric("cost")})
setMethod("cost",signature=list("cpop.class"),
          function(object)
          {
            df <- fitted(object)
	    return(sum(residuals(object)^2/object@sd^2)+object@beta*(length(object@changepoints)-2))
          })	      


#' A function for generating simulated data
#'
#' @name simulate
#'
#' @description Generates simulated data for use with \code{\link{cpop}}.
#' 
#' @param x A numeric vector containing the locations of the data.
#' @param changepoints A numeric vector of changepoint locations.
#' @param change.slope A numeric vector indicating the change in slope at each changepoint. The initial slope is assumed to be 0.
#' @param sigma The residual standard deviation. Can be a single numerical value or a vector of values for the case of varying residual standard deviation. Default value is 1.
#' @return A vector of simulated y values.
#'
#' @examples
#' library(cpop)
#' library(pacman)
#' p_load(ggplot2)
#' 
#' # simulate changepoints with constsant sigma
#' set.seed(1)
#' changepoints <- c(0,25,50,100)
#' change.slope <- c(0.2,-0.3,0.2,-0.1)
#' x <- 1:200
#' sig <- 0.2
#' y <- simulate(x,changepoints,change.slope,sig)
#' df <- data.frame("x"=x,"y"=y)
#' p <- ggplot(data=df,aes(x=x,y=y))
#' p <- p + geom_point()
#' p <- p + geom_vline(xintercept=changepoints,
#'                     color="red",
#'                     linetype="dashed")
#' print(p)
#' 
#' # simulate changepoints with varying sigma
#' sig <- 0.2 + x/200
#' y <- simulate(x,changepoints,change.slope,sig)
#' df$y <- y
#' p <- ggplot(data=df,aes(x=x,y=y))
#' p <- p + geom_point()
#' p <- p + geom_vline(xintercept=changepoints,
#'                     color="red",
#'                     linetype="dashed")
#' print(p)
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

#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # analyse data
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # generate plot object
#' p <- plot(res)
#'
#' # visualise
#' print(p)
#'
#' @rdname plot-methods
#'
#' @aliases plot,cpop.class,ANY-method
#' 
#' @export 
setMethod("plot",signature=list("cpop.class"),function(x)
{
  object <- x
  df <- data.frame("x"=object@x,"y"=object@y)
  cpts<-object@changepoints
  # appease ggplot2
  y <- NULL
  p <- ggplot(data=df, aes(x=x, y=y))
  p <- p + geom_point(alpha=0.3)
  if(length(cpts) > 2)
  {
   for(cpt in cpts[2:(length(cpts)-1)])
   {
     p <- p + geom_vline(xintercept = cpt,color="red")
   }
  }
  # appease ggplot2
  x0 <- y0 <- x1 <- y1 <- NULL  
  #df<-data.frame("xs"=obj@x[cpts][1:(length(obj@x[cpts])-1)],
  #	         "ys"=obj@y_hat[cpts][1:(length(obj@y_hat[cpts])-1)],
  #		 "xends"=obj@x[cpts][2:length(obj@x[cpts])],
  #		 "yends"=obj@y_hat[cpts][2:length(obj@y_hat[cpts])])
  p <- p + geom_segment(data=fitted(object),aes(x=x0,y=y0,xend=x1,yend=y1))
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
#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # determine changepoints
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # display a summary of the results
#' summary(res)
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
     msg<-paste(msg," detected at x = ",'\n',sep="")
     cat(msg)
     msg<-""
     for(cpt in object@changepoints[2:(length(object@changepoints)-1)])
     {
       msg<-paste(msg,cpt,sep=" ")
     }
     msg<-paste(msg,'\n',sep="")
     cat(msg)
  }
  df <- fitted(object)
  cat("fitted values : ",'\n',sep="")
  print(df)
  cat('\n',"overall RSS = ",sum(df$RSS),'\n',sep="")
  cat("cost = ",cost(object),'\n\n',sep="")
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
#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # determine changepoints
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # display a summary of the results using show
#' show(res)
#'
#' @export
setGeneric("show",function(object) {standardGeneric("show")})
setMethod("show",signature=list("cpop.class"),function(object)
{
    summary(object)
    invisible()
})


#' Extract Model Fitted Values
#'
#' @name fitted
#'
#' @description Extracts the fitted values produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @return A data frame containing the endpoint coordinates for each line segment fitted between
#' the detected changepoints. The data frame also contains the gradient and intercept values
#' for each segment and the corresponding residual sum of squares (RSS). 
#'
#' @rdname fitted-methods
#'
#' @aliases fitted,cpop.class-method
#'
#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # determine changepoints
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # calculate segments
#' fitted(res)

#' @export
setGeneric("fitted",function(object) {standardGeneric("fitted")})
setMethod("fitted",signature=list("cpop.class"),
          function(object)
          {
	  x<-sort(unique(c(object@x,object@changepoints)))
	  y_hat<-estimate(object,x)$y_hat
	  cpts<-object@changepoints
	  y_0<-unlist(Map(function(val) y_hat[which(x==val)],cpts))
  	  df<-data.frame("x0"=cpts[1:(length(cpts)-1)],
			 "y0"=y_0[1:(length(cpts)-1)],
		         "x1"=cpts[2:length(cpts)],
		         "y1"=y_0[2:length(y_0)])
          df <- cbind(df,data.frame("gradient"=(df$y1 - df$y0)/(df$x1 - df$x0)))
	  df <- cbind(df,data.frame("intercept"=df$y1 - df$gradient * df$x1))
	  residuals<-residuals(object)
          rss<-unlist(Map(function(a,b) sum((residuals*residuals)[which(object@x >= a & object@x < b)]),cpts[1:(length(cpts)-1)],cpts[2:length(cpts)]))
	  rss[length(rss)]<-rss[length(rss)]+(residuals*residuals)[length(residuals)]
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
#' @param object  An instance of an cpop S4 class produced by \code{\link{cpop}}.
#' 
#' @return A data frame.
#' 
#' @aliases changepoints,cpop.class-method
#'
#' @examples
#' library(cpop)
#'
#' # generate some test data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sigma <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sigma)
#'
#' # use the locations in x
#' res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
#' changepoints(res)
#'
#' @export
if(!isGeneric("changepoints")) {setGeneric("changepoints",function(object) {standardGeneric("changepoints")})}
setGeneric("changepoints",function(object) {standardGeneric("changepoints")})
setMethod("changepoints",signature=list("cpop.class"),
          function(object)
          {
	      if(length(object@changepoints) > 2)
	      {
	      	      df <- data.frame("location"=object@changepoints[2:(length(object@changepoints)-1)])	
	      }
	      else
	      {
	      	      df <- data.frame("location"=numeric(0))	
	      }
	      return(df)
          })	      
#' cpop
#'
#'  Algorithm for finding the best segmentation of data for a change-in-slope model.
#' 
#' @param y A vector of length n containing the data.
#' @param x A vector of length n containing the locations of y. Default value is NULL, in which case the locations \code{x = 1:length(y)-1} are assumed.
#' @param grid An ordered vector of possible locations for the change points. If this is NULL, then this is set to x, the vector of times/locations of the data points.
#' @param beta A positive real value for the penalty incurred for adding a changepoint (prevents over-fitting).
#' @param sd Estimate of residual standard deviation. Can be a single numerical value or a vector of values for the case of varying standard deviation. Default value is 1.
#' @param minseglen The minimum allowable segment length, i.e. distance between successive changepoints. Default is 0.
#' @param prune.approx Only relevant if a minimum segment length is set. If True, cpop will use an approximate pruning algorithm that will speed up computation but may
#' occasionally lead to a sub-optimal solution in terms of the estimate change point locations. If the minimum segment length is 0, then an exact pruning algorithm is possible and is used.
#'
#' @return An instance of an S4 class of type cpop.class.
#'
#' @description \loadmathjax{} The CPOP algorithm fits a change-in-slope model to data. It assumes that we have we have data points \mjeqn{(y_1,x_1),\ldots,(y_n,x_n)}{(y_1,x_1),...,(y_n,x_n)}, ordered
#' so that \mjeqn{x_1 < x_2 < \cdots < x_n}{x_1 < x_2 < ... < x_n}. For example  \mjeqn{x_i}{x_i} could be a time-stamp of when response \mjeqn{y_i}{y_i} is obtained. We model the response, \mjeqn{y}{y}, as a signal
#' plus noise where the signal is modelled as a continuous piecewise linear function of \mjeqn{x}{x}. That is \mjdeqn{y_i=f(x_i)+\epsilon_i}{y_i=f(x_i)+e_i} where \mjeqn{f(x)}{f(x)} is a continuous
#' piecewise linear function.
#'
#' To estimate the function \mjeqn{f(x)}{f(x)} we specify a set of \mjeqn{N}{N} grid points, \mjeqn{g_{1:N}}{g_{1:N}} with these ordered so that \mjeqn{g_i < g_j}{g_i < g_j} if and only if \mjeqn{i < j}{i < j},
#' and allow the slope of \mjeqn{f(x)}{f(x)} to only change at these grid points. We then estimate the number of changes, the location of the changes, and hence the resulting
#' function \mjeqn{f(x)}{f(x)} by minimising a penalised weighted least squares criteria. This criteria includes a term which measures fit to the data, which is calculated as the sum of the
#'square residuals of the fitted function, scaled by the variance of the noise. The penalty is proportional to the number of changes.
#'
#' That is our estimated function will depend on \mjeqn{K}{K}, the number of changes in slope, their locations, \mjeqn{\tau_1,\ldots,\tau_K}{tau_1,...,tau_K}, and the value of the function
#' \mjeqn{f(x)}{f(x)} at these change points, \mjeqn{\alpha_1,\ldots,\alpha_K}{alpha_1,...,alpha_K}, and its values, \mjeqn{\alpha_0}{alpha_0} at  \mjeqn{\tau_0 < x_1}{tau_0 < x_1} and
#' \mjeqn{\alpha_{K+1}}{alpha_{K+1}} at some \mjeqn{\tau_{K+1} > x_N}{tau_{K+1}>x_N}. The CPOP algorithm then estimates \mjeqn{K}{K}, \mjeqn{\tau_{1:K}}{tau_{1:K}} and \mjeqn{\alpha_{0:K+1}}{alpha_{0:K+1}} 
#' as the values that solve the following minimisation problem
#'\mjdeqn{\min_{K,\tau_{1:K}\in g_{1:N}, \alpha_{0:K+1} }\left\lbrace\sum_{i=1}^n \frac{1}{\sigma^2_i} \left(y_i -  \alpha_{j(i)}-(\alpha_{j(i)+1}-  \alpha_{j(i)})\frac{x_i-\tau_{j(i)}}{\tau_{j(i)+1}-\tau_{j(i)}}\right)^2+K\beta\right\rbrace}{\min_{K,\tau_{1:K}\in g_{1:N}, \alpha_{0:K+1} }\left\lbrace\sum_{i=1}^n \frac{1}{\sigma^2_i} \left(y_i -  \alpha_{j(i)}-(\alpha_{j(i)+1}-  \alpha_{j(i)})\frac{x_i-\tau_{j(i)}}{\tau_{j(i)+1}-\tau_{j(i)}}\right)^2+K\beta\right\rbrace}
#'
#' where \mjeqn{\sigma^2_1,\ldots,\sigma^2_n}{sigma^2_1,...,sigma^2_n} are the variances of the noise \mjeqn{\epsilon_i}{\epsilon_i} for \mjeqn{i=1,\ldots,n}{i=1,...,n}, and \mjeqn{\beta}{beta} is the
#' penalty for adding a changepoint. The sum in this expression is the weighted residual sum of squares, and the \mjeqn{K\beta}{K*beta} term is the penalty for having \mjeqn{K}{K} changes. 
#'
#' If we know, or have a good estimate of, the residual variances, and the noise is (close to) independent over time then an appropriate choice for the penalty is
#' \mjeqn{\beta=2 \log n}{beta=2log(n)}, and this is the default for CPOP. However in many applications these assumptions will not hold and it is advised to look at segmentations for
#' different value of \mjeqn{\beta}{beta} -- this is possible
#' using CPOP with the CROPS algorithm \code{\link{cpop.crops}}. Larger values of \mjeqn{\beta}{beta} will lead to functions with fewer changes. Also there is a trade-off between the variances of the residuals
#' and \mjeqn{\beta}{beta}: e.g. if we double the variances and half the value of \mjeqn{\beta}{beta} we will obtain the same estimates for the number and location of the changes and the
#' underlying function.
#'
#' @references \insertRef{doi:10.1080/10618600.2018.1512868}{cpop}
#'
#' @examples
#' library(cpop)
#' 
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' # analyse using cpop
#' res <- cpop(y,x)
#' p <- plot(res)
#' print(p)
#'  
#' # generate the "true" mean
#' mu <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=0)
#' # add the true mean to the plot
#' library(pacman)
#' p_load(ggplot2)
#' p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dotted")
#' print(p)
#'
#' # heterogeneous data
#' set.seed(1)
#' sigma <- (1:50)/25
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=sigma)
#' 
#' # analysis assuming constant noise standard deviation
#' res <- cpop(y,x,beta=2*log(length(y)),sd=sqrt(mean(sigma^2)))
#' p <- plot(res)
#' print(p)
#'
#' # analysis with the true noise standard deviation
#' res.true <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
#' p <- plot(res.true)
#' print(p)
#' 
#' # add the true mean to the plot
#' p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dotted")
#' print(p)
#'
#' @export
cpop<-function(y,x=1:length(y)-1,grid=x,beta=2*log(length(y)),sd=1,minseglen=0,prune.approx=FALSE)
{
    if(length(sd)!=length(y))
    {
      sd=rep(sd[1],length(y))
    }
    sigsquared<-sd^2
    if(minseglen != 0)
    {
      res<-cpop.grid.minseglen(y,x,grid,beta,sigsquared,minseg=minseglen,FALSE,prune.approx)
    }
    else
    {
      res<-cpop.grid(y,x,grid,beta,sigsquared)
    }
    return(cpop.class(y,x,beta,sd,res$changepoints))
}


design<-function(object,x=object@x)
{
  n=length(x)
  cpts<-object@changepoints[-1]
  p=length(cpts)
  X=matrix(NA,nrow=n,ncol=p+1)
  X[,1]=1
  X[,2]=x-x[1]
  if(p>1)
  {
    for(i in 1:(p-1))
    {
      X[,i+2]=pmax(rep(0,n),x-cpts[i])
    }
  }
  return(X)
}

parameters<-function(object)
{
  n=length(object@y)
  cpts<-object@changepoints[-1]	
  p<-length(cpts)
  W<-diag(object@sd^-2)
  X<-design(object)
  XTX<-t(X)%*%W%*%X
  pars<-as.vector(solve(XTX)%*%t(X)%*%W%*%object@y)
  return(pars)
}

#' A function for estimating the fit of a cpop model
#'
#' @name estimate
#'
#' @description Estimates the fit of a cpop model at the specified locations. If no locations are specified it evaluates the estimates
#' at the locations specified when calling \code{\link{cpop}}. 
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#' @param x Locations at which the fit is to be estimated. Default value is the x locations at which the cpop object was defined.
#' @param ... Additional arguments.
#'
#' @return A data frame with two columns containing the locations x and the corresponding estimates y_hat. 
#'
#' @rdname estimate-methods
#'
#' @aliases estimate,cpop.class-method
#'
#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # determine changepoints
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # estimate fit at points used in call to cpop
#' estimate(res)
#'
#' # estimate fit at specified locations
#' estimate(res,seq(0,100,10))
#' 
#' #extrapolate fit
#' estimate(res,seq(-20,140,20))
#' 
#' @export
setGeneric("estimate",function(object,x=object@x,...) {standardGeneric("estimate")})
setMethod("estimate",signature=list("cpop.class"),
          function(object,x)
          {
             return(data.frame("x"=x,"y_hat"=design(object,x)%*%parameters(object)))
          })	      


#' Extracts residuals from a cpop model
#'
#' @name residuals
#'
#' @description Extracts residuals from a cpop model.
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @return A single column matrix containing the residuals. 
#'
#' @rdname residuals-methods 
#'
#' @examples
#' library(cpop)
#'
#' # simulate data with change in gradient
#' set.seed(1)
#' x <- (1:50/5)^2
#' y <- simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sigma=1)
#' 
#' # determine changepoints
#' res <- cpop(y,x,beta=2*log(length(y)))
#'
#' # calculate the residuals
#' residuals(res)
#'
#' @aliases residuals,cpop.class-method
#'
#' @export
setMethod("residuals",signature=list("cpop.class"),
          function(object)
          {
	    object@y-design(object)%*%parameters(object)
          })	      


#residuals<-function(object)
#{
#  object@y-design(object)%*%parameters(object)
#}


cpop.fit<-function(y,x,out.changepoints,sigsquared)
{
  n=length(y)
  if(length(sigsquared)!=n)
  {
     sigsquared=rep(sigsquared[1],n)
  }
  p=length(out.changepoints)
  W=diag(sigsquared^-1)
  X=matrix(NA,nrow=n,ncol=p+1)
  X[,1]=1
  X[,2]=x-x[1]
  if(p>1)
  {
    for(i in 1:(p-1))
    {
      X[,i+2]=pmax(rep(0,n),x-out.changepoints[i])

    }
  }
  XTX=t(X)%*%W%*%X
  beta=as.vector(solve(XTX)%*%t(X)%*%W%*%y)
  fit=X%*%beta
  residuals=y-fit
  return(list(fit=fit,residuals=residuals,X=X,pars=beta))
}



CPOP.uneven_impl<-function(y,x,beta,sigsquared=1)
{  
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
       keep2=prune2.c(new.coeffs.p)
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




##########################################################################################################
##second pruning
## again x is matrix of quadratics
##version to avoid nested loops
##########################################################################################################





###########################################################################################################
###################################PELT pruning function###################################################
###########################################################################################################

peltprune=function(x,beta){
  minx<-x[,5]-x[,4]^2/(4*x[,3])
  return(which(minx<=(min(minx)+2*beta)))  
}


prune2.c<-function(x){
  nrows<-dim(x)[[1]]
  # coutput<- .C("prune2R",as.double(as.vector((t(x)))),as.integer(nrows),Sets=integer(nrows))
  coutput<-prune(as.vector((t(x))),nrows)
  # result<-which(coutput$Sets==1)
  result<-which(coutput==1)
  return(result)
}

################end#######################




