
.cpop.crops.class<-setClass("cpop.crops.class",contains="crops.class",representation(y="numeric",x="numeric",sd="numeric"))

cpop.crops.class<-function(crops.result,y,x)
{
   .cpop.crops.class(crops.result,y=y,x=x)    
}

#' Calculate changepoint locations over a range of penalty values
#'
#' @name cpop.crops
#'
#' @description Runs the Changepoints for a Range of Penalties (CROPS) algorithm of Haynes et al. (2017) to find all of the optimal segmentations for multiple penalty values over a continuous range.
#' For details of the CROPS method see the \pkg{\link[crops]{crops}} package. To obtain details for each segmentation determined by \code{\link[cpop]{cpop.crops}} use the 
#' \code{\link[crops]{segmentations}} method (see example below).    
#'
#' @param y A vector of length n containing the data.
#' @param x A vector of length n containing the locations of y. Default value is NULL, in which case the locations \code{x = 1:length(y)-1} are assumed.
#' @param grid An ordered vector of possible locations for the change points. If this is NULL, then this is set to x, the vector of times/locations of the data points.
#' @param beta_min Minimum value of the penalty range to be searched. Default is \code{1.5*log(length(y))}
#' @param beta_max Maximum value of the penalty range to be searched. Default is \code{2.5*log(length(y))} 
#' @param sd Estimate of residual standard deviation. Can be a single numerical value or a vector of values for the case of varying standard deviation.
#' Default value is \code{sd = sqrt(mean(diff(diff(y))^2)/6)}.
#' @param minseglen The minimum allowable segment length, i.e. distance between successive changepoints. Default is 0.
#' @param prune.approx Only relevant if a minimum segment length is set. If TRUE, cpop will use an approximate pruning algorithm that will speed up computation but may
#' occasionally lead to a sub-optimal solution in terms of the estimate change point locations. If the minimum segment length is 0, then an exact pruning algorithm is possible and is used. 
#'
#' @return An instance of an S4 class of type cpop.crops.class which extends the crops.class in \pkg{crops}.
#'
#' @rdname cpop.crops-methods
#'
#' @aliases crops.cpop,cpop.crops.class-method
#'
#' @examples
#' \donttest{
#' library(cpop)
#'
#' # generate some data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sd <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sd)
#' # calculate the changepoint locations and cost over a range of penalty values
#' res <- cpop.crops(y,x,sd=sd,beta_min=0.4*log(length(y)),beta_max=2.5*log(length(y)))
#' summary(res)
#'
#' # plot the results
#' plot(res) 
#'
#' # show the segmentations
#'  segmentations(res)
#' }
#'
#' @references \insertRef{crops-article}{cpop}
#' @references \insertRef{crops-package}{cpop}
#' @export
cpop.crops<-function(y,x = 1:length(y),grid = x, beta_min = 1.5 * log(length(y)),
                     beta_max = 2.5 * log(length(y)),sd = sqrt(mean(diff(diff(y))^2)/6),
		     minseglen = 0,prune.approx = FALSE)
{
    if(base::missing(sd))
    {
       cat("No value set for sd. An estimate for the noise standard deviation based on the variance of the second differences of","\n",
           "the data has been used. If this estimate is too small it may lead to over-estimation of changepoints. You are advised","\n",
	   "to check this by comparing the standard deviation of the residuals to the estimated value used for sd.","\n",sep="")
    }
    CPT.init <- function(y,x,grid,sd,minseglen,prune.approx)
    {
       return(
           function(beta)
               {
                  res <- cpop(y=y,x=x,grid=grid,beta=beta,sd=sd,minseglen=minseglen,
                              prune.approx=prune.approx)
                  cpts <- changepoints(res)
                  ncpts <- nrow(cpts)
                  locations <- cpts$location
                  cost <- cost(res) - ncpts*beta
                  return(list(cost,locations,res))
               }
              )
    }
    CPT<-CPT.init(y,x,grid,sd,minseglen,prune.approx)
    res.crops <- crops(CPT,beta_min,beta_max)
    res.crops <- unique(res.crops)
    return(cpop.crops.class(res.crops,y,x))
}


#' Extract the cpop models created by cpop.crops
#'
#' @name cpop.crops.models
#'
#' @description Obtains a list of models corresponding to each of the beta (penalty) values considered during a \code{cpop.crops} analysis. 
#'
#' @param object An S4 object of type \code{cpop.crops.class} produced by \code{cpop.crops}.
#'
#' @return A list of S4 \code{cpop.class} objects corresponding to the beta values considered by a \code{cpop.crops} analysis.
#'
#' @rdname cpop.crops.models-methods
#'
#' @aliases cpop.crops.models,cpop.crops.class-method
#'
#' @examples
#' \donttest{
#' library(cpop)
#'
#' set.seed(1)
#' n <- 500
#' x <- 1:n
#' m <- 10
#' mu <- simchangeslope(x,changepoints=(n/(m+1))*0:m,change.slope=c(0.1,0.2*(-1)^(1:m)),sd=0)
#' epsilon <- rnorm(n+2)
#' y <- mu+(epsilon[1:n]+epsilon[2:(n+1)]+epsilon[3:(n+2)])/sqrt(3)
#' res.crops <- cpop.crops(y,x,beta_min=0.5*log(length(y)),beta_max=40*log(length(y)))
#' models <- cpop.crops.models(res.crops)
#' for(m in models)
#' {
#'   plot(m)
#' }
#'}
#'
#' @seealso \code{\link{cpop.crops}},\code{\link[crops]{crops}}
#'
#' @references \insertRef{crops-article}{cpop}
#' @references \insertRef{crops-package}{cpop}
#' @export
setGeneric("cpop.crops.models",function(object) {standardGeneric("cpop.crops.models")})
setMethod("cpop.crops.models",signature=list("cpop.crops.class"),
function(object)
{
   return(Map(function(beta) as.list(object@method(beta))[[3]] ,as.list(object@betas)))
})
