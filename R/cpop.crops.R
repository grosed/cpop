
.cpop.crops.class<-setClass("cpop.crops.class",contains="crops.class",representation(y="numeric",x="numeric",sd="numeric"))

cpop.crops.class<-function(crops.result,y,x)
{
   .cpop.crops.class(crops.result,y=y,x=x)    
}

#' Calculate changepoint locations over a range of penalty values
#'
#' @name cpop.crops
#'
#' @description Runs the Changepoints for a Range of Penalties (CROPS) algorithm of Haynes et al. (2014) to find all of the optimal segmentations for multiple penalty values over a continuous range.
#' For details of the \code{crops}function see the \pkg{crops package}. Methods from the \pkg{crops} package that are applicable to crops are
#' \code{\link[crops]{segmentations}}, \code{\link[crops]{subset}}, and \code{\link[crops]{unique}}.    
#'
#' @param y A vector of length n containing the data.
#' @param x A vector of length n containing the locations of y. Default value is NULL, in which case the locations \code{x = 1:length(y)-1} are assumed.
#' @param grid An ordered vector of possible locations for the change points. If this is NULL, then this is set to x, the vector of times/locations of the data points.
#' @param beta_min Minimum value of the penalty range to be searched. Default is \code{1.5*log(length(y))}
#' @param beta_max Maximum value of the penalty range to be searched. Default is \code{2.5*log(length(y))} 
#' @param sd Estimate of residual standard deviation. Can be a single numerical value or a vector of values for the case of varying standard deviation. Default value is 1.
#' @param minseglen The minimum allowable segment length, i.e. distance between successive changepoints. Default is 0.
#' @param prune.approx Only relevant if a minimum segment length is set. If True, cpop will use an approximate pruning algorithm that will speed up computation but may
#' occasionally lead to a sub-optimal solution in terms of the estimate change point locations. If the minimum segment length is 0, then an exact pruning algorithm is possible and is used.
#'
#' @return An instance of an S4 class of type cpop.crops.class which extends the crops.class in \pkg{crops}.
#'
#' @rdname cpop.crops-methods
#'
#' @aliases crops.cpop,cpop.crops.class-method
#'
#' @examples
#'
#' # generate some data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sigma <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sigma)
#' # calculate the changepoint locations and cost over a range of penalty values
#' res <- cpop.crops(y,x,sd=sigma,beta_min=0.4*log(length(y)),beta_max=2.5*log(length(y)))
#' summary(res)
#'
#' # plot the results
#' plot(res)
#' plot(unique(res)) 
#'
#' # show the segmentations
#  segmentations(res)
#'
#' @references \insertRef{crops-article}{crops}
#' @references \insertRef{crops-package}{cpop}
#' @export
cpop.crops<-function(y,x = 1:length(y),grid = x, beta_min = 1.5 * log(length(y)),
                     beta_max = 2.5 * log(length(y)),sd = 1,minseglen = 0,prune.approx = FALSE)
{
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
                  return(list(cost,locations))
               }
              )
    }
    CPT<-CPT.init(y,x,grid,sd,minseglen,prune.approx)
    res.crops <- crops(CPT,beta_min,beta_max)
    return(cpop.crops.class(res.crops,y,x))
}