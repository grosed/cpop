
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
#' \code{\link[cpop]{segmentations}} method (see example below).    
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
#' @references \insertRef{cpop-jss-article-2024}{cpop}
#'
#' @export
cpop.crops<-function(y,x = 1:length(y),grid = x, beta_min = 1.5 * log(length(y)),
                     beta_max = 2.5 * log(length(y)),sd = sqrt(mean(diff(diff(y))^2)/6),
		     minseglen = 0,prune.approx = FALSE)
{
#    if(base::missing(sd))
#    {
#       message("No value set for sd. An estimate for the noise standard deviation based on the variance of the second differences of the data has been used. If this estimate is too small it may lead to over-estimation of changepoints. You are advised to check this by comparing the standard deviation of the residuals to the estimated value used for sd.")
#    }
    if(length(sd)!=length(y))
    {
#      message("Length of sd and y differ. Applying first value of sd to all values of y.")
      sd <- rep(sd[1],length(y))
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
                  # locations <- cpts$location
                  cost <- cost(res) - ncpts*beta
                  # return(list(cost,locations,res))
		  return(list(cost,ncpts,res))
               }
              )
    }
    CPT<-CPT.init(y,x,grid,sd,minseglen,prune.approx)
    res.crops <- crops(CPT,beta_min,beta_max)
    return(unique(cpop.crops.class(res.crops,y,x)))
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
#' @seealso \code{\link[cpop]{cpop.crops}},\code{\link[crops]{crops}}
#'
#' @references \insertRef{crops-article}{cpop}
#' @references \insertRef{crops-package}{cpop}
#' @references \insertRef{cpop-jss-article-2024}{cpop}
#'
#' @export
setGeneric("cpop.crops.models",function(object) {standardGeneric("cpop.crops.models")})
setMethod("cpop.crops.models",signature=list("cpop.crops.class"),
function(object)
{
   return(Map(function(beta) as.list(object@method(beta))[[3]] ,as.list(object@betas)))
})



#' Summary of segmentations by penalty value
#'
#' @name segmentations
#'
#' @description Produces a summary of the segmentations for each penalty value in the form of a data frame.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link[crops]{crops}}.
#' 
#' @return A data frame containing the penalties, costs, penalised costs, and changepoint locations. Note - if no changepoints are detected in
#' the penalty interval [beta_min,beta_max], then the value returned is NULL. 
#'
#' @rdname segmentations-methods
#'
#' @aliases segmentations,cpop.crops.class-method
#' 
#' @seealso \code{\link[crops]{crops}}.
#'
#' @examples
#' # see the crops example
setGeneric("segmentations",function(object) {standardGeneric("segmentations")})
setMethod("segmentations",signature=list("cpop.crops.class"),
          function(object)
          {
	    # appease package checks
            . <- NULL
            results <- Map(object@method,unlist(object@betas))
            segs <- Map(function(res) (changepoints(res[[3]]))$location, results) # assumes slot 3 contains a cpop result
	    if(length(Filter(function(seg) length(seg) > 1, segs)) == 0)
	    {
	       return(NULL)
	    }
	    n <- segs %>% Map(length,.) %>% unlist %>% max
            mat <- segs %>%  
                   Map(function(.) c(.,rep(NA,n-length(.))),.) %>%
                   Reduce(rbind,.,matrix(nrow=0,ncol=n),right=TRUE)
                   colnames(mat) <- Map(function(.) paste("cpt.",.,sep=""),1:n) %>% unlist      
            return(      
                    Map(function(beta,seg) tibble(beta=beta,Qm=seg[[1]],Q=seg[[1]]+beta*length(seg),
                        m=length(seg)),
                        unlist(object@betas),
                        segs) %>%                                          
                    Reduce(add_row,.,tibble(beta=numeric(),Qm=numeric(),Q=numeric(),m=numeric())) %>%
                    cbind(.,mat)
                   )           
})


#' Pretty printing for crops results
#'
#' @name print
#'
#' @description Pretty prints a summary of a crops result
#'
#' @docType methods
#'
#' @param  x An instance of an S4 class produced by \code{\link[crops]{crops}}.
#'
#' @rdname print-methods
#'
#' @aliases print,cpop.crops.class-method
#' 
#' @seealso \code{\link[crops]{crops}}.
#'
#' @examples
#' # see the crops example
#'
setMethod("print",signature=list("cpop.crops.class"),
          function(x)
          {
	    summary(x)            
})


#' Visualisation of data, costs, penalty values and changepoint locations.
#'
#' @name plot
#'
#' @description Plot methods for an S4 object returned by \code{\link[crops]{crops}}. The plot can also be combined with the original data if required.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link[crops]{crops}}.
#' @param y A dataframe containing the locations and values of the data points. The data plot is plotted below, and is aligned with, the changepoint locations.
#' 
#' @return A ggplot object. Note - if no changepoints are detected in the penalty interval [beta_min,beta_max], then the value returned is NULL. 
#'
#' @rdname plot-methods
#'
#' @aliases plot,cpop.crops.class,data.frame-method
#' 
#' @seealso \code{\link[crops]{crops}}.
#'
#' @examples
#' # see the crops example
#'
#' @export  
setMethod("plot",signature=list("cpop.crops.class","data.frame"),
          function(x,y)
          {
	    p1 <- plot(x)
	    p2 <- ggplot(data=y) + geom_line(aes(x=x,y=y))
            return(plot_grid(p1,p2,ncol=1))
          })

#' @name plot
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link[crops]{crops}}.
#' 
#' @rdname plot-methods
#'
#' @aliases plot,crops.class,missing-method
#'
#' @examples
#' # see the crops example
#'
#' @export
setMethod("plot",signature=list("crops.class","missing"),
          function(x)
          {
	     # appease ggplot and tidyverse
	     . <- Q <- Qm <- m <- value <- dummy <- NULL
	     object <- x
             df <- segmentations(object)
	     if(is.null(df))
	     {
	        return(NULL)
	     }
	     df <- cbind(df,data.frame("dummy"=1:nrow(df)))
             p <- df %>%
	          subset(.,select = -c(beta,Q,Qm,m)) %>%
                  melt(., id=c("dummy")) %>% 
                  .[complete.cases(.), ] %>%
                  ggplot(.,aes(x=value,y=dummy)) %>% 
                  add(geom_point()) %>%
                  add(labs(x="location",y="penalty")) %>%
                  add(geom_hline(aes(yintercept=dummy))) %>%
                  add(scale_y_continuous(breaks = seq(1:nrow(df)),labels=signif(df$beta,digits=3),sec.axis = sec_axis( ~.,breaks = seq(1:nrow(df)),labels=signif(df$Qm,digits=4),name="unpenalised cost"))) %>%
                  add(theme_bw())
             return(p)       
          })


#' Summary of crops result
#'
#' @name summary
#'
#' @description Prints a short summary of a crops result.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link[crops]{crops}}.
#'
#' @rdname summary-methods
#'
#' @aliases summary,cpop.crops.class-method
#' 
#' @seealso \code{\link[crops]{crops}}.
#'
#' @examples
#' # see the crops example
#'
setMethod("summary",signature=list("cpop.crops.class"),
          function(object)
          {
            cat("crops analysis",sep="")
	    cat('\n',sep="")
	    cat('\n',sep="")
	    cat("minimum penalty value = ",min(object@betas)," : maximum penalty value = ",max(object@betas),sep="")
	    cat('\n',sep="")
	    segs <- segmentations(object)
	    if(is.null(segs))
	    {
	       cat("no segmentations found in the penalty interval [",min(object@betas),",",max(object@betas),"]",sep="")
	       cat('\n',sep="")
	    }
	    else
	    {
	       cat("number of segmentations calculated : ",nrow(segs),sep="")
	       cat('\n',sep="")	    
	       cat("least number of changepoints  = ",min(segs$m), " : maximum number of changepoints = ",max(segs$m),sep="")
	       cat('\n',sep="")
	    }
	  invisible()
})


#' Remove duplicate entries from a crops result
#'
#' @name unique
#'
#' @description Removes duplicate entries from a crops result. A duplicate entry is one having the same number of changepoints as another entry.
#' Note that the changepoint locations and the associated penalty and cost values are not taken into consideration. The \code{unique} function can be useful
#' for simplifying plots and the details produced by \code{segmentations}.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link[crops]{crops}}.
#'
#' @return An instance of the S4 class type \code{crops.class}. This is the same type as produced by the \code{\link[crops]{crops}} function.
#'
#' @rdname unique-methods
#'
#' @aliases unique,cpop.crops.class-method
#' 
#'
setMethod("unique",signature=list("cpop.crops.class"),
         function(x)
             {
	     	# appease package checks
                . <- NULL
                object<-x
                hash_map <- new.env()
                keys <- object@betas %>% 
                unlist %>% 
                Map(object@method,.) %>% 
                Map(function(.) .[2],.) %>% 
                Map(as.character,.)
                key_value_pairs <- Map(tuple,keys,object@betas %>% unlist)
                hash_map <-    
                key_value_pairs %>%  
                Reduce(function(pair,map) {map[[pair[[1]]]] <- pair[[2]]; return(map);},
                   .,
                   hash_map,right=TRUE)
                object@betas <-    
                Map(function(key) hash_map[[key]],
                    hash_map %>% ls) %>%
                unname %>% 
                as.set
                return(object)
             })


#' Subset crops results based on penalty values
#'
#' @name subset
#'
#' @description Removes entries from a crops result that fall outside a specified range of penalty values. 
#' The \code{subset} function can be useful for simplifying plots and the details produced by \code{segmentations}.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link[crops]{crops}}.
#' @param beta_min A positive numeric value specifying the minimum penalty value for entries in the crops result. Default value is 0.
#' @param beta_max A positive numeric value specifying the maximum penalty value for entries in the crops result. Default value is Inf.
#'
#' @return An instance of the S4 class type \code{crops.class}. This is the same type as produced by the \code{\link[crops]{crops}} function.
#'
#' @rdname subset-methods
#'
#' @aliases subset,crops.class-method
#' 
#'
setMethod("subset",signature=list("crops.class"),
         function(x,beta_min=0,beta_max=Inf)
             {
	     	# appease package checks
                . <- NULL
	        object <- x
                object@betas %<>% 
                unlist %>% 
                Filter(function(.) . <= beta_max & . >= beta_min,.) %>% 
                as.set
                return(object)            
             })

