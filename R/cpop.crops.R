
.cpop.crops.class<-setClass("cpop.crops.class",contains="crops.class",representation(y="numeric",x="numeric",sd="numeric"))

cpop.crops.class<-function(crops.result,y,x)
{
   .cpop.crops.class(crops.result,y=y,x=x)    
}

setMethod("plot",signature=list("cpop.crops.class","missing"),
          function(x)
          {
            object <- as(x,"crops.class")
            return(plot(object,data.frame("x"=x@x,"y"=x@y)))
          })

cpop.crops<-function(y,x = 1:length(y),grid = x, beta_min = 1.5 * log(length(y)),
                     beta_max = 2.5 * log(length(y)),sd = 1,minseglen = 0,prune.approx = FALSE)
{
    library(crops)
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
                  return(list(cost,ncpts,locations))
               }
              )
    }
    CPT<-CPT.init(y,x,grid,sd,minseglen,prune.approx)
    res.crops <- crops(CPT,beta_min,beta_max)
    return(cpop.crops.class(res.crops,y,x))
}