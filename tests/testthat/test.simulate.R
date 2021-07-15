
context("testing simulate")

test_that("test that simulate produces data consistently",
{
   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   sig=1+x/200
   y<-simulate(x,changepoints,change.slope,sig)
   load("test.RData")
   expect_equal(y,test.simulate.y)
})