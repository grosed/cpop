
context("testing simchangeslope")

test_that("test that simchangeslope produces data consistently",
{
   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   sd=1+x/200
   y<-simchangeslope(x,changepoints,change.slope,sd)
   load("test.RData")
   expect_equal(y,test.simulate.y)
})