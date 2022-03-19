

test_that("Final list of test results", {
  data(airquality)
  results <- BNP.test(x=airquality[airquality$Month==5,3],y=airquality[airquality$Month==8,3], n.mcm=1000)
  expect_type(results, "list")
})


test_that("Length of paired data", {
  data.pre <- rnorm(10)
  data.post <- rnorm(15)
  expect_error(BNP.test(x=data.pre,y=data.post, n.mcm=1000), "the length of the vectors x, y is not equal")
})



test_that("Length of results", {
  data.pre <- rnorm(10)
  data.post <- rnorm(10)
  results <- BNP.test(x=data.pre,y=data.post, n.mcm=1000)
  expect(length(results), 3)
})



test_that("Length of sampling paramethers", {
  data.pre <- rnorm(10)
  data.post <- rnorm(10)
  results <- BNP.test(x=data.pre,y=data.post, n.mcm=1000)
  expect(length(results$sampling.parameters),100 )
})




test_that("value of probability", {
  set.seed(100)
  data1 <- rnorm(n = 300, mean = 3, sd = 0.5)
  set.seed(400)
  data2 <- rnorm(n = 300, mean = 3, sd = 0.5)
  results <- BNP.test(x=data1,y=data2, n.mcm=1000)
  expect_true(results$posterior.probability.H1>=0)
  expect_true(results$posterior.probability.H1<=1)
})


