
test_that("result of shift", {
  set.seed(100)
  data1 <- rnorm(n = 300, mean = 3, sd = 0.5)
  set.seed(400)
  data2 <- rnorm(n = 300, mean = 3, sd = 0.5)
  results <- BNP.test(x=data1,y=data2, n.mcm=1000)
  expect_s3_class(plotshift.function(results_BNP = results), "plotly")
})


