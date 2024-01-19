test_that("new_mixture returns expected error when dist and pdf_func are not given", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.2,0.2)
  
  params = c(eta = p, mu = mu, sigma = sigma)
  
  expect_error(new_Mixture(params),  "you have to specify either dist or pdf_func")
})