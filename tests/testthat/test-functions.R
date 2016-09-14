################################################################################
# yarn function unit tests
################################################################################
library("yarn"); library("testthat"); 

test_that("`calcNormFactors` function provides expected values", {
  data(skin)
  expect_equal(nrow(skin),40824)
  skin = filterLowGenes(skin)
  expect_equal(nrow(skin),19933)
})