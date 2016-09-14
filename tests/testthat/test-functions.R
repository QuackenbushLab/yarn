################################################################################
# yarn function unit tests
################################################################################
library("yarn"); library("testthat"); 

test_that("`filterLowGenes` function provides expected values", {
  data(skin)
  expect_equal(as.numeric(nrow(skin)),40824)
  skin = filterLowGenes(skin,"SMTSD")
  expect_equal(as.numeric(nrow(skin)),19933)
})