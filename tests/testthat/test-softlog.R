# test_softlog.R

# Source the softlog function
source("R/softlog.R")

# Load the testthat library
library(testthat)

test_that("softlog works", {
  expect_equal(c(0,log(1/2),log(1e-20)), softlog(c(1,1/2,0)))
})

test_that("softlog catches errors", {
  expect_error(softlog(c(1,0,-1)), "Error: softlog attempt to compute log of a negative value.")
})
