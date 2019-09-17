library(testthat)

context("Test test for now.")

test_that("The test test works.",{
  expect_equivalent(0, 0, tolerance=1e-5)
})