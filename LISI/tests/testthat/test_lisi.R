library(testthat)
library(lisi)
context('Test basic LISI functions')

lisi_res <- lisi::compute_lisi(lisi::X, lisi::meta_data, c('label1', 'label2'))

test_that('LISI values are between 1 and total number of labels', {
  expect_gte(min(lisi_res), 1)
  expect_lte(max(lisi_res), 2)
})

test_that('dimensions are correct', {
  expect_equal(dim(lisi_res), dim(lisi::meta_data))
})

test_that('Label 2 gives higher mean LISI than label 1', {
  expect_gt(mean(lisi_res$label2), mean(lisi_res$label1))
})
