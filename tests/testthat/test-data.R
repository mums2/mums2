test_that("example returns proper paths to test data", {
  path <- example("")
  expect_true("character" %in% class(example("squid_peak_table.csv")))
  expect_true("character" %in% class(example("squid_meta_data.csv")))
  expect_true("character" %in% class(example()))
  expect_true(length(example()) == 7)
  expect_error(example("dat"))
})
