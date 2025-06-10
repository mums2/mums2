test_that("mums2_example returns proper paths to test data", {
  path <- mums2_example("")
  expect_true("character" %in% class(mums2_example("squid_peak_table.csv")))
  expect_true("character" %in% class(mums2_example("squid_meta_data.csv")))
  expect_true("character" %in% class(mums2_example()))
  expect_true(length(mums2_example()) == 7)
  expect_error(mums2_example("dat"))
})
