test_that("mums2_example returns proper paths to test data", {
  path <- mums2_example("")
  expect_true("character" %in%
                class(mums2_example("230112_botryllus_peaktable.csv")))
  expect_true("character" %in% class(mums2_example("botryllus_v2.gnps.mgf")))
  expect_true("character" %in% class(mums2_example()))
  expect_true(length(mums2_example()) == 8)
  expect_error(mums2_example("dat"))
})
