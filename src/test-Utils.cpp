// /*
//  * This file uses the Catch unit testing library, alongside
//  * testthat's simple bindings, to test a C++ function.
//  *
//  * For your own packages, ensure that your test files are
//  * placed within the `src/` folder, and that you include
//  * `LinkingTo: testthat` within your DESCRIPTION file.
//  */
//
// // All test files should include the <testthat.h>
// // header file.
// #include <testthat.h>
// #include "Utils/Utils.h"
//
//
// // Initialize a unit test context. This is similar to how you
// // might begin an R test file with 'context()', expect the
// // associated context should be wrapped in braced.
// context("Test cpp Utils") {
//
//   // The format for specifying tests is similar to that of
//   // testthat's R functions. Use 'test_that()' to define a
//   // unit test, and use 'expect_true()' and 'expect_false()'
//   // to test the desired conditions.
//   test_that("my_grep works") {
//     // c("happy", "place", "Happy", "nothappy"), "happy"
//     std::vector<std::string> vec{"happy", "place", "Happy", "nothappy"};
//     Utils util;
//     std::vector<int> result = util.my_grep(vec, "happy");
//
//     std::vector<int> expected{0, 2};
//
//     expect_true(result == expected);
//   }
//
//     test_that("testing my_rep works") {
//     std::vector<int> values(3);
//     values[0] = 1;
//     values[1] = 2;
//     values[2] = 3;
//
//     std::vector<int> times(3);
//     times[0] = 2;
//     times[1] = 3;
//     times[2] = 2;
//
//     // out: 1, 1, 2, 2, 2, 3, 3
//     std::vector<int> out{1, 1, 2, 2, 2, 3, 3};
//     Utils util;
//
//     std::vector<int> rep = util.my_rep(values, times);
//     Rcpp::Rcout << rep[0] << std::endl;
//     expect_true(rep == out);
//   }
//
// }
