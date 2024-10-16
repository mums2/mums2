/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include "Utils/Utils.h"


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Test cpp Utils") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("my_grep works") {
    // c("happy", "place", "Happy", "nothappy"), "happy"
    std::vector<std::string> vec{"happy", "place", "Happy", "nothappy"};
    Utils util;
    std::vector<int> result = util.my_grep(vec, "happy");
    
    std::vector<int> expected{0, 2};
    
    expect_true(result == expected);
  }

}
