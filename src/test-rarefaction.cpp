//
// Created by Gregory Johnson on 1/29/25.
//

#include <testthat.h>

#include "Rarefy/Rarefaction.h"


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Test Rarefaction") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Rarefaction will rarefy the numeric matrix properly") {
        Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
        const Rcpp::Function setSeed = base["set.seed"];
        setSeed(2);
        Rarefaction rarefaction;
        std::vector<int64_t> v = {1, 2, 3, 4, 5, 6};
        std::vector<int64_t> eligbile = {0, 1, 2, 3, 4, 5};
        const auto vec = rarefaction.Rarefy({}, v,
                eligbile, v, 10, 2);
        const auto sum = std::accumulate(vec.begin(), vec.end(), 0LL);
        std::vector<int64_t> expected = {0, 0, 3, 0, 3, 4};
        expect_true(sum == 10);
        expect_true(vec.size() == 6);
        expect_true(vec == expected);
    }

}