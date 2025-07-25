//
// Created by Gregory Johnson on 1/29/25.
//

#include <testthat.h>

#include "Rarefy/Rarefaction.h"
#include "Math/ParallelRandomNumberSitmo.h"

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
        std::vector<uint32_t> abundances = {3, 2, 1, 8};
        std::vector<uint32_t> abundanceRanges = {3, 5, 6, 14};
        const std::vector<uint32_t> eligible = {0, 1, 2, 3};
        const uint32_t size = 10;
        const uint32_t threshold = 2;
        const uint32_t currentSum = std::accumulate(abundances.begin(), abundances.end(), 0L);
        ParallelRandomNumberSitmo rngEngine(2120);
        const auto vec = rarefaction.Rarefy(abundances,
                eligible, abundanceRanges, rngEngine, size, currentSum, threshold);
        const auto sum = std::accumulate(vec.begin(), vec.end(), 0LL);
        std::vector<uint32_t> expected = {4,0,0,6};
        expect_true(sum == size);
        expect_true(vec.size() == 4);
        expect_true(vec == expected);
    }

}