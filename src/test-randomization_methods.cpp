//
// Created by Gregory Johnson on 1/29/25.
//
#include <testthat.h>
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"

context("Randomization Methods") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Test GetRandomVectorWithoutReplacement works as expected") {
        Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
        const Rcpp::Function setSeed = base["set.seed"];
        setSeed(25);
        // auto range = std::vector<int64_t>{1,4,6};
        // const auto res = RandomizationMethods::GetRandomVectorWithoutReplacement(range, 3, 11);
        // std::vector<size_t> result = {1,2,0};
        // expect_true(res == result);
    }

}