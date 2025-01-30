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
    test_that("Test GetRandomIndex works as expected") {
        Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
        const Rcpp::Function setSeed = base["set.seed"];
        setSeed(10);
        const auto res = RandomizationMethods::GetRandomNumberIndex({1,2,3}, 3, 6);
        expect_true(res == 2);
    }

}