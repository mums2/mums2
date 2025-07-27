//
// Created by Gregory Johnson on 1/29/25.
//
#include <testthat.h>
#include "DiversityMetrics/BetaDiversityCalculators/BetaDiversity.h"
#include "Utils/Utils.h"


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Test Beta Diversity") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Calculate Bray Curtis returns proper results") {
//        BetaDiversity test;
//        Rcpp::NumericVector v = {1,2,3,4};
//        // Set the number of rows and columns to attribute dim of the vector object.
//        v.attr("dim") = Rcpp::Dimension(2, 2);
//        const Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(v);
//        const Rcpp::NumericMatrix res = test.CalculateDiversity(mat, "bray");
//        bool isGreaterThanZero = (res(0,1) > 0 && res(1,0) > 0);
//        bool isEqualThanZero = (res(0,0) == 0 && res(1,1) == 0);
//        expect_true(res.ncol() == 2);
//        expect_true(res.nrow() == 2);
//        expect_true(res(0,0) == res(1,1));
//        expect_true(res(0,1) == res(1,0));
//        expect_true(isGreaterThanZero);
//        expect_true(isEqualThanZero);
    }

}