//
// Created by Gregory Johnson on 1/29/25.
//
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
#include "DiversityMetrics/AlphaDiversityCalculators/AlphaDiversity.h"
#include "Utils/Utils.h"


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Test Alpha Diversity") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Calculate Shannon Diversity returns proper results") {
        AlphaDiversity test;
        Rcpp::NumericVector v = {1,2,3,4};
        // Set the number of rows and columns to attribute dim of the vector object.
        v.attr("dim") = Rcpp::Dimension(2, 2);
        const Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(v);
        const Rcpp::NumericMatrix res = test.CalculateDiversity(mat, "shannon");
        bool resultBoolean = res(0,0) > 0 && res(0,1) > 0;
        expect_true(res.ncol() == 2);
        expect_true(res.nrow() == 1);
        expect_true(resultBoolean);
    }
    test_that("Calculate Simpson Diversity returns proper results") {
        AlphaDiversity test;
        Rcpp::NumericVector v = {1,1,1,1};
        // Set the number of rows and columns to attribute dim of the vector object.
        v.attr("dim") = Rcpp::Dimension(2, 2);
        const Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(v);
        const Rcpp::NumericMatrix res = test.CalculateDiversity(mat, "simpson");
        expect_true(res.ncol() == 2);
        expect_true(res.nrow() == 1);
        expect_true(res(0,0) == 1);
        expect_true(res(0,1) == 1);
    }

}
