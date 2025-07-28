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
    test_that("Calculate Simpson Diversity returns proper results") {
        std::vector<std::vector<double>> mat(1);
        mat[0] = {1,1,1};
        CppMatrix matrix(mat);
        AlphaDiversity test;
        std::vector<std::vector<double>> expectedResult(1);
        expectedResult[0] = {1};
        CppMatrix expectedMatrix(expectedResult);
        const CppMatrix res = test.CalculateDiversity(matrix, "simpson");
        bool val = expectedMatrix == res;
        expect_true(val == true);
    }
    test_that("Calculate Shannon Diversity returns proper results") {
        std::vector<std::vector<double>> mat(1);
        mat[0] = {1,1,1};
        CppMatrix matrix(mat);
        AlphaDiversity test;
        std::vector<std::vector<double>> expectedResult(1);
        expectedResult[0] = {1};
        CppMatrix expectedMatrix(expectedResult);
        const CppMatrix res = test.CalculateDiversity(matrix, "shannon");
        const Rcpp::NumericMatrix resultMatrix = res.ToRcppMatrix();
        bool val = resultMatrix(0,0) > 1;
        expect_true(val == true);
    }

}
