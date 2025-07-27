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
        std::vector<std::vector<double>> mat(3);
    	mat[0] = {1,2,3};
    	mat[1] = {1,2,3};
    	mat[2] = {1,2,3};
    	const CppMatrix matrix(mat);
    	std::vector<std::vector<double>> mat2(3);
    	mat2[0] = {0,0,0};
    	mat2[1] = {0,0,0};
    	mat2[2] = {0,0,0};
        const CppMatrix expectedResult(mat2);
    	BetaDiversity diversity;
    	const auto result = diversity.CalculateDiversity(matrix, "bray");
        expect_true(expectedResult == result);
    }

}