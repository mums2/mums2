//
// Created by gregj on 7/26/2025.
//
#include <testthat.h>
#include "DataStructures/CppMatrix.h"

context("Test CppMatrix Operations") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("You are able to add CppMatrix together") {
        std::vector<std::vector<double>> mat(3);
        mat[0] = {1,2,3};
        mat[1] = {4,5,6};
        mat[2] = {7,8,9};
        std::vector<std::vector<double>> mat2(3);
        mat2[0] = {1,2,3};
        mat2[1] = {4,5,6};
        mat2[2] = {7,8,9};
        const CppMatrix matrix(mat);
        const CppMatrix matrix2(mat2);
        const auto res = matrix + matrix2;
        std::vector<std::vector<double>> expectedResult(3);
        expectedResult[0] = {2,4,6};
        expectedResult[1] = {8,10,12};
        expectedResult[2] = {14,16,18};
        const CppMatrix expectedMatrix(expectedResult);
        expect_true(res == expectedMatrix);
    }
    test_that("CppMatrix Equivalence compares matrices properly") {
        std::vector<std::vector<double>> mat(3);
        mat[0] = {1,2,3};
        mat[1] = {4,5,6};
        mat[2] = {7,8,9};
        std::vector<std::vector<double>> mat2(3);
        mat2[0] = {1,2,3};
        mat2[1] = {4,5,6};
        mat2[2] = {7,8,9};
        std::vector<std::vector<double>> mat3(3);
        mat3[0] = {1,2,3};
        mat3[1] = {4,9,6};
        mat3[2] = {7,8,9};
        const CppMatrix matrix(mat);
        const CppMatrix matrix2(mat2);
        const CppMatrix matrix3(mat3);
        expect_true(matrix == matrix2);
        expect_false(matrix == matrix3);
        expect_false(matrix2 == matrix3);
    }
    test_that("You are able to subtract CppMatrix together") {
        std::vector<std::vector<double>> mat(3);
        mat[0] = {1,2,3};
        mat[1] = {4,5,6};
        mat[2] = {7,8,9};
        std::vector<std::vector<double>> mat2(3);
        mat2[0] = {1,2,3};
        mat2[1] = {4,5,6};
        mat2[2] = {7,8,9};
        const CppMatrix matrix(mat);
        const CppMatrix matrix2(mat2);
        const auto res = matrix - matrix2;
        std::vector<std::vector<double>> expectedResult(3);
        expectedResult[0] = {0,0,0};
        expectedResult[1] = {0,0,0};
        expectedResult[2] = {0,0,0};
        const CppMatrix expectedMatrix(expectedResult);
        expect_true
        (res == expectedMatrix);
    }

    test_that("Cpp Matrix errors when you perform operations on matrices of different sizes") {
        std::vector<std::vector<double>> mat(3);
        mat[0] = {1,2,3};
        mat[1] = {4,5,6};
        mat[2] = {7,8,9};
        std::vector<std::vector<double>> mat2(2);
        mat2[0] = {1,2,3};
        mat2[1] = {4,5,6};
        CppMatrix matrix(mat);
        const CppMatrix matrix2(mat2);
        expect_error(matrix - matrix2);
        expect_error(matrix + matrix2);
        expect_error(matrix += matrix2);
        expect_error(matrix == matrix2);
    }
    test_that("You are able to divide CppMatrix by a scaler") {
        std::vector<std::vector<double>> mat(3);
        mat[0] = {5,5,5};
        mat[1] = {10,10,10};
        mat[2] = {15,15,15};
        CppMatrix matrix(mat);
        std::vector<std::vector<double>> expectedResult(3);
        expectedResult[0] = {1,1,1};
        expectedResult[1] = {2,2,2};
        expectedResult[2] = {3,3,3};
        const CppMatrix expectedMatrix(expectedResult);
        matrix/=5;
        expect_true(matrix == expectedMatrix);
    }

}