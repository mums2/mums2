// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
#include <testthat.h>
#include "Spectra/ReadSpectra.h"
context("Test Read Spectra") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Read Spectra will fail if given a file that does not exist") {
        expect_error(ReadSpectra::ReadSpectraFile(""));
    }

}