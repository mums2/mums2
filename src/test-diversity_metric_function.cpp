//
// Created by Gregory Johnson on 1/29/25.
//
#include <testthat.h>
#include <Rcpp.h>
#include "DiversityMetrics/DiversityMetricFactory.h"

context("Diversity Metrics") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.


    test_that("Test Choose Diversity Metric Based on Name returns a non-null index") {
        //either shannon, simpson, or bray
        DiversityCalculator* shannon = DiversityMetricFactory::ChooseDiversityMetricBasedOnName("shannon");
        DiversityCalculator* simpson = DiversityMetricFactory::ChooseDiversityMetricBasedOnName("simpson");
        DiversityCalculator* bray = DiversityMetricFactory::ChooseDiversityMetricBasedOnName("bray");
        DiversityCalculator* isNull = DiversityMetricFactory::ChooseDiversityMetricBasedOnName("");

        expect_true(shannon != nullptr);
        expect_true(simpson != nullptr);
        expect_true(bray != nullptr);
        expect_true(isNull == nullptr);

        delete isNull;
        delete shannon;
        delete simpson;
        delete bray;
    }
    test_that("Test Choose Diversity Based on Index returns a non-null index") {
        // either beta or alpha diversity
        Diversity* shannon = DiversityMetricFactory::ChooseDiversityBasedOnIndex("shannon");
        Diversity* simpson = DiversityMetricFactory::ChooseDiversityBasedOnIndex("simpson");
        Diversity* bray = DiversityMetricFactory::ChooseDiversityBasedOnIndex("bray");
        Diversity* isNull = DiversityMetricFactory::ChooseDiversityBasedOnIndex("no");

        expect_true(shannon != nullptr);
        expect_true(simpson != nullptr);
        expect_true(bray != nullptr);
        expect_true(isNull == nullptr);

        delete isNull;
        delete shannon;
        delete simpson;
        delete bray;

    }

}