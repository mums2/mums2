//
// Created by Gregory Johnson on 5/30/25.
//
#include <testthat.h>
#include <Rcpp.h>
#include "FragmentationTree/GreedyHeuristic.h"

context("Greedy Heuristic") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.


    test_that("Greedy Heuristic returns the correct formula") {
        std::vector<double> score{0.3, 0.23, 0.1, 0.3};
        std::vector<double> mass{155.23, 155.22, 100.8, 120.12};
        std::vector<int> color{0,0,1,2};
        std::vector<std::string> formula{"C6H12O6", "C5H24O6", "C2H10", "C4H5O5"};
        DecompResult decompResult;
        decompResult.formula = {formula[0], formula[1]};
        decompResult.exactmass = {mass[0], mass[1]};
        decompResult.score = {score[0], score[1]};

        DecompResult decompResult1;
        decompResult1.formula = {formula[2]};
        decompResult1.exactmass = {mass[2]};
        decompResult1.score = {score[2]};

        DecompResult decompResult2;
        decompResult2.formula = {formula[3]};
        decompResult2.exactmass = {mass[3]};
        decompResult2.score = {score[3]};

        FragmentationTree tree({decompResult, decompResult1, decompResult2},
            155.24);

        tree.AddMolecularFormulaToGraph(0);
        tree.AddMolecularFormulaToGraph(1);
        tree.AddMolecularFormulaToGraph(2);
        tree.AddMolecularFormulaToGraph(3);
        expect_true(GreedyHeuristic::CalculateHeuristic(tree) == "C5H24O6");
    }

}