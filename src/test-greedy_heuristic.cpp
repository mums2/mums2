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
        Rcpp::List data = Rcpp::List::create(Rcpp::Named("formula") = formula,
            Rcpp::Named("mass") = mass,
            Rcpp::Named("color") = color,
            Rcpp::Named("score") = score);
        FragmentationTree tree(data, 155.24);

        tree.AddMolecularFormulaToGraph(0);
        tree.AddMolecularFormulaToGraph(1);
        tree.AddMolecularFormulaToGraph(2);
        tree.AddMolecularFormulaToGraph(3);
        Rcpp::Rcout << tree.GetFragmentationNodes()[0].subTreeScore << std::endl;
        Rcpp::Rcout << tree.GetFragmentationNodes()[1].subTreeScore << std::endl;
        expect_true(GreedyHeuristic::CalculateHeuristic(tree) == "C5H24O6");
    }

}