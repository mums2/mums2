//
// Created by Gregory Johnson on 5/30/25.
//
#include <testthat.h>
#include <Rcpp.h>
#include "FragmentationTree/FragmentationTree.h"

context("Fragmentation Tree") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.


    test_that("Fragmentation Tree initialize correctly") {
        std::vector<double> score{0.3, 0.23, 0.1, 0.3};
        std::vector<double> mass{155.23, 155.22, 100.8, 120.12};
        std::vector<int> color{0,0,1,2};
        std::vector<std::string> formula{"C6H12O6", "C5H24O6", "C2H10", "C4H5O5"};
        Rcpp::List data = Rcpp::List::create(Rcpp::Named("formula") = formula,
            Rcpp::Named("mass") = mass,
            Rcpp::Named("color") = color,
            Rcpp::Named("score") = score);
        FragmentationTree tree(data, 155.24);
        const auto nodes = tree.GetFragmentationNodes();
        expect_true(nodes.size() == 4);
        expect_true(tree.GetFragmentationNodes()[0].score == score[0]);
        expect_true(tree.GetFragmentationNodes()[1].score == score[1]);
        expect_true(tree.GetFragmentationNodes()[2].score == score[2]);
        expect_true(tree.GetFragmentationNodes()[3].score == score[3]);
    }
    test_that("Fragment Tree algorithm is computing properly") {
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
        const auto nodes = tree.GetFragmentationNodes();
        expect_true(nodes[0].subTreeScore != 0);
        expect_true(nodes[1].subTreeScore != 0);
        expect_true(nodes[2].subTreeScore == 0);
        expect_true(nodes[3].subTreeScore == 0);
    }
    test_that("Fragmentation tree is sorting the nodes correctly") {
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
        const std::string result = tree.GetBestFormula();
        expect_true(result == "C6H12O6");
    }


}