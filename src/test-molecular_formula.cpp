//
// Created by Gregory Johnson on 1/29/25.
//
#include <testthat.h>
#include <Rcpp.h>
#include "Chemicals/MolecularFormula/MolecularFormula.h"

context("Molecular Formula") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.

    test_that("Test that MolecularFormula creates the table for a given string and outputs it properly") {
        const std::string formula = "C6H12O6";
        const MolecularFormula molecularFormula(formula);
        bool result = molecularFormula.GetMolecularFormula() == "C6H12O6";
        expect_true(result);
        const std::string formula1 = "C20H15O7N3";
        const MolecularFormula molecularFormula2(formula1);
        result = molecularFormula2.GetMolecularFormula() == "C20H15N3O7";
        expect_true(result);
    }
    test_that("Get Atoms returns proper amount of atoms or 0 if not found") {
        const std::string formula = "C6H12O6";
        const MolecularFormula molecularFormula(formula);
        expect_true(molecularFormula.GetAtomsForElement('C') == 6);
        expect_true(molecularFormula.GetAtomsForElement('H') == 12);
        expect_true(molecularFormula.GetAtomsForElement('O') == 6);
        expect_error(molecularFormula.GetAtomsForElement('A') == 0);
    }
    test_that("We can subtract molecular formulas and get a resultant difference") {
        const std::string formula = "C6H12O6";
        const std::string formula1 = "C20H15O7N3";
        const MolecularFormula molecularFormula(formula);
        const MolecularFormula otherFormula(formula1);
        const std::vector<int> difference = molecularFormula - otherFormula;
        // C14H3N3O
        const std::vector<int> expected{14, 3, 3, 1, 0, 0};
        expect_true(difference == expected);

        const std::string formula2 = "C3H2";
        const MolecularFormula molecularFormula2(formula);
        const MolecularFormula otherFormula2(formula2);
        const std::vector<int> difference2 = molecularFormula2 - otherFormula2;
        // "C3H10O6"
        const std::vector<int> expected1{3, 10, 0, 6, 0, 0};
        expect_true(difference2 == expected1);
    }

    test_that("We can check if formulas are subformulas of another") {

        const std::string formula = "C6H12O6";
        const std::string formula1 = "C20H15O7N3";
        const std::string formula2 = "CH2";
        const std::string formula3 = "H10O5";
        const MolecularFormula molecularFormula(formula);
        const MolecularFormula otherFormula(formula1);
        expect_true(molecularFormula.CheckIfOtherIsSubFormula(otherFormula) == false);

        const MolecularFormula molecularFormula2(formula);
        const MolecularFormula otherFormula2(formula2);
        expect_true(molecularFormula2.CheckIfOtherIsSubFormula(otherFormula2) == true);

        const MolecularFormula molecularFormula3(formula);
        const MolecularFormula otherFormula3(formula3);
        expect_true(molecularFormula3.CheckIfOtherIsSubFormula(otherFormula3) == true);
    }


}