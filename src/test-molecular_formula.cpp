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
//        const MolecularFormula molecularFormula("C6H12O6");
//        bool result = molecularFormula.GetMolecularFormula() == "C6H12O6";
//        expect_true(result);
//        const MolecularFormula molecularFormula2("C20H15O7N3");
//        result = molecularFormula2.GetMolecularFormula() == "C20H15N3O7";
//        expect_true(result);
    }
    test_that("Get Atoms returns proper amount of atoms or 0 if not found") {
//        const MolecularFormula molecularFormula("C6H12O6");
//        expect_true(molecularFormula.GetAtomsForElement('C') == 6);
//        expect_true(molecularFormula.GetAtomsForElement('H') == 12);
//        expect_true(molecularFormula.GetAtomsForElement('O') == 6);
//        expect_error(molecularFormula.GetAtomsForElement('A') == 0);
    }
    test_that("We can subtract molecular formulas and get a resultant difference") {
//        const MolecularFormula molecularFormula("C6H12O6");
//        const MolecularFormula otherFormula("C20H15O7N3");
//        const std::string difference = molecularFormula - otherFormula;
//        expect_true(difference == "C14H3N3O");
//
//        const MolecularFormula molecularFormula2("C6H12O6");
//        const MolecularFormula otherFormula2("C3H2");
//        const std::string difference2 = molecularFormula2 - otherFormula2;
//        expect_true(difference2 == "C3H10O6");
    }

    test_that("We can check if formulas are subformulas of another") {
//        const MolecularFormula molecularFormula("C6H12O6");
//        const MolecularFormula otherFormula("C20H15O7N3");
//        expect_true(molecularFormula.CheckIfOtherIsSubFormula(otherFormula) == false);
//
//        const MolecularFormula molecularFormula2("C6H12O6");
//        const MolecularFormula otherFormula2("CH2");
//        expect_true(molecularFormula2.CheckIfOtherIsSubFormula(otherFormula2) == true);
//
//        const MolecularFormula molecularFormula3("C6H12O6");
//        const MolecularFormula otherFormula3("H10O5");
//        expect_true(molecularFormula3.CheckIfOtherIsSubFormula(otherFormula3) == true);
    }


}