//
// Created by Gregory Johnson on 7/11/25.
//

#include <testthat.h>
#include <Rcpp.h>
#include "Chemicals/MolecularFormula/MolecularMakeup.h"

context("Molecular Makeup") {

    test_that("Get Alaphabet returns the proper element.") {
        MolecularMakeup makeup("C6H12O6");
        const auto& alphabet = makeup.GetAlphabet();
        auto it = alphabet.begin();
        expect_true(*it++ == "C");
        expect_true(*it++ == "H");
        expect_true(*it == "O");
    }

    test_that("GetAtomsForElement returns the proper amount of atoms per element.") {
        MolecularMakeup makeup("C6H12O6");
        expect_true(makeup.GetAtomsForElement("C") == 6);
        expect_true(makeup.GetAtomsForElement("H") == 12);
        expect_true(makeup.GetAtomsForElement("O") == 6);
    }

}
