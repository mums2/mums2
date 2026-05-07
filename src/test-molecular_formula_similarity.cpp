//
// Created by Gregory Johnson on 7/11/25.
//

#include <testthat.h>
#include <Rcpp.h>
#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"

context("Molecular Formula Similarity") {

    test_that("Molecular Formula Similarity Computes Cosine Score as expected.") {
        expect_true(MolecularFormulaSimilarity::ComputeSimilarity("C6H12O6", "C6H12O6") == 1);
        expect_true(MolecularFormulaSimilarity::ComputeSimilarity("C6H12O6", "C3H6O3") == 1);
        expect_true(MolecularFormulaSimilarity::ComputeSimilarity("C6H12O6", "") == 0);
        expect_true(static_cast<int>(std::floor(MolecularFormulaSimilarity::ComputeSimilarity("C4H8",
            "C4H2") * 100)) == 80);
    }

}
