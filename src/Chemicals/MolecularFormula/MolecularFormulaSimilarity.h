//
// Created by Gregory Johnson on 5/27/25.
//

#ifndef MOLECULARFORMULASIMILARITY_H
#define MOLECULARFORMULASIMILARITY_H
#include <Rcpp.h>
#include "MolecularMakeup.h"
#include "../../Math/VectorMath.h"

class MolecularFormulaSimilarity {
    public:
    static double ComputeSimilarity(const Rcpp::String& formula, const Rcpp::String& other) {
        if (formula.get_cstring()[0] == '\0' || other.get_cstring()[0] == '\0') return 0;
        
        const MolecularMakeup makeup(formula);
        const MolecularMakeup otherMolecularMakeup(other);
        const auto& alphabet = makeup.GetAlphabet();
        const auto& otherAlphabet =  otherMolecularMakeup.GetAlphabet();
        std::unordered_set<std::string> uniqueAlphabet;
        for (const auto& element : alphabet) {
            uniqueAlphabet.insert(element);
        }
        for (const auto& element : otherAlphabet) {
            uniqueAlphabet.insert(element);
        }
        const int size = static_cast<int>(uniqueAlphabet.size());
        Rcpp::NumericVector formulaAtoms(size);
        Rcpp::NumericVector otherFormulaAtoms(size);
        int count = 0;
        for (const auto& element : uniqueAlphabet) {
            formulaAtoms[count] = makeup.GetAtomsForElement(element);
            otherFormulaAtoms[count++] = otherMolecularMakeup.GetAtomsForElement(element);
        }
        return VectorMath::CosineScore(formulaAtoms, otherFormulaAtoms);
    }
};
#endif //MOLECULARFORMULASIMILARITY_H
