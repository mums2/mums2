//
// Created by Gregory Johnson on 5/27/25.
//

#ifndef MOLECULARFORMULASIMILARITY_H
#define MOLECULARFORMULASIMILARITY_H
#include <Rcpp.h>
#include "MolecularMakeup.h"
#include "../../Math/VectorMath.h"
#include <string>

class MolecularFormulaSimilarity {
    public:
    static double ComputeSimilarity(const std::string& formula, const std::string& other) {
        if (formula.empty() || other.empty()) return 0;
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
        std::vector<double> formulaAtoms(size);
        std::vector<double> otherFormulaAtoms(size);
        int count = 0;
        for (const auto& element : uniqueAlphabet) {
            formulaAtoms[count] = makeup.GetAtomsForElement(element);
            otherFormulaAtoms[count++] = otherMolecularMakeup.GetAtomsForElement(element);
        }
        return VectorMath::CosineScore(formulaAtoms, otherFormulaAtoms);
    }
};
#endif //MOLECULARFORMULASIMILARITY_H
