//
// Created by gregj on 5/7/2025.
//

#ifndef MOLECULARFORMULA_H
#define MOLECULARFORMULA_H
#include <vector>
#include <string>
#include <unordered_map>
#include <Rcpp.h>

class MolecularFormula {
public:
    MolecularFormula() = default;
    explicit MolecularFormula(const Rcpp::String& molecularFormula, double molecularMass = 0);
    // MolecularFormula(const std::unordered_map<std::string, int>&,
    //     const std::vector<std::string>&);
    [[nodiscard]] double GetLossMass(const MolecularFormula& other) const;
    [[nodiscard]] int GetAtomsForElement(const char&) const;
    [[nodiscard]] std::string GetMolecularFormula() const;
    std::string operator-(const MolecularFormula& other) const;
    [[nodiscard]] int CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const;
    [[nodiscard]] double GetMass() const;
protected:
    double molecularMass;
    std::vector<int> chemicalAtomAmounts;
    static std::vector<char> chemicalAtomNamesOrder;
    static std::vector<double> chemicalAtomMassVector;
private:
    static size_t ConvertASCIIElementToIndex(int num);
};



#endif //MOLECULARFORMULA_H
