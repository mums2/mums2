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
    explicit MolecularFormula(const std::string& molecularFormula, double molecularMass = 0);
    double GetLossMass(const MolecularFormula& other) const;
    int GetAtomsForElement(const char&) const;
    std::string GetMolecularFormula() const;
    std::vector<int> operator-(const MolecularFormula& other) const;
    bool CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const;
    double GetMass() const;
protected:
    double molecularMass{};
    std::vector<int> chemicalAtomAmounts;
    static std::vector<char> chemicalAtomNamesOrder;
    static std::vector<double> chemicalAtomMassVector;
    std::vector<int8_t> chemicalAtomsIndexTest;
    int carbon;
    int hydrogen;
    int nitrogen;
    int oxygen;
    int phosphorus;
    int sulfur;
};



#endif //MOLECULARFORMULA_H
