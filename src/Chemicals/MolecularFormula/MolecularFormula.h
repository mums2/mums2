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
    double GetLossMass(const MolecularFormula& other) const;
    int GetAtomsForElement(const char&) const;
    std::string GetMolecularFormula() const;
    std::string operator-(const MolecularFormula& other) const;
    int CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const;
    double GetMass() const;
protected:
    double molecularMass;
    std::vector<int> chemicalAtomAmounts;
    static std::vector<char> chemicalAtomNamesOrder;
    static std::vector<double> chemicalAtomMassVector;
private:
    static size_t ConvertASCIIElementToIndex(int num);
};



#endif //MOLECULARFORMULA_H
