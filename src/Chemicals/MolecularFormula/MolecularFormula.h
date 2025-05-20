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
    explicit MolecularFormula(const Rcpp::String& molecularFormula);
    MolecularFormula(const std::unordered_map<std::string, int>&,
        const std::vector<std::string>&);
    double GetLossMass(const MolecularFormula& other) const;
    [[nodiscard]] int GetAtomsForElement(const std::string&) const;
    [[nodiscard]] std::string GetMolecularFormula() const;
    std::string operator-(const MolecularFormula& other) const;
    bool CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const;
    double GetMass() const;
protected:
    std::unordered_map<std::string, int> chemicalAtomMap;
    std::vector<std::string> chemicalAtomNamesOrder;
    std::unordered_map<std::string, double> chemicalAtomMassMap;
};



#endif //MOLECULARFORMULA_H
