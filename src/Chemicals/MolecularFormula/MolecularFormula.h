//
// Created by gregj on 5/7/2025.
//

#ifndef MOLECULARFORMULA_H
#define MOLECULARFORMULA_H
#include <vector>
#include <string>
#include "../PeriodicElements/Element.h"
#include <unordered_map>


class MolecularFormula {
public:
    explicit MolecularFormula(const std::string& molecularFormula);
    MolecularFormula(const std::unordered_map<std::string, int>&,
        const std::vector<std::string>&);
    int GetAtomsForElement(const std::string&) const;
    [[nodiscard]] std::string GetMolecularFormula() const;
    MolecularFormula operator-(const MolecularFormula& other) const;
private:
    std::unordered_map<std::string, int> chemicalAtomMap;
    std::vector<std::string> chemicalAtomNamesOrder;
};



#endif //MOLECULARFORMULA_H
