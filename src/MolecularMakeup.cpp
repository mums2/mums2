//
// Created by Gregory Johnson on 5/27/25.
//

#include "Chemicals/MolecularFormula/MolecularMakeup.h"

MolecularMakeup::MolecularMakeup(const std::string& molecularFormula) {

    std::string formula = molecularFormula;
    formula.erase(std::remove_if(formula.begin(), formula.end(),
        ::isspace), formula.end());
    const size_t size = formula.size();
    std::string chemicalSymbol;
    std::string amountOfAtoms;
    bool hasCheckedAtomAmount = false;
    for (size_t i = 0; i < size; i++) {
        if (!std::isdigit(formula[i])) {
            //If you have checked to see how many atoms there are or
            //The chemical does not have a lowercase letter, then the chemical before
            // Is standalone and only has 1 atom.
            // or it has checked the atoms and the chemical before has been created
            if(!chemicalSymbol.empty()) {
                if (hasCheckedAtomAmount || !std::islower(formula[i])) {
                    if (amountOfAtoms.empty()) amountOfAtoms = "1";
                    fullFormula += amountOfAtoms;
                    alphabet.emplace_back(chemicalSymbol);
                    elementsAtomAmountMap[chemicalSymbol] = std::stoi(amountOfAtoms);
                    hasCheckedAtomAmount = false;
                    chemicalSymbol = "";
                    amountOfAtoms = "";
                }
            }
            // Only elements are CHNOPS
            // if it is a chemical character
            chemicalSymbol += formula[i];
            fullFormula += formula[i];
            continue;
        }
        amountOfAtoms += formula[i];
        hasCheckedAtomAmount = true;
    }
    if (amountOfAtoms.empty()) amountOfAtoms = "1";
    elementsAtomAmountMap[chemicalSymbol] = std::stoi(amountOfAtoms);
    fullFormula += amountOfAtoms;
    alphabet.emplace_back(chemicalSymbol);
}

int MolecularMakeup::GetAtomsForElement(const std::string &element) const {
    if (elementsAtomAmountMap.find(element) == elementsAtomAmountMap.end()) return 0;
    return elementsAtomAmountMap.at(element);
}
