//
// Created by gregj on 5/7/2025.
//


#include "Chemicals/MolecularFormula/MolecularFormula.h"

MolecularFormula::MolecularFormula(const std::string &molecularFormula) {
    const size_t size = molecularFormula.size();
    std::string chemicalSymbol;
    std::string amountOfAtoms;
    bool hasCheckedAtomAmount = false;
    for (size_t i = 0; i < size; i++) {
        if (!std::isdigit(molecularFormula[i])) {
            //If you have checked to see how many atoms there are or
            //The chemical does not have a lowercase letter, then the chemical before
            // Is standalone and only has 1 atom.
            // or it has checked the atoms and the chemical before has been created
            if (hasCheckedAtomAmount || !std::islower(molecularFormula[i])) {
                if (amountOfAtoms.empty()) amountOfAtoms = "1";
                elements.emplace_back(Element{chemicalSymbol, std::stoi(amountOfAtoms)});
                hasCheckedAtomAmount = false;
                chemicalSymbol = "";
                amountOfAtoms = "";
            }
            // If you are lowercase, you are apart of the current equation
            // if it is a chemical character
            chemicalSymbol += molecularFormula[i];
            continue;
        }
        amountOfAtoms += molecularFormula[i];
        hasCheckedAtomAmount = true;
    }
    if (amountOfAtoms.empty()) amountOfAtoms = "1";
    elements.emplace_back(Element{chemicalSymbol, std::stoi(amountOfAtoms)});
}

std::string MolecularFormula::GetMolecularFormula() const{
    std::string formula;
    for (const auto& element : elements) {
        formula += element.chemicalSymbol;
        if (element.atomCount <= 1) continue;
        formula += std::to_string(element.atomCount);
    }
    return formula;
}
