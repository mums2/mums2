//
// Created by gregj on 5/7/2025.
//


#include "Chemicals/MolecularFormula/MolecularFormula.h"

#include <map>
#include <unordered_map>
#include <unordered_set>

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
            if(!chemicalSymbol.empty()) {
                if (hasCheckedAtomAmount || !std::islower(molecularFormula[i])) {
                    if (amountOfAtoms.empty()) amountOfAtoms = "1";
                    chemicalAtomNamesOrder.emplace_back(chemicalSymbol);
                    chemicalAtomMap[chemicalSymbol] = std::stoi(amountOfAtoms);
                    hasCheckedAtomAmount = false;
                    chemicalSymbol = "";
                    amountOfAtoms = "";
                }
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
    chemicalAtomMap[chemicalSymbol] = std::stoi(amountOfAtoms);
    chemicalAtomNamesOrder.emplace_back(chemicalSymbol);
}

MolecularFormula::MolecularFormula(const std::unordered_map<std::string, int> &elementMap,
    const std::vector<std::string> &elementNamesOrder):chemicalAtomMap(elementMap),
chemicalAtomNamesOrder(elementNamesOrder) {}

int MolecularFormula::GetAtomsForElement(const std::string &chemicalElement) const {
    if(chemicalAtomMap.find(chemicalElement) == chemicalAtomMap.end()) return 0;
    return chemicalAtomMap.at(chemicalElement);
}

std::string MolecularFormula::GetMolecularFormula() const {
    std::string formula;
    for (const auto& element : chemicalAtomNamesOrder) {
        formula += element;
        const int amountOfAtoms = chemicalAtomMap.at(element);
        if (amountOfAtoms <= 1) continue;
        formula += std::to_string(amountOfAtoms) ;
    }
    return formula;
}

MolecularFormula MolecularFormula::operator-(const MolecularFormula &other) const {
    std::unordered_map<std::string, int> newElementMap;
    std::vector<std::string> tempElementNamesOrder;
    std::unordered_set<std::string> uniqueElementNames;
    std::vector<std::string> newElementNamesOrder;
    for(const auto& element : chemicalAtomNamesOrder) { // First one will have all uniques
        tempElementNamesOrder.emplace_back(element);
        uniqueElementNames.insert(element);

    }
    for (const auto& element : other.chemicalAtomNamesOrder) {
        if(uniqueElementNames.find(element) != uniqueElementNames.end()) continue;
        tempElementNamesOrder.emplace_back(element);
        uniqueElementNames.insert(element);
    }
    size_t count = 0;
    newElementNamesOrder.reserve(uniqueElementNames.size());
    for(const auto& currentChemicalSymbol : tempElementNamesOrder) {
        const int currentAtoms = GetAtomsForElement(currentChemicalSymbol);
        const int otherAtoms = other.GetAtomsForElement(currentChemicalSymbol);
        const int newAtomAmount = std::abs(currentAtoms - otherAtoms);
        if(newAtomAmount <= 0) continue;
        newElementMap[currentChemicalSymbol] = newAtomAmount;
        newElementNamesOrder.emplace_back(currentChemicalSymbol);
    }
    return MolecularFormula{newElementMap, newElementNamesOrder};
}

bool MolecularFormula::CheckIfSubformula(const MolecularFormula &subFormulaCandidate) const {
    // Cant have more elements than the main formula
    if (subFormulaCandidate.chemicalAtomNamesOrder.size() > chemicalAtomNamesOrder.size()) return false;
    std::unordered_set<std::string> uniqueElementNames;
    std::vector<std::string> newElementNamesOrder;
    for(const auto& element : chemicalAtomNamesOrder) { // First one will have all uniques
        uniqueElementNames.insert(element);
    }
    for (const auto& element : subFormulaCandidate.chemicalAtomNamesOrder) {
        if(uniqueElementNames.find(element) != uniqueElementNames.end()) continue;
        uniqueElementNames.insert(element);
    }
    return std::all_of(uniqueElementNames.begin(), uniqueElementNames.end(),
        [this, &subFormulaCandidate](const std::string& uniqueElement) {
            return GetAtomsForElement(uniqueElement) >= subFormulaCandidate.GetAtomsForElement(uniqueElement);
        });
}
