//
// Created by gregj on 5/7/2025.
//


#include "Chemicals/MolecularFormula/MolecularFormula.h"

#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>

MolecularFormula::MolecularFormula(const Rcpp::String &molecularFormula) {
    const char* formula = molecularFormula.get_cstring();
    const size_t size = std::strlen(formula);
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
                    chemicalAtomNamesOrder.emplace_back(chemicalSymbol);
                    chemicalAtomMap[chemicalSymbol] = std::stoi(amountOfAtoms);
                    hasCheckedAtomAmount = false;
                    chemicalSymbol = "";
                    amountOfAtoms = "";
                }
            }
            // If you are lowercase, you are apart of the current equation
            // if it is a chemical character
            chemicalSymbol += formula[i];
            continue;
        }
        amountOfAtoms += formula[i];
        hasCheckedAtomAmount = true;
    }
    if (amountOfAtoms.empty()) amountOfAtoms = "1";
    chemicalAtomMap[chemicalSymbol] = std::stoi(amountOfAtoms);
    chemicalAtomNamesOrder.emplace_back(chemicalSymbol);
    delete formula;
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

std::string MolecularFormula::operator-(const MolecularFormula &other) const {
    std::vector<std::string> tempElementNamesOrder;
    std::unordered_set<std::string> uniqueElementNames;

    for(const auto& element : chemicalAtomNamesOrder) { // First one will have all uniques
        tempElementNamesOrder.emplace_back(element);
        uniqueElementNames.insert(element);
    }
    for (const auto& element : other.chemicalAtomNamesOrder) {
        if(uniqueElementNames.find(element) != uniqueElementNames.end()) continue;
        tempElementNamesOrder.emplace_back(element);
        uniqueElementNames.insert(element);
    }
    std::string formula;
    for(const auto& currentChemicalSymbol : tempElementNamesOrder) {
        const int currentAtoms = GetAtomsForElement(currentChemicalSymbol);
        const int otherAtoms = other.GetAtomsForElement(currentChemicalSymbol);
        const int newAtomAmount = std::abs(currentAtoms - otherAtoms);
        if(newAtomAmount <= 0) continue;
        if(newAtomAmount == 1){
            formula += currentChemicalSymbol;
            continue;
        }
        formula += (currentChemicalSymbol + std::to_string(newAtomAmount));
    }
    return formula;
}

bool MolecularFormula::CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const {
    // Cant have more elements than the main formula
    if (subFormulaCandidate.chemicalAtomNamesOrder.size() > chemicalAtomNamesOrder.size()) return false;
    // for (const auto& element : chemicalAtomMap) {
    //     int atoms1 = GetAtomsForElement(element.first);
    //     int atoms2 = subFormulaCandidate.GetAtomsForElement(element.first);
    //     if (GetAtomsForElement(element.first) < subFormulaCandidate.GetAtomsForElement(element.first)) return false;
    // }
    return std::all_of(chemicalAtomMap.begin(), chemicalAtomMap.end(),
        [this, &subFormulaCandidate](const std::pair<const std::string, int>& element) {
            return GetAtomsForElement(element.first) >= subFormulaCandidate.GetAtomsForElement(element.first);
        });
}
