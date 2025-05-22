//
// Created by gregj on 5/7/2025.
//


#include "Chemicals/MolecularFormula/MolecularFormula.h"

#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>


std::vector<double> MolecularFormula::chemicalAtomMassVector = {12.011,
    1.0078, 14.007, 15.999, 30.974, 32.065}; // C H N O P S
std::vector<char> MolecularFormula::chemicalAtomNamesOrder = {'C', 'H', 'N', 'O', 'P', 'S'};

MolecularFormula::MolecularFormula(const Rcpp::String &molecularFormula) {
    chemicalAtomAmounts = std::vector<int>(6, 0);
    const char* formula = molecularFormula.get_cstring();
    const size_t size = std::strlen(formula);
    char chemicalSymbol = ' ';
    std::string amountOfAtoms;
    bool hasCheckedAtomAmount = false;
    for (size_t i = 0; i < size; i++) {
        if (!std::isdigit(formula[i])) {
            //If you have checked to see how many atoms there are or
            //The chemical does not have a lowercase letter, then the chemical before
            // Is standalone and only has 1 atom.
            // or it has checked the atoms and the chemical before has been created
            if(chemicalSymbol != ' ') {
                if (hasCheckedAtomAmount || !std::islower(formula[i])) {
                    if (amountOfAtoms.empty()) amountOfAtoms = "1";
                    chemicalAtomAmounts[ConvertASCIIElementToIndex(chemicalSymbol)] = std::stoi(amountOfAtoms);
                    hasCheckedAtomAmount = false;
                    chemicalSymbol = ' ';
                    amountOfAtoms = "";
                }
            }
            // Only elements are CHNOPS
            // if it is a chemical character
            chemicalSymbol = formula[i];
            continue;
        }
        amountOfAtoms += formula[i];
        hasCheckedAtomAmount = true;
    }
    if (amountOfAtoms.empty()) amountOfAtoms = "1";
    chemicalAtomAmounts[ConvertASCIIElementToIndex(chemicalSymbol)] = std::stoi(amountOfAtoms);
}

// MolecularFormula::MolecularFormula(const std::unordered_map<std::string, int> &elementMap,
//     const std::vector<std::string> &elementNamesOrder):chemicalAtomMap(elementMap),
// chemicalAtomNamesOrder(elementNamesOrder) {}

double MolecularFormula::GetLossMass(const MolecularFormula &other) const {
    return std::abs(GetMass() - other.GetMass());
}

int MolecularFormula::GetAtomsForElement(const char &chemicalElement) const {
    //return chemicalAtomMap.at(chemicalElement);
    return chemicalAtomAmounts[ConvertASCIIElementToIndex(chemicalElement)];
}

std::string MolecularFormula::GetMolecularFormula() const {
    std::string formula;
    for (const auto& element : chemicalAtomNamesOrder) {

        const int amountOfAtoms = GetAtomsForElement(element);
        if (amountOfAtoms <= 0) continue;
        formula += element;
        if (amountOfAtoms == 1) continue;
        formula += std::to_string(amountOfAtoms) ;
    }
    return formula;
}

std::string MolecularFormula::operator-(const MolecularFormula &other) const {
    std::string formula;
    for (const auto& element : chemicalAtomNamesOrder) {
        const int atoms = std::abs(chemicalAtomAmounts[ConvertASCIIElementToIndex(element)] -
            other.chemicalAtomAmounts[ConvertASCIIElementToIndex(element)]);
        if (atoms <= 0) continue;
        formula += element;
        if (atoms == 1) continue;
        formula += std::to_string(atoms);
    }
    return formula;
}

int MolecularFormula::CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const {
    if (subFormulaCandidate.chemicalAtomNamesOrder.size() > chemicalAtomNamesOrder.size()) return false;
    bool thisFormula = true;
    bool otherFormula = true;
    for (const auto& element : chemicalAtomNamesOrder) {
        if (!thisFormula && !otherFormula) return 0; // Neither is a subformula
        if (GetAtomsForElement(element) == subFormulaCandidate.GetAtomsForElement(element)) continue;
        if (GetAtomsForElement(element) > subFormulaCandidate.GetAtomsForElement(element)) otherFormula = false;
        else thisFormula = false;
    }
    // so if thisformula and otherformula are both true, it just returns one, since other formula and this formula
    // are subformulas of eachother
    // But we only represent the first formula with the link so its not cyclic
    if (thisFormula) return 1;
    // if it reached this far, the only option is that the other formula is true so we return a 2
    return 2;

}

double MolecularFormula::GetMass() const {
    double mass = 0;
    for (size_t i = 0; i < chemicalAtomAmounts.size(); i++) {
        mass += chemicalAtomAmounts[i] * chemicalAtomMassVector[i];
    }
    return mass;
}

size_t MolecularFormula::ConvertASCIIElementToIndex(int num) const {
    if (num == 67) return 0; // C
    if (num == 72) return 1; // H
    if (num == 78) return 2; // N
    if (num == 79) return 3; // O
    if (num == 80) return 4; // P
    if (num == 83) return 5; // S
    //if 83 which it all it can be
    // We screen before-hand and the alphabet only contains CHNOPS
    Rcpp::stop("Chemical Element is Not CHNOPS");
    return -1; // error
}
