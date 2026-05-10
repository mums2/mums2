//
// Created by gregj on 5/7/2025.
//


#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include <unordered_map>


std::vector<char> MolecularFormula::chemicalAtomNamesOrder = {'C', 'H', 'N', 'O', 'P', 'S'};

MolecularFormula::MolecularFormula(const std::string& molecularFormula, const double molecularMass):
molecularMass(molecularMass) {
    chemicalAtomsIndexTest = std::vector<int8_t>(100, -1);
    chemicalAtomsIndexTest[67] = 0;
    chemicalAtomsIndexTest[72] = 1;
    chemicalAtomsIndexTest[78] = 2;
    chemicalAtomsIndexTest[79] = 3;
    chemicalAtomsIndexTest[80] = 4;
    chemicalAtomsIndexTest[83] = 5;
    chemicalAtomAmounts = std::vector<int>(6, 0);
    // const char* formula = molecularFormula.c_str();
    const size_t size = molecularFormula.size();
    char chemicalSymbol = ' ';
    std::string amountOfAtoms;
    bool hasCheckedAtomAmount = false;
    for (size_t i = 0; i < size; i++) {
        if (!std::isdigit(molecularFormula[i])) {
            //If you have checked to see how many atoms there are or
            //The chemical does not have a lowercase letter, then the chemical before
            // Is standalone and only has 1 atom.
            // or it has checked the atoms and the chemical before has been created
            if(chemicalSymbol != ' ') {
                if (hasCheckedAtomAmount || !std::islower(molecularFormula[i])) {
                    if (amountOfAtoms.empty()) amountOfAtoms = "1";
                    // chemicalAtomsIndexTest[static_cast<int>(chemicalSymbol)]
                    chemicalAtomAmounts[chemicalAtomsIndexTest[static_cast<int>(chemicalSymbol)]] = std::stoi(amountOfAtoms);
                    hasCheckedAtomAmount = false;
                    chemicalSymbol = ' ';
                    amountOfAtoms = "";
                }
            }
            // Only elements are CHNOPS
            // if it is a chemical character
            chemicalSymbol = molecularFormula[i];
            continue;
        }
        amountOfAtoms += molecularFormula[i];
        hasCheckedAtomAmount = true;
    }
    if (amountOfAtoms.empty()) amountOfAtoms = "1";
    chemicalAtomAmounts[chemicalAtomsIndexTest[static_cast<int>(chemicalSymbol)]] = std::stoi(amountOfAtoms);
    carbon = chemicalAtomAmounts[0];
    hydrogen = chemicalAtomAmounts[1];
    nitrogen = chemicalAtomAmounts[2];
    oxygen = chemicalAtomAmounts[3];
    phosphorus = chemicalAtomAmounts[4];
    sulfur = chemicalAtomAmounts[5];
}


double MolecularFormula::GetLossMass(const MolecularFormula &other) const {
    return std::abs(GetMass() - other.GetMass());
}

int MolecularFormula::GetAtomsForElement(const char &chemicalElement) const {
    //return chemicalAtomMap.at(chemicalElement);
    const int8_t index = chemicalAtomsIndexTest[static_cast<int>(chemicalElement)];
    if (index < 0)
        Rcpp::stop("Chemical Element is Not CHNOPS");
    return chemicalAtomAmounts[chemicalAtomsIndexTest[static_cast<int>(chemicalElement)]];
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


std::vector<int> MolecularFormula::operator-(const MolecularFormula &other) const {
    std::vector<int> results(chemicalAtomNamesOrder.size());
    // C H N O P S
    //{'C', 'H', 'N', 'O', 'P', 'S'};
    results[0] = std::abs(carbon - other.carbon);
    results[1] = std::abs(hydrogen - other.hydrogen);
    results[2] = std::abs(nitrogen - other.nitrogen);
    results[3] = std::abs(oxygen - other.oxygen);
    results[4] = std::abs(phosphorus - other.phosphorus);
    results[5] = std::abs(sulfur - other.sulfur);
    return results;
}

bool MolecularFormula::CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const {
    // To be a subformula the atoms of the given subformula candidate must be
    // less than or equal to the atoms of the parent formula candidate.
    if (carbon < subFormulaCandidate.carbon) return false;
    if (hydrogen < subFormulaCandidate.hydrogen) return false;
    if (nitrogen < subFormulaCandidate.nitrogen) return false;
    if (oxygen < subFormulaCandidate.oxygen) return false;
    if (phosphorus < subFormulaCandidate.phosphorus) return false;
    if (sulfur < subFormulaCandidate.sulfur) return false;
    return true;
}


double MolecularFormula::GetMass() const {
    return molecularMass;
}