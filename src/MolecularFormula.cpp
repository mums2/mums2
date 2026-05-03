//
// Created by gregj on 5/7/2025.
//


#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include <unordered_map>


std::vector<double> MolecularFormula::chemicalAtomMassVector = {12.011,
    1.0078, 14.007, 15.999, 30.974, 32.065}; // C H N O P S
std::vector<char> MolecularFormula::chemicalAtomNamesOrder = {'C', 'H', 'N', 'O', 'P', 'S'};

MolecularFormula::MolecularFormula(const Rcpp::String &molecularFormula, const double molecularMass):
molecularMass(molecularMass) {
    chemicalAtomsIndexTest = std::vector<int8_t>(100, -1);
    chemicalAtomsIndexTest[67] = 0;
    chemicalAtomsIndexTest[72] = 1;
    chemicalAtomsIndexTest[78] = 2;
    chemicalAtomsIndexTest[79] = 3;
    chemicalAtomsIndexTest[80] = 4;
    chemicalAtomsIndexTest[83] = 5;
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
                    // chemicalAtomsIndexTest[static_cast<int>(chemicalSymbol)]
                    chemicalAtomAmounts[chemicalAtomsIndexTest[static_cast<int>(chemicalSymbol)]] = std::stoi(amountOfAtoms);
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
    chemicalAtomAmounts[chemicalAtomsIndexTest[static_cast<int>(chemicalSymbol)]] = std::stoi(amountOfAtoms);
}

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
}


double MolecularFormula::GetLossMass(const MolecularFormula &other) const {
    return std::abs(GetMass() - other.GetMass());
}

int MolecularFormula::GetAtomsForElement(const char &chemicalElement) const {
    //return chemicalAtomMap.at(chemicalElement);
    const int8_t index = chemicalAtomsIndexTest[static_cast<int>(chemicalElement)];
    if (index < 0)
        Rcpp::stop("Chemical Element is Not CHNOPS");
    return chemicalAtomAmounts[index];
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
        const int8_t index = chemicalAtomsIndexTest[static_cast<int>(element)];
        const int atoms = std::abs(chemicalAtomAmounts[index] -
            other.chemicalAtomAmounts[index]);
        if (atoms <= 0) continue;
        formula += element;
        if (atoms == 1) continue;
        formula += std::to_string(atoms);
    }
    return formula;
}

bool MolecularFormula::CheckIfOtherIsSubFormula(const MolecularFormula &subFormulaCandidate) const {
    // O(M) function
    // With M being equal to the number of elements (in this case CHNOPS) M = 6
    return std::all_of(chemicalAtomNamesOrder.cbegin(), chemicalAtomNamesOrder.cend(),
        [&subFormulaCandidate, this](const char& element) {
            const int currentAtoms = GetAtomsForElement(element);
            const int otherAtoms = subFormulaCandidate.GetAtomsForElement(element);
            return currentAtoms >= otherAtoms;
    });

}

double MolecularFormula::GetMass() const {
    return molecularMass;
}
