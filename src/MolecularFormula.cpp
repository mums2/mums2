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
    chemicalAtomsIndexTest = std::vector<size_t>(100, -1);
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

double MolecularFormula::GetLossMass(const MolecularFormula &other) const {
    return std::abs(GetMass() - other.GetMass());
}

int MolecularFormula::GetAtomsForElement(const char &chemicalElement) const {
    //return chemicalAtomMap.at(chemicalElement);
    const size_t index = chemicalAtomsIndexTest[static_cast<int>(chemicalElement)];
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
        const size_t index = chemicalAtomsIndexTest[static_cast<int>(element)];
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
    // for (const auto& element : chemicalAtomNamesOrder) {
    //     const int currentAtoms = GetAtomsForElement(element);
    //     const int otherAtoms = subFormulaCandidate.GetAtomsForElement(element);
    //     if (currentAtoms < otherAtoms) return false;
    // }
    return std::all_of(chemicalAtomNamesOrder.cbegin(), chemicalAtomNamesOrder.cend(),
        [&subFormulaCandidate, this](const char& element) {
            const int currentAtoms = GetAtomsForElement(element);
            const int otherAtoms = subFormulaCandidate.GetAtomsForElement(element);
            return currentAtoms >= otherAtoms;
    });
    // // so if thisformula and otherformula are both true, it just returns one, since other formula and this formula
    // // are subformulas of each other
    // // But we only represent the first formula with the link so its not cyclic
    // if (thisFormula) return 1;
    // // if it reached this far, the only option is that the other formula is true so we return a 2
    // return 2;

}

double MolecularFormula::GetMass() const {
    return molecularMass;
}

size_t MolecularFormula::ConvertASCIIElementToIndex(const int num) {
    if (num == 67) return 0; // C
    if (num == 72) return 1; // H
    if (num == 78) return 2; // N
    if (num == 79) return 3; // O
    if (num == 80) return 4; // P
    if (num == 83) return 5; // S
    //if 83 which it all it can be
    // We screen before-hand and the alphabet only contains CHNOPS
    Rcpp::stop("Chemical Element is Not CHNOPS");
}
