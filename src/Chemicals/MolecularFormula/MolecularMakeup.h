//
// Created by Gregory Johnson on 5/27/25.
//

#ifndef MOLECULARMAKEUP_H
#define MOLECULARMAKEUP_H
#include <string>
#include <Rcpp.h>

class MolecularMakeup {
public:
    explicit MolecularMakeup(const std::string& molecularFormula);
    int GetAtomsForElement(const std::string &element) const;
    const std::list<std::string>& GetAlphabet() const {return alphabet;}
private:
    std::string fullFormula;
    std::list<std::string> alphabet;
    std::unordered_map<std::string, int> elementsAtomAmountMap;
};


#endif //MOLECULARMAKEUP_H
