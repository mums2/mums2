//
// Created by gregj on 5/7/2025.
//

#ifndef MOLECULARFORMULA_H
#define MOLECULARFORMULA_H
#include <vector>
#include <string>
#include "../PeriodicElements/Element.h"


class MolecularFormula {
public:
    explicit MolecularFormula(const std::string& molecularFormula);
    const std::vector<Element>& GetElements() {return elements;}
    std::string GetMolecularFormula() const;
private:
    std::vector<Element> elements;
};



#endif //MOLECULARFORMULA_H
