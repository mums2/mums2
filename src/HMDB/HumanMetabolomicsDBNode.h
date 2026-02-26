//
// Created by Gregory Johnson on 11/11/25.
//

#ifndef MUMS2_HUMANMETABOLOMICSDBNODE_H
#define MUMS2_HUMANMETABOLOMICSDBNODE_H
#include <string>
#include <vector>
#include <Rcpp.h>
#include "../Spectra/Spectra.h"

struct HumanMetabolomicsDBNode {
     HumanMetabolomicsDBNode() = default;
     HumanMetabolomicsDBNode(const std::vector<std::string>& names,
          const std::vector<std::string>& dataValues) {
          int count = 0;
          int precursorMassIndex = -1;
          for (size_t i = 0; i < names.size(); i++) {
               if (count >= 4)
                    break;
               if (names[i] == "accession") {
                    databaseName = dataValues[i];
                    count++;
                    continue;
               }
               if (names[i] == "monisotopic_molecular_weight") {
                    precursorMassIndex = i;
                    if (dataValues[i] != "NA" && dataValues[i] != "NULL" && !dataValues[i].empty())
                         precursorMz = std::stod(dataValues[i]);
                    count++;
                    continue;
               }
               if (names[i] == "chemical_formula") {
                    chemicalFormula = dataValues[i];
                    count++;
                    continue;
               }
               if (names[i] == "name") {
                    annoName = dataValues[i];
                    count++;
               }
          }
          keys = names;
          keys[precursorMassIndex] = "precursormz";
          values = dataValues;
          if (values[precursorMassIndex].empty()) {
               values[precursorMassIndex] = "NA";
          }
     }
     std::vector<std::string> keys;
     std::vector<std::string> values;
     std::string databaseName;
     std::list<Spectra> spectraList;
     std::string annoName;
     std::string chemicalFormula;
     double precursorMz = -1;
};
#endif //MUMS2_HUMANMETABOLOMICSDBNODE_H