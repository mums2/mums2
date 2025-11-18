//
// Created by Gregory Johnson on 11/11/25.
//

#ifndef MUMS2_HUMANMETABOLOMICSDBNODE_H
#define MUMS2_HUMANMETABOLOMICSDBNODE_H
#include <string>
#include <vector>
#include <Rcpp.h>
#include "../Distance/DataStructures/Spectra.h"

struct HumanMetabolomicsDBNode {
     HumanMetabolomicsDBNode() = default;
     HumanMetabolomicsDBNode(const std::vector<std::string>& names,
          const std::vector<std::string>& dataValues) {
          int count = 0;
          int precursorMassIndex = 0;
          for (int i = 0; i < names.size(); i++) {
               if (count >= 2)
                    break;
               if (names[i] == "accession") {
                    databaseName = dataValues[i];
                    count++;
                    continue;
               }
               if (names[i] == "average_molecular_weight") {
                    precursorMassIndex = i;
                    if (values[i] == "NA" || values[i] == "NULL")
                         precursorMz = std::stod(values[i]);
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
     double precursorMz;
};
#endif //MUMS2_HUMANMETABOLOMICSDBNODE_H