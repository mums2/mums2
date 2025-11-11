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
          for (int i = 0; i < names.size(); i++) {
               if (count >= 2)
                    break;
               if (names[i] == "accession") {
                    databaseName = dataValues[i];
                    count++;
                    continue;
               }
               if (names[i] == "average_molecular_weight") {
                    precursorMz = std::stod(dataValues[i]);
                    count++;
               }
          }
          keys = names;
          values = dataValues;
     }
     Rcpp::List ToList() {
          Rcpp::List list = Rcpp::List::create(
               Rcpp::Named("info") = Rcpp::List::create(
                    Rcpp::Named("key") = keys,
                    Rcpp::Named("value") = values),
               Rcpp::Named("spec") = Rcpp::List::create(
                    Rcpp::Named("mz") = spectraList.front().mz,
                    Rcpp::Named("intensity") = spectraList.front().intensity));
          return list;
     }
     std::vector<std::string> keys;
     std::vector<std::string> values;
     std::string databaseName;
     double precursorMz;
     std::list<Spectra> spectraList;
};
#endif //MUMS2_HUMANMETABOLOMICSDBNODE_H