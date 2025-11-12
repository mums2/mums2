//
// Created by Gregory Johnson on 11/11/25.
//

#ifndef MUMS2_HUMANMETABOLOMICSDB_H
#define MUMS2_HUMANMETABOLOMICSDB_H
#include <list>
#include <unordered_map>

#include "HumanMetabolomicsDBNode.h"
#include "../Spectra/Spectra.h"
#include <vector>


class HumanMetabolomicsDB {
public:
    HumanMetabolomicsDB() = default;
    explicit HumanMetabolomicsDB(const std::vector<HumanMetabolomicsDBNode>& nodes);
    void AddHumanMetabolomicNode(const HumanMetabolomicsDBNode& node);
    void AddSpectraFiles(const std::string& spectraFiles, const std::string& databaseName);
    void ProcessSpectraFiles();
    Rcpp::List ConstructDataBase();
    void PrintHumanMetabolomicsDB();
private: // We could probably remove the unordered map and use a list. Add the spectra to the node
    std::unordered_map<std::string, HumanMetabolomicsDBNode> nodeMap;
    std::unordered_map<std::string, std::list<std::string>> spectraMap;
};


#endif //MUMS2_HUMANMETABOLOMICSDB_H
