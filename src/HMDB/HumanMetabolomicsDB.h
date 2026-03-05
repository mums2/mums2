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

#include "../AnnotationStructure/AnnotationController.h"


class HumanMetabolomicsDB {
public:
    explicit HumanMetabolomicsDB(size_t nodeSize);
    void AddHumanMetabolomicNode(const HumanMetabolomicsDBNode& node, size_t index);
    void AddSpectraFiles(const std::string &spectraFile, const std::string& databaseName);
    void ProcessSpectraFiles();
    AnnotationController* ConstructDataBase() const;
private: // We could probably remove the unordered map and use a list. Add the spectra to the node
    std::vector<HumanMetabolomicsDBNode> nodes;
    std::vector<std::vector<std::string>> spectraFiles;
    std::unordered_map<std::string, size_t> nameToIndexMap;
};


#endif //MUMS2_HUMANMETABOLOMICSDB_H
