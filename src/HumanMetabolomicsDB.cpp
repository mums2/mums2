//
// Created by Gregory Johnson on 11/11/25.
//

#include "HMDB/HumanMetabolomicsDB.h"
#include <ostream>
#include <Rcpp.h>

#include "CustomProgressBar/CliProgressBar.h"
#include "Spectra/ReadSpectra.h"

HumanMetabolomicsDB::HumanMetabolomicsDB(const size_t nodeSize) {
    nodes = std::vector<HumanMetabolomicsDBNode>(nodeSize);
    spectraFiles = std::vector<std::vector<std::string>>(nodeSize);
}

void HumanMetabolomicsDB::AddHumanMetabolomicNode(const HumanMetabolomicsDBNode &node, const size_t index) {
    nameToIndexMap[node.databaseName] = index;
    nodes[index] = node;
}

void HumanMetabolomicsDB::AddSpectraFiles(const std::string &spectraFile, const std::string& databaseName) {
    if (nameToIndexMap.find(databaseName) == nameToIndexMap.end()) return;
    spectraFiles[nameToIndexMap[databaseName]].emplace_back(spectraFile);
}

void HumanMetabolomicsDB::ProcessSpectraFiles() {
    CliProgressBar progressBar;
    const auto nodeMapCount = static_cast<float>(spectraFiles.size());
    float counter = 0;
    for (size_t i = 0; i < spectraFiles.size(); i++) {
        for (const auto &spectraFile : spectraFiles[i]) {
            Spectra spectra = ReadSpectra::ReadSpectraFile(spectraFile);
            nodes[i].spectraList.emplace_back(spectra);
        }
        counter++;
        progressBar.update(counter/nodeMapCount);
    }
}

AnnotationController* HumanMetabolomicsDB::ConstructDataBase() const {
    size_t referenceSize = 0;
    for (const auto& node : nodes) {
        if (node.spectraList.empty())
            continue;
        referenceSize += node.spectraList.size();
    }
    size_t count = 0;
    std::vector<AnnotationNode> databaseReferences(referenceSize);
    for (const auto& node : nodes) {
        if (node.spectraList.empty())
            continue;
        AnnotationNode annotation;
        annotation.name = node.annoName;
        annotation.precursorMz = node.precursorMz;
        annotation.chemicalFormula = node.chemicalFormula;
        for (size_t i = 0; i < node.keys.size(); i++) {
            annotation.keyValues.emplace_back(KeyValues{node.keys[i], node.values[i]});
            annotation.keys.emplace_back(node.keys[i]);
            annotation.values.emplace_back(node.values[i]);
        }
        for (const auto& spectra : node.spectraList) {
            annotation.spectra = spectra;
            databaseReferences[count++] = annotation;
        }
    }
    return new AnnotationController(databaseReferences);
}
