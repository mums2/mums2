//
// Created by Gregory Johnson on 11/11/25.
//

#include "HMDB/HumanMetabolomicsDB.h"
#include <ostream>
#include <Rcpp.h>

#include "CustomProgressBar/CliProgressBar.h"
#include "Spectra/ReadSpectra.h"

HumanMetabolomicsDB::HumanMetabolomicsDB(const std::vector<HumanMetabolomicsDBNode> &nodes) {
    for (const auto &node : nodes) {
        nodeMap[node.databaseName] = node;
    }
}

void HumanMetabolomicsDB::AddHumanMetabolomicNode(const HumanMetabolomicsDBNode &node) {
    nodeMap[node.databaseName] = node;
}

void HumanMetabolomicsDB::AddSpectraFiles(const std::string &spectraFiles, const std::string& databaseName) {
    spectraMap[databaseName].emplace_back(spectraFiles);
}

void HumanMetabolomicsDB::ProcessSpectraFiles() {
    CliProgressBar progressBar;
    Rcpp::Rcout << "Processing spectra files..." << std::endl;
    const auto nodeMapCount = static_cast<float>(nodeMap.size());
    float counter = 0;
    for (auto &nodes : nodeMap) {
        for (const auto &spectraFile : spectraMap[nodes.first]) {
            Spectra spectra = ReadSpectra::ReadSpectraFile(spectraFile);
            nodes.second.spectraList.emplace_back(spectra);
        }
        counter++;
        progressBar.update(counter/nodeMapCount);
    }
}

void HumanMetabolomicsDB::PrintHumanMetabolomicsDB() {
    for (const auto &node : nodeMap) {
        Rcpp::Rcout << "Database " << node.first << std::endl;
        Rcpp::Rcout << "Spectras " << node.second.spectraList.size() << std::endl;
    }
}

AnnotationController* HumanMetabolomicsDB::ConstructDataBase() const {
    size_t referenceSize = 0;
    for (const auto& node : nodeMap) {
        if (node.second.spectraList.empty())
            continue;
        referenceSize += node.second.spectraList.size();
    }
    size_t count = 0;
    std::vector<AnnotationNode> databaseReferences(referenceSize);
    for (const auto& node : nodeMap) {
        if (node.second.spectraList.empty())
            continue;
        AnnotationNode annotation;
        annotation.name = node.second.annoName;
        annotation.precursorMz = node.second.precursorMz;
        annotation.chemicalFormula = node.second.chemicalFormula;
        for (size_t i = 0; i < node.second.keys.size(); i++) {
            annotation.keyValues.emplace_back(KeyValues{node.second.keys[i], node.second.values[i]});
            annotation.keys.emplace_back(node.second.keys[i]);
            annotation.values.emplace_back(node.second.values[i]);
        }
        for (const auto& spectra : node.second.spectraList) {
            annotation.spectra = spectra;
            databaseReferences[count++] = annotation;
        }
    }
    return new AnnotationController(databaseReferences);
}
