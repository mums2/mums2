//
// Created by Gregory Johnson on 11/11/25.
//

#include "HMDB/HumanMetabolomicsDB.h"
#include <ostream>
#include <Rcpp.h>
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
    for (auto &nodes : nodeMap) {
        const double precursorMz = nodes.second.precursorMz;
        for (const auto &spectraFile : spectraMap[nodes.first]) {
            Spectra spectra = ReadSpectra::ReadSpectraFile(spectraFile);
            spectra.precursorMz = precursorMz;
            nodes.second.spectraList.emplace_back(spectra);
        }
    }
}

void HumanMetabolomicsDB::PrintHumanMetabolomicsDB() {
    for (const auto &node : nodeMap) {
        Rcpp::Rcout << "Database " << node.first << std::endl;
        Rcpp::Rcout << "Spectras " << node.second.spectraList.size() << std::endl;
    }
}
