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
        for (const auto &spectraFile : spectraMap[nodes.first]) {
            Spectra spectra = ReadSpectra::ReadSpectraFile(spectraFile);
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

Rcpp::List HumanMetabolomicsDB::ConstructDataBase() {

    size_t size = 0;
    for (const auto& node : nodeMap) {
        if (node.second.spectraList.empty())
            continue;
        size += node.second.spectraList.size();
    }

    Rcpp::List database(static_cast<int>(size));
    int counter = 0;


    for (const auto& node : nodeMap) {
        if (node.second.spectraList.empty())
            continue;
        for (const auto& spectra : node.second.spectraList) {
            database[counter++]  = Rcpp::List::create(
              Rcpp::Named("info") = Rcpp::List::create(
                   Rcpp::Named("key") = node.second.keys,
                   Rcpp::Named("value") = node.second.values),
              Rcpp::Named("spec") = Rcpp::List::create(
                   Rcpp::Named("mz") = spectra.mz,
                   Rcpp::Named("intensity") = spectra.intensity));
        }
    }
    return database;
}
