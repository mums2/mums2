//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <algorithm>

#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

// [[Rcpp::export]]
double CalculateDiversity(Rcpp::List abundanceList, std::string& diversityIndex) {
    // Beta diversity requires two vectors
    // Alpha diversity requires one
    // Rarefaction requires one but also other parameters
    // Should I make a class to deal with the parametrization?
    std::transform(diversityIndex.begin(), diversityIndex.end(), diversityIndex.begin(), tolower);
    const DiversityCalculator* diversity = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(diversityIndex);
    const int size = abundanceList.size();
    std::vector<std::vector<double>> listOfAbundances(abundanceList.size(), std::vector<double>());
    if(size > 1) {
        for(int i = 0; i < size; i++) {
            listOfAbundances[i] = Rcpp::as<std::vector<double>>(abundanceList[i]);
        }
    }
    else {
        listOfAbundances[0] = Rcpp::as<std::vector<double>>(abundanceList[0]);
    }
    return diversity->Calculate(listOfAbundances);
}

// [[Rcpp::export]]
Rcpp::DataFrame Rarefaction(const Rcpp::NumericMatrix& communityMatrix, const int64_t size,
    const int64_t threshold) {

    for(int i = 0; i < communityMatrix.nrow(); i++) {
        Rcpp::NumericVector community = communityMatrix[i];
        std::vector<double> communityVector = Rcpp::as<std::vector<double>>(community);
        // Rarefy
        // Then compute the alpha diversity
    }
    return {};
}