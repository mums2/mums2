//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>

#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

// [[Rcpp::export]]
double CalculateDiversity(Rcpp::List abundanceList, const std::string& diversityIndex) {
    // Beta diversity requires two vectors
    // Alpha diversity requires one
    // Rarefaction requires one but also other parameters
    // Should I make a class to deal with the parametrization?
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
