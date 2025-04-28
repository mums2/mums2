//
// Created by gregj on 10/2/2024.
//

#include "ScoringMethods/ScoringFactory.h"

#include <utility>

#include "ScoringMethods/SpectralEntropy/entropy.h"
#include "ScoringMethods/GNPS/GNPSScoringDynamicPriorityQueue.h"

ScoringFactory::ScoringFactory(const Rcpp::List &parameters) {
    scoringMethod = Rcpp::as<std::string>(parameters["method"]);
    if(scoringMethod == "gnps") {
        currentScoringAlgorithm = new GNPSScoringDynamicPriorityQueue(parameters);
    }
    else if(scoringMethod == "entropy") {
        currentScoringAlgorithm = new Entropy(parameters);
    }
}

double ScoringFactory::CalculateScore(const Spectra& firstSpectra, const Spectra& secondSpectra) const {
    return currentScoringAlgorithm->CalculateScore(firstSpectra, secondSpectra);
}
