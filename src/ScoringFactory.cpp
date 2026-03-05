//
// Created by gregj on 10/2/2024.
//

#include "ScoringMethods/ScoringFactory.h"
#include "ScoringMethods/SpectralEntropy/entropy.h"
#include "ScoringMethods/GNPS/ModifiedCosineScore.h"

ScoringFactory::ScoringFactory(const Rcpp::List &parameters) {
    scoringMethod = Rcpp::as<std::string>(parameters["method"]);
    if(scoringMethod == "gnps") {
        currentScoringAlgorithm = new ModifiedCosineScore(parameters);
    }
    else if(scoringMethod == "entropy") {
        currentScoringAlgorithm = new Entropy(parameters);
    }
}

double ScoringFactory::CalculateScore(const Spectra& firstSpectra, const Spectra& secondSpectra,
    const size_t minPeaks) const {
    // If there are not enough peaks, we return 0. As we cannot properly score
    if(firstSpectra.mz.size() < minPeaks || secondSpectra.mz.size() < minPeaks) return 0;
    return currentScoringAlgorithm->CalculateScore(firstSpectra, secondSpectra);
}
