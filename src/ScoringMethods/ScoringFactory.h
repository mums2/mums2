//
// Created by gregj on 10/2/2024.
//

#ifndef SCORINGFACTORY_H
#define SCORINGFACTORY_H
#include "Score.h"
#include <string>
#include <Rcpp.h>

class ScoringFactory {
public:
    explicit ScoringFactory(const Rcpp::List& parameters);
    double CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra, size_t minPeaks) const;
    ~ScoringFactory() {
        delete currentScoringAlgorithm;
    }
private:
    Score* currentScoringAlgorithm = nullptr;
    std::string scoringMethod;

};



#endif //SCORINGFACTORY_H
