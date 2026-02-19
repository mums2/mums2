//
// Created by gregj on 2/21/2024.
//

#ifndef GNPSSCORINGDYNAMICBST_H
#define GNPSSCORINGDYNAMICBST_H
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <Rcpp.h>
#include "../Score.h"
#include "../ScoreValues.h"
#include "../../NormalizeMs2/SquareRootNormalize.h"

class GNPSScoringDynamicPriorityQueue: public Score {
private:
    std::unordered_map<int, std::unordered_set<int>> ConstructPeaks(const std::vector<double>&, const std::vector<double>&,
        double, double, int&);
    std::vector<ScoreValues> ConstructPriorityQueue(std::unordered_map<int, std::unordered_set<int>>&, const std::vector<double>&,
        const std::vector<double>&, int);
    double ScoreMatches(std::vector<ScoreValues>&, size_t, int&);
    double tolerance;
public:
    GNPSScoringDynamicPriorityQueue() = default;
    explicit GNPSScoringDynamicPriorityQueue(const Rcpp::List&);
    std::vector<double> ScoreRData(const std::vector<double>&, std::vector<double>&, const std::vector<double>&, std::vector<double>&,
       double, double);

    double CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) override;
};



#endif // GNPSSCORINGDYNAMICBST_H
