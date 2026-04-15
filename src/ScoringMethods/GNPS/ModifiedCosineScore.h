//
// Created by gregj on 2/21/2024.
//

#ifndef MODIFEDCOSINESCORE_H
#define MODIFEDCOSINESCORE_H
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <Rcpp.h>
#include "../Score.h"
#include "../ScoreValues.h"

class ModifiedCosineScore final : public Score {
private:
    std::vector<ScoreValues> ConstructPeaks(const std::vector<double>&, const std::vector<double>&,
        const std::vector<double>&, const std::vector<double>&, double, double);
    std::vector<ScoreValues> ConstructPriorityQueue(std::unordered_map<size_t, std::unordered_set<size_t>>&, const std::vector<double>&,
        const std::vector<double>&, size_t);
    double ScoreMatches(std::vector<ScoreValues>&, size_t, size_t&);
    double tolerance = 0;
public:
    ModifiedCosineScore() = default;
    explicit ModifiedCosineScore(const Rcpp::List&);
    std::vector<double> ScoreRData(const std::vector<double>&, std::vector<double>&, const std::vector<double>&, std::vector<double>&,
       double, double);

    double CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) override;
};



#endif // MODIFEDCOSINESCORE_H
