
#include "Rarefy/rarefyMs.h"
#include <algorithm>
#include <random>
#include <map>
#include <unordered_map>
#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"
#include "DiversityMetrics/AlphaDiversityCalculators/SimpsonsDiversityIndex.h"
#include "Utils/Utils.h"
#include "Rarefy/AbundanceMap.h"
#include <Rcpp.h>
#include <chrono>
#include "DiversityMetrics/BetaDiversityCalculators/BrayCurtisDissimilarity.h"
// #include <vector>
// using namespace Rcpp;

size_t GetRandomNumberIndex(const std::vector<int64_t>& weightedToPull, const int64_t sum) {
    auto randomNumber = static_cast<int64_t>(R::runif(0, static_cast<double>(sum)));
    for (size_t i = 0; i < weightedToPull.size(); i++) {
        if (randomNumber < weightedToPull[i]) {
            return i;
        }
        randomNumber -= weightedToPull[i];
    }
    return 0;
}

// [[Rcpp::export]]
Rcpp::DataFrame rarefyMs(const std::vector<int>& feature, std::vector<int64_t>& abund,
    const int64_t size, const int threshold) {
    const int64_t sum = std::accumulate(abund.begin(), abund.end(), 0LL);
    int64_t grandTotal = 0;
    int64_t incrementer = size;
    const size_t abundSize = abund.size();
    std::unordered_map<int, int64_t> filtered;
    std::vector<int> counter(abundSize, 0);
    while(grandTotal <= size) {
        std::unordered_map<int, int64_t> finalMap;
        const int64_t currentSize = incrementer;
        for(int64_t i = 0; i < currentSize; i++) {
            const size_t index = GetRandomNumberIndex(abund, sum - i);
            abund[index]--;
            counter[index]++;
        }
        int count = 0;
        for(const auto& feat: feature) {
            const int value = counter[count++];
            if(value < threshold)
                continue;
            grandTotal += value;
            finalMap[feat] = value;
        }
        if(grandTotal >= size) {
            filtered = finalMap;
            break;
        }
        incrementer = size - grandTotal;
        grandTotal = 0;
    }
    std::vector<int> rare_mz(filtered.size(), 0);
    std::vector<int64_t> rare_abund(filtered.size(), 0);
    int pos = 0;
    for (const auto & it : filtered) {
        rare_mz[pos] = it.first;
        rare_abund[pos] = it.second;
        pos ++;
    }

    // delete calculator;
    // This is where we compute alpha and beta diversity
    return DataFrame::create(Named("mz") = rare_mz,
                             _["abund"] = rare_abund);
}

// [[Rcpp::export]]
double CalculateAlphaDiversitySimpson(const std::vector<int>& feature, const std::vector<int64_t>& abund,
    const int size, const int threshold, const int iterations = 1000) {
    const DiversityCalculator* calculator = new SimpsonsDiversityIndex();
    double simpsonDiversity = 0;
    for(int i = 0; i < iterations; i++) {
        auto abundCopy = std::vector<int64_t>(abund.begin(), abund.end());
        const auto df = rarefyMs(feature,abundCopy, size, threshold);
        const std::vector<int> rare_abund = df["abund"];
        simpsonDiversity += calculator->Calculate({std::vector<double>(rare_abund.begin(),
            rare_abund.end())});
    }
    delete calculator;
    return simpsonDiversity/iterations;
}

// [[Rcpp::export]]
double CalculateAlphaDiversityShannon(const std::vector<int>& feature, std::vector<int64_t>& abund,
    const int size, const int threshold, const int iterations = 1000) {
    const DiversityCalculator* calculator = new ShannonDiversityIndex();
    double diversity = 0;
    for(int i = 0; i < iterations; i++) {
        auto abundCopy = std::vector<int64_t>(abund.begin(), abund.end());
        const auto df = rarefyMs(feature,abundCopy, size, threshold);
        const std::vector<int> rare_abund = df["abund"];
        diversity += calculator->Calculate({std::vector<double>(rare_abund.begin(),
            rare_abund.end())});
    }
    delete calculator;
    return diversity/iterations;
}

// [[Rcpp::export]]
NumericMatrix CalculateBrayCurtisDissimilarity(const Rcpp::List &features, Rcpp::List& abund,
    const int size, const int threshold, const int iterations = 1000) {
    const DiversityCalculator* calculator = new BrayCurtisDissimilarity();

    const size_t abundSize = abund.size();
    NumericMatrix brayCurtisMatrix(abundSize, abundSize);
    for(int i = 0; i < iterations; i++) {
        std::vector<std::vector<double>> rarefyAbundanceVector(abundSize, std::vector<double>());
        for(size_t j = 0; j < abundSize; j++) {
            std::vector<int> feat = features[j];
            std::vector<int64_t> ab = abund[j];
            const Rcpp::DataFrame val = rarefyMs(feat,
                ab, size, threshold);
            rarefyAbundanceVector[j] = as<std::vector<double>>(val["abund"]);
        }
        for(size_t j = 0; j < abundSize; j++) {
            /* TODO: Create a list of dataframes then use those to calculate bray between each of them.
            // Make this work for all of them using the Rcpp::list notation
            // We have to ensure that we rarefy each one of the samples and then compare them all*/
            const std::vector<double>& currentAbundance = rarefyAbundanceVector[j];
            std::vector<double> brayResults = std::vector<double>(abundSize, 0);
            for(size_t k = j; k < abundSize; k++) {
                if(j == k) continue;
                const std::vector<double>& otherAbundance = rarefyAbundanceVector[k];
                    const double result = calculator->Calculate({std::vector<double>(currentAbundance.begin(),
                        currentAbundance.end()),
                            std::vector<double>(otherAbundance.begin(),otherAbundance.end())});
                brayCurtisMatrix(j,k) += result;
                brayCurtisMatrix(k,j) += result;
            }
        }
    }
    delete calculator;
    return brayCurtisMatrix/iterations;
}
