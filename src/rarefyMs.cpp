
#include "Rarefy/rarefyMs.h"
#include <algorithm>
#include <random>
#include <map>

#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/AlphaDiversityCalculators/SimpsonsDiversityIndex.h"
#include "Utils/Utils.h"
#include "Rarefy/AbundanceMap.h"
// #include <Rcpp.h>
// #include <vector>
// using namespace Rcpp;

void RandomShuffle(std::vector<int>& randomize) {
    Rcpp::IntegerVector randomValues = Rcpp::wrap(randomize);
    const int size = static_cast<int>(randomize.size());
    randomValues = Rcpp::sample(randomValues, size);
    randomize = Rcpp::as<std::vector<int>>(randomValues);
}

// [[Rcpp::export]]
DataFrame rarefyMs(IntegerVector feature, IntegerVector abund, int size, int threshold) {
    size -= 1;

    std::vector<int> feature1 = Rcpp::as<std::vector<int>>(feature);
    std::vector<int> abund1 = Rcpp::as<std::vector<int>>(abund); 

    Utils util;
    std::vector<int> pool = util.my_rep(feature1, abund1);

    RandomShuffle(pool);

    int x = size;
    // bool done = false;
    int grand_total = 0;
    
    std::map<int, int> filtered;
    for (int i = 0; i < size; i++) {
        auto it = filtered.find(pool[i]);
        
        if (it == filtered.end()) {
            filtered[pool[i]] = 1;
            } else {
                it->second++;
            }
        }
   
    // std::unorderd_set<int> uniqueFeatures;
    x += 1;
    while(grand_total <= size) {
        grand_total = 0;
 
        // correctedFilter.clear();
        // for (int i = x; i < x; i++) {
        //     auto it = filtered.find(pool[i]);
        //     if (it == filtered.end()) {
        //         filtered[pool[i]] = 1;
        //     } else {
        //         it->second++;
        //     }
        // }
        filtered[pool[x++]]++;
        for (auto it = filtered.begin(); it != filtered.end();) {
            if (it->second < threshold) {
                filtered.erase(it++);
            } else {
                grand_total += it->second;
                it++;
            }
        }
        // std::cout << "grand total: " << grand_total << std::endl;
        // std::cout << "x: " << x << std::endl;
 }


    std::vector<int> rare_mz(filtered.size(), 0);
    std::vector<int> rare_abund(filtered.size(), 0);
    int pos = 0;

    for (auto it = filtered.begin(); it != filtered.end(); it++) {
        rare_mz[pos] = it->first;
        rare_abund[pos] = it->second;
        pos ++;
    }
    
    return DataFrame::create(Named("mz") = rare_mz,
                             _["abund"] = rare_abund);
}

// [[Rcpp::export]]
Rcpp::DataFrame rarefyMs_2(const std::vector<int>& feature, const std::vector<int>& abund, int size, int threshold) {
    const IntegerVector featureSugar = Rcpp::wrap(feature);

    const int sizeOfPool = static_cast<int>(feature.size());
    std::vector<AbundanceMap> pool(sizeOfPool);
    std::vector<double> percentageToPull(sizeOfPool);
    int currentMaxAbundance = 0;
    for(int i = 0; i < sizeOfPool; i++) {
        currentMaxAbundance += abund[i];
    }
    for(int i = 0; i < sizeOfPool; i++) {
        const int abundance = abund[i];
        pool[i] = AbundanceMap{feature[i], abundance};
        percentageToPull[i] = static_cast<double>(abundance)/currentMaxAbundance;

    }
    const NumericVector percentage = Rcpp::wrap(percentageToPull);
    int grandTotal = 0;
    int incrementer = 0;
    std::unordered_map<int, int> filtered;
    while(grandTotal <= size) {
        std::unordered_map<int, int> finalMap;
        std::unordered_map<int, int> counter;
        for(int i = 0; i < size+incrementer; i++) {
            // const double randomNumber = R::runif(0, 1);
            // const int index = CloserIndex(percentageToPull, randomNumber);
            int value = sample(featureSugar, 1, false, percentage)[0];
            counter[value]++;
        }
        for(const auto& values: counter) {
            if(values.second < threshold)
                continue;

            grandTotal += values.second;
            finalMap[values.first] = values.second;
        }

        if(grandTotal >= size) {
            filtered = finalMap;
            break;
        }
        incrementer++;
    }

    std::vector<int> rare_mz(filtered.size(), 0);
    std::vector<int> rare_abund(filtered.size(), 0);
    int pos = 0;
    for (const auto & it : filtered) {
        rare_mz[pos] = it.first;
        rare_abund[pos] = it.second;
        pos ++;
    }
    DiversityCalculator* calculator = new SimpsonsDiversityIndex();
    const double simpsonDiversity = calculator->Calculate(std::vector<double>(rare_abund.begin(), rare_abund.end()));
    Rcpp::Rcout << "The Simpson's Diversity is: " << simpsonDiversity << std::endl;
    // delete calculator;
    // This is where we compute alpha and beta diversity
    return DataFrame::create(Named("mz") = rare_mz,
                             _["abund"] = rare_abund);
}

Rcpp::DataFrame rarefyMs_3(const std::vector<int>& feature, const std::vector<int>& abund, int size, int threshold) {
    const IntegerVector featureSugar = Rcpp::wrap(feature);

    const int sizeOfPool = static_cast<int>(feature.size());
    std::vector<AbundanceMap> pool(sizeOfPool);
    std::vector<double> percentageToPull(sizeOfPool);
    int currentMaxAbundance = 0;
    for(int i = 0; i < sizeOfPool; i++) {
        currentMaxAbundance += abund[i];
    }
    for(int i = 0; i < sizeOfPool; i++) {
        const int abundance = abund[i];
        pool[i] = AbundanceMap{feature[i], abundance};
        percentageToPull[i] = static_cast<double>(abundance)/currentMaxAbundance;

    }
    const NumericVector percentage = Rcpp::wrap(percentageToPull);
    int grandTotal = 0;
    int incrementer = 0;
    std::unordered_map<int, int> filtered;
    while(grandTotal <= size) {
        std::unordered_map<int, int> finalMap;
        std::unordered_map<int, int> counter;
        for(int i = 0; i < size+incrementer; i++) {
            // const double randomNumber = R::runif(0, 1);
            // const int index = CloserIndex(percentageToPull, randomNumber);
            int value = sample(featureSugar, 1, false, percentage)[0];
            counter[value]++;
        }
        for(const auto& values: counter) {
            if(values.second < threshold)
                continue;

            grandTotal += values.second;
            finalMap[values.first] = values.second;
        }

        if(grandTotal >= size) {
            filtered = finalMap;
            break;
        }
        incrementer++;
    }

    std::vector<int> rare_mz(filtered.size(), 0);
    std::vector<int> rare_abund(filtered.size(), 0);
    int pos = 0;
    for (const auto & it : filtered) {
        rare_mz[pos] = it.first;
        rare_abund[pos] = it.second;
        pos ++;
    }
    DiversityCalculator* calculator = new SimpsonsDiversityIndex();
    const double simpsonDiversity = calculator->Calculate(std::vector<double>(rare_abund.begin(), rare_abund.end()));
    Rcpp::Rcout << "The Simpson's Diversity is: " << simpsonDiversity << std::endl;
    // delete calculator;
    // This is where we compute alpha and beta diversity
    return DataFrame::create(Named("mz") = rare_mz,
                             _["abund"] = rare_abund);
}