
#include "Rarefy/rarefyMs.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <map>
#include <unordered_map>
#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"
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

size_t GetRandomNumberIndex(const std::vector<int>& weightedToPull, const int sum) {
    auto randomNumber = static_cast<int>(R::runif(0, static_cast<double>(sum)));
    for (size_t i = 0; i < weightedToPull.size(); i++) {
        if (randomNumber <= weightedToPull[i]) {
            return i;
        }
        randomNumber -= weightedToPull[i];
    }
    return 0;
}


// [[Rcpp::export]]
Rcpp::DataFrame rarefyMs_2(const std::vector<int>& feature, const std::vector<int>& abund,
    const int size, const int threshold) {
    const auto startTime = std::chrono::system_clock::now();
    const int sum = std::accumulate(abund.begin(), abund.end(), 0);
    int grandTotal = 0;
    int incrementer = 0;
    std::unordered_map<int, int> filtered;
    while(grandTotal <= size) {
        const auto startTimeLoop = std::chrono::system_clock::now();
        std::unordered_map<int, int> finalMap;
        std::unordered_map<int, int> counter;
        std::vector<int> abundanceCopy = abund;
        const int currentSize = size + incrementer;
        const auto firstLoop = std::chrono::system_clock::now();
        for(int i = 0; i < currentSize; i++) {
            const size_t index = GetRandomNumberIndex(abundanceCopy, sum - i);
            abundanceCopy[index]--;
            counter[feature[index]]++;
        }
        const auto secondTimeLoop = std::chrono::system_clock::now();
        std::chrono::duration<double> currentTime = secondTimeLoop - firstLoop;
        Rcpp::Rcout << "First Loop Timing: " << currentTime.count() << std::endl;
        const auto secondLoop = std::chrono::system_clock::now();
        for(const auto& values: counter) {
            if(values.second < threshold)
                continue;

            grandTotal += values.second;
            finalMap[values.first] = values.second;
        }
        const auto secondLoopEndTime = std::chrono::system_clock::now();
        currentTime = secondLoopEndTime - secondLoop;
        Rcpp::Rcout << "Second Loop Timing: " << currentTime.count() << std::endl;
        if(grandTotal >= size) {
            filtered = finalMap;
            break;
        }
        incrementer++;
        const auto endTime = std::chrono::system_clock::now();
        currentTime = endTime - startTimeLoop;
        Rcpp::Rcout << "Loop Timing: " << currentTime.count() << std::endl;
    }

    std::vector<int> rare_mz(filtered.size(), 0);
    std::vector<int> rare_abund(filtered.size(), 0);
    int pos = 0;
    for (const auto & it : filtered) {
        rare_mz[pos] = it.first;
        rare_abund[pos] = it.second;
        pos ++;
    }

    const auto endTime = std::chrono::system_clock::now();
    const std::chrono::duration<double> currentTime = endTime - startTime;
    Rcpp::Rcout << "Rarefy Time: " << currentTime.count() << std::endl;
    // delete calculator;
    // This is where we compute alpha and beta diversity
    return DataFrame::create(Named("mz") = rare_mz,
                             _["abund"] = rare_abund);
}

// [[Rcpp::export]]
void Test(const std::vector<int>& abund) {
    constexpr auto size = 25011;
    const auto startTime = std::chrono::system_clock::now();
    const int sum = std::accumulate(abund.begin(), abund.end(), 0);
    //auto vec = std::vector<size_t>(size);
    for(int i = 0; i < 25011; i++) {
        GetRandomNumberIndex(abund, sum - i);
    }
    const auto endTime = std::chrono::system_clock::now();
    const std::chrono::duration<double> currentTime = endTime - startTime;
    Rcpp::Rcout << "Time: " << currentTime.count() << std::endl;
}


// [[Rcpp::export]]
Rcpp::DataFrame rarefyMs_3(const std::vector<std::string>& feature, const std::vector<int>& abund,
    const int size, const int threshold) {
    int sum = 0;
    for (const auto& value : abund) {
        sum += value;
    }
    int grandTotal = 0;
    int incrementer = 0;
    std::unordered_map<std::string, int> filtered;
    while(grandTotal <= size) {
        std::unordered_map<std::string, int> finalMap;
        std::unordered_map<std::string, int> counter;
        std::vector<int> abundanceCopy = abund;
        const int currentSize = size + incrementer;
        for(int i = 0; i < currentSize; i++) {
            const size_t index = GetRandomNumberIndex(abundanceCopy, sum - i);
            abundanceCopy[index]--;
            counter[feature[index]]++;
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

    std::vector<std::string> rare_mz(filtered.size(), "");
    std::vector<int> rare_abund(filtered.size(), 0);
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
double CalculateAlphaDiverstiy(const std::vector<std::string>& feature, const std::vector<int>& abund,
    const int size, const int threshold, const int iterations = 1000) {
    DiversityCalculator* calculator = new SimpsonsDiversityIndex();
    double simpsonDiversity = 0;
    for(int i = 0; i < iterations; i++) {
        const auto df = rarefyMs_3(feature, abund, size, threshold);
        // const std::vector<int> rare_abund = df["abund"];
        // simpsonDiversity += calculator->Calculate(std::vector<double>(rare_abund.begin(), rare_abund.end()));
    }
    delete calculator;
    return simpsonDiversity/iterations;
}

// [[Rcpp::export]]
double CalculateAlphaDiversityInt(const std::vector<int>& feature, const std::vector<int>& abund,
    const int size, const int threshold, const int iterations = 1000) {
    const DiversityCalculator* calculator = new SimpsonsDiversityIndex();
    double simpsonDiversity = 0;
    for(int i = 0; i < iterations; i++) {
        const auto df = rarefyMs_2(feature, abund, size, threshold);
        const std::vector<int> rare_abund = df["abund"];
        simpsonDiversity += calculator->Calculate(std::vector<double>(rare_abund.begin(), rare_abund.end()));
    }
    delete calculator;
    return simpsonDiversity/iterations;
}

// [[Rcpp::export]]
double CalculateAlphaDiversityShannon(const std::vector<int>& feature, const std::vector<int>& abund,
    const int size, const int threshold, const int iterations = 1000) {
    const DiversityCalculator* calculator = new ShannonDiversityIndex();
    double diversity = 0;
    for(int i = 0; i < iterations; i++) {
        const auto df = rarefyMs_2(feature, abund, size, threshold);
        const std::vector<int> rare_abund = df["abund"];
        diversity += calculator->Calculate(std::vector<double>(rare_abund.begin(), rare_abund.end()));
    }
    delete calculator;
    return diversity/iterations;
}

