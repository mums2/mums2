//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <algorithm>
#include <complex>
#include <cstdint>
#include <numeric>
#include "DiversityMetrics/Diversity.h"
// #include "../../../../Downloads/gperftools-2.15/src/gperftools/profiler.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& abundances, const std::string& diversityIndex) {
    std::string index = diversityIndex;
    std::transform(index.begin(), index.end(), index.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    Rcpp::NumericMatrix results = diversity->CalculateDiversity(abundances, index);
    delete diversity;
    return results;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const Rcpp::NumericMatrix& communityMatrix, const int64_t size,
    const int64_t threshold) {
    const int row = communityMatrix.nrow();
    const int col = communityMatrix.ncol();

    std::vector<int> indexToName(col);
    std::iota(indexToName.begin(), indexToName.end(), 0);
    Rarefaction rarefaction;
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
        // Rcpp::Rcout << communityVector << std::endl;
        // Rarefy
        // We are going to have to switch between the transpose of the community matrix and
        // the original matrix.
        // Samples are represented by rows and columns represent species.
        // So ex..
        //          FD39  FD09 <- is backwards            species1     <- Correct
        // Species1  0      1                       FD39    1
        //                                          DF09    0
        const auto results = rarefaction.Rarefy(indexToName, communityVector, size, threshold);
        for(int j = 0; j < col; j++) {
            resultantMatrix(i, j) = results[j];
        }
    }
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const Rcpp::NumericMatrix& communityMatrix, const std::string& index,
    const int64_t size, const int64_t threshold, const int iterations = 1000) {
    Rcpp::NumericMatrix diversity = CalculateDiversity(RarefactionCalculation(communityMatrix,
        size, threshold), index);
    for(int i = 1; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = RarefactionCalculation(communityMatrix, size, threshold);
        diversity += CalculateDiversity(rarefyMatrix, index);
    }
    return diversity/iterations;
}


// Test Sampling algorithms:

void AlgorithmOne(int n, int s, const std::vector<int>& weights, bool replace = true) {
    if(s <= 0)
        return;

}

struct ReservoirPairs {
    double key;
    int value;
};
struct CompareReservoir {
    bool operator()(const ReservoirPairs& first, const ReservoirPairs& other) const {
        return first.key > other.key;
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector SomePaper(int V, int m, const std::vector<double>& weights) {
    // Population = V
    // weighted items
    // m is the sample size

    // Step one the first m items are inserted into R
    const size_t size = weights.size();
    std::vector<ReservoirPairs> reservoirPairs(size);
    auto it = reservoirPairs.begin();
    ++it;
    // Step two: for each item calculate the key, k_i = Ui^(1/weight[i])
    // Ui = random(0,1)
    for(int i = 0; i < size; i++) {
        const double u_i = R::runif(0, 1);
        reservoirPairs[i] = {std::pow(u_i, 1.0 / weights[i]), i};
    }
    std::make_heap(reservoirPairs.begin(), reservoirPairs.end(), CompareReservoir());
    // Step three: minimum threshold, T_w, the min key
    double minKey = reservoirPairs.front().key;
    // Step 4 Repeat step 5 - 10 until exhausted
    for(int i = m; i < V; i++){
        // Step 5 Let r = random(0,1) and Xw = log(r)/ log(T_w)
        const double r = R::runif(0, 1);
        const double x_w = std::log(r) / std::log(minKey);
        // x_w is a weight constant
        // Step 6 & 7, from current item, vc, skip items until item vi,
        // such that wc + w(c+1) + w(c + 2) + ... w(c + n) < X_w
        // X_w <= wc + w(c + 1) + w(c + 2) + ... w(c + n)
        double weight_wc = 0;
        for(int j = i; j < V; j++) {
            weight_wc += weights[j];
            if(weight_wc > x_w)
                break;
            i++;
        }

        if(i > V) break;
        // Step 9, lets t_w = Tw^(w_i)
        // r_2 = random (t_w, 1)
        // v_i(key): ki = (r_2)^(1/w_i)
        // T_w is the current minkey
        //r_2 is the second random value
        const double t_w = std::pow(minKey, 1.0 / weights[i]);
        const double r_2 = R::runif(t_w, 1);
        const double k_i = std::pow(r_2, 1.0 / weights[i]);
        // k_i is the new minKey
        // Step 8, item in R is replaced with item v_i
        std::pop_heap(reservoirPairs.begin(), reservoirPairs.end(), CompareReservoir());
        reservoirPairs.pop_back();
        reservoirPairs.push_back({k_i, i});
        std::push_heap(reservoirPairs.begin(), reservoirPairs.end(), CompareReservoir());
        minKey = reservoirPairs.front().key;
    }

    Rcpp::NumericVector vals(static_cast<int>(m));
    int count = 0;
    while(!reservoirPairs.empty() && count < m) {
        vals[count++] = reservoirPairs.back().value;
        reservoirPairs.pop_back();
    }
    return vals;
}

// [[Rcpp::export]]
size_t GetRandomNumberIndex(const std::vector<double> &weightedToPull, const size_t vectorSize,
    const double sum) {
    auto randomNumber = static_cast<double>(R::runif(0, static_cast<double>(sum)));
    for (size_t i = 0; i < vectorSize; i++) {
        const double currentWeight = weightedToPull[i];
        if (randomNumber < currentWeight) {
            return i;
        }
        randomNumber -= currentWeight;
    }
    return 0;
}
