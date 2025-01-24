//
// Created by gregj on 1/14/2025.
//
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"
#include "Random/ReservoirPairs.h"
#include <complex>
#include <algorithm>
size_t RandomizationMethods::GetRandomNumberIndex(const std::vector<int64_t> &weightedToPull, const size_t vectorSize,
    const int64_t sum) {
    auto randomNumber = static_cast<int64_t>(R::runif(0, static_cast<double>(sum)));
    for (size_t i = 0; i < vectorSize; i++) {
        const int64_t currentWeight = weightedToPull[i];
        if (randomNumber < currentWeight) {
            return i;
        }
        randomNumber -= currentWeight;
    }
    return 0;
}

std::vector<size_t> RandomizationMethods::ReservoirSampling(const int V, const size_t m, const std::vector<int64_t>& weights) {
    // Population = V
    // weighted items
    // m is the sample size

    // Step one the first m items are inserted into R
    const size_t size = weights.size();
    std::vector<ReservoirPairs> reservoirPairs(size);
    // Step two: for each item calculate the key, k_i = Ui^(1/weight[i])
    // Ui = random(0,1)
    for(size_t i = 0; i < size; i++) {
        const double u_i = R::runif(0, 1);
        reservoirPairs[i] = {std::pow(u_i, 1.0 / weights[i]), i};
    }
    std::make_heap(reservoirPairs.begin(), reservoirPairs.end(), CompareReservoir());
    // Step three: minimum threshold, T_w, the min key
    double minKey = reservoirPairs.front().key;
    // Step 4 Repeat step 5 - 10 until exhausted
    for(size_t i = m; i < V; i++){
        // Step 5 Let r = random(0,1) and Xw = log(r)/ log(T_w)
        const double r = R::runif(0, 1);
        const double x_w = std::log(r) / std::log(minKey);
        // x_w is a weight constant
        // Step 6 & 7, from current item, vc, skip items until item vi,
        // such that wc + w(c+1) + w(c + 2) + ... w(c + n) < X_w
        // X_w <= wc + w(c + 1) + w(c + 2) + ... w(c + n)
        double weight_wc = 0;
        for(size_t j = i; j < V; j++) {
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

    std::vector<size_t> res(m);
    int count = 0;
    while(!reservoirPairs.empty() && count < m) {
        res[count++] = reservoirPairs.back().value;
        reservoirPairs.pop_back();
    }
    return res;
}
