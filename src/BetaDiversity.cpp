//
// Created by Gregory Johnson on 1/16/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/BetaDiversity.h"

#include "DiversityMetrics/DiversityMetricFactory.h"

Rcpp::NumericMatrix BetaDiversity::CalculateDiversity(const Rcpp::NumericMatrix &communityMatrix,
                                                      const std::string &index) {
    const DiversityCalculator* calculator = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(index);
    const auto sampleSize = static_cast<size_t>(communityMatrix.nrow());
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rcpp::NumericMatrix brayCurtisMatrix(sampleSize, sampleSize);
    Rcpp::List diversityList(2);
    for(size_t i = 0; i < sampleSize; i++) {
        // Rcpp::NumericVector abundance = communityMatrix(i, Rcpp::_);
        diversityList[0] = communityMatrix(i, Rcpp::_);
        // std::vector<double> currentCommunity = Rcpp::as<std::vector<double>>(abundance);
        for(size_t j = i; j < sampleSize; j++) {
            if(j == i) continue;
            diversityList[1] = communityMatrix(j, Rcpp::_);
            // const std::vector<double> otherCommunity = Rcpp::as<std::vector<double>>(tempCommunity);
            const double result = calculator->Calculate(diversityList);
            brayCurtisMatrix(i,j) += result;
            brayCurtisMatrix(j,i) += result;
        }
    }
    Rcpp::colnames(brayCurtisMatrix) = samples;
    Rcpp::rownames(brayCurtisMatrix) = samples;
    delete calculator;
    return brayCurtisMatrix;
}
