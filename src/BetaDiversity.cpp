//
// Created by Gregory Johnson on 1/16/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/BetaDiversity.h"

#include "DiversityMetrics/DiversityMetricFactory.h"

CppMatrix BetaDiversity::CalculateDiversity(const CppMatrix &communityMatrix,
                                                      const std::string &index) {
    const DiversityCalculator* calculator = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(index);
    const auto sampleSize = communityMatrix.GetRowSize();
    // const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rcpp::NumericMatrix brayCurtisMatrix(sampleSize, sampleSize);
    std::vector<double> brayCurtis(sampleSize * sampleSize);
    // Rcpp::List diversityList(2);
    std::vector<std::vector<double>> diversityList(2);
    for(size_t i = 0; i < sampleSize; i++) {
        // Rcpp::NumericVector abundance = communityMatrix(i, Rcpp::_);
        diversityList[0] = communityMatrix.GetRow(i);
        // std::vector<double> currentCommunity = Rcpp::as<std::vector<double>>(abundance);
        for(size_t j = i; j < sampleSize; j++) {
            if(j == i) continue;
            diversityList[1] = communityMatrix.GetRow(j);
            // const std::vector<double> otherCommunity = Rcpp::as<std::vector<double>>(tempCommunity);
            const double result = calculator->Calculate(diversityList);
            // brayCurtisMatrix(i,j) += result;
            // brayCurtisMatrix(j,i) += result;
            brayCurtis[i * sampleSize + j] = result;
            brayCurtis[j * sampleSize + i] = result;
        }
    }
    // Rcpp::colnames(brayCurtisMatrix) = samples;
    // Rcpp::rownames(brayCurtisMatrix) = samples;
    delete calculator;
    return CppMatrix(brayCurtis, sampleSize, sampleSize);
}
