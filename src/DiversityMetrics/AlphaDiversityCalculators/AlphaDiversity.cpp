//
// Created by Gregory Johnson on 1/16/25.
//

#include "AlphaDiversity.h"

#include "../DiversityMetricFactory.h"


Rcpp::NumericMatrix AlphaDiversity::CalculateDiversity(const Rcpp::NumericMatrix &communityMatrix,
                                                       const std::string &index) {
    const DiversityCalculator* calculator = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(index);
    const int size = communityMatrix.nrow();
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    std::vector<double> results(size);
    for(int i = 0; i < size; i++) {
        Rcpp::NumericVector abundance = communityMatrix(i, Rcpp::_);
        std::vector<double> diversities = Rcpp::as<std::vector<double>>(abundance);
        results[i] = calculator->Calculate({diversities});
    }
    delete calculator;


}
