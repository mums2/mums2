//
// Created by Gregory Johnson on 1/16/25.
//

#include "DiversityMetrics/AlphaDiversityCalculators/AlphaDiversity.h"
#include "DiversityMetrics/DiversityMetricFactory.h"


Rcpp::NumericMatrix AlphaDiversity::CalculateDiversity(const Rcpp::NumericMatrix &communityMatrix,
                                                       const std::string &index) {
    const DiversityCalculator* calculator = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(index);
    const int size = communityMatrix.nrow();
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rcpp::NumericMatrix results(1, size);
    for(int i = 0; i < size; i++) {
        Rcpp::NumericVector abundance = communityMatrix(i, Rcpp::_);
        std::vector<double> diversities = Rcpp::as<std::vector<double>>(abundance);
        results(0, i) = calculator->Calculate({diversities});
    }
    Rcpp::colnames(results) = samples;
    delete calculator;
    return results;
}
