//
// Created by Gregory Johnson on 1/16/25.
//

#include "DiversityMetrics/AlphaDiversityCalculators/AlphaDiversity.h"

#include "DataStructures/CppMatrix.h"
#include "DiversityMetrics/DiversityMetricFactory.h"


CppMatrix AlphaDiversity::CalculateDiversity(const CppMatrix &communityMatrix,
                                                       const std::string &index) {
    const DiversityCalculator* calculator = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(index);
    const size_t rowSize = communityMatrix.GetRowSize();
    // const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    // Rcpp::NumericMatrix results(1, size);
    std::vector<double> results(rowSize);
    for(int i = 0; i < rowSize; i++) {
        std::vector<std::vector<double>> temp(1);
        temp[0] = communityMatrix.GetRow(i);
        results[i] = calculator->Calculate(temp);
    }
    // Rcpp::colnames(results) = samples;
    delete calculator;
    return CppMatrix(results, 1, rowSize);
}
