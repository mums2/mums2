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
    std::vector<double> results(rowSize);
    for(size_t i = 0; i < rowSize; i++) {
        std::vector<std::vector<double>> temp(1);
        temp[0] = communityMatrix.GetRow(i);
        results[i] = calculator->Calculate(temp);
    }
    delete calculator;
    return CppMatrix(results, rowSize, rowSize);
}
