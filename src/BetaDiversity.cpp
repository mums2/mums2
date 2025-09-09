//
// Created by Gregory Johnson on 1/16/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/BetaDiversity.h"

#include "DiversityMetrics/DiversityMetricFactory.h"

CppMatrix BetaDiversity::CalculateDiversity(const CppMatrix &communityMatrix,
                                                      const std::string &index) {
    const DiversityCalculator* calculator = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(index);
    const auto sampleSize = communityMatrix.GetRowSize();
    std::vector<double> betaDiversityMatrix(sampleSize * sampleSize);
    std::vector<std::vector<double>> diversityList(2);
    for(size_t i = 0; i < sampleSize; i++) {
        diversityList[0] = communityMatrix.GetRow(i);
        for(size_t j = i; j < sampleSize; j++) {
            if(j == i) continue;
            diversityList[1] = communityMatrix.GetRow(j);
            const double result = calculator->Calculate(diversityList);
            betaDiversityMatrix[i * sampleSize + j] = result;
            betaDiversityMatrix[j * sampleSize + i] = result;
        }
    }
    delete calculator;
    return CppMatrix(betaDiversityMatrix, sampleSize, sampleSize);
}
