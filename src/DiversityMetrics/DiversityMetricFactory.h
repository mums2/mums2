//
// Created by Gregory Johnson on 12/20/24.
//

#ifndef DIVERSITYMETRICFACTORY_H
#define DIVERSITYMETRICFACTORY_H
#include <string>

#include "Diversity.h"
#include "DiversityCalculator.h"


class DiversityMetricFactory final {
public:
    DiversityMetricFactory() = default;
    static DiversityCalculator* ChooseDiversityMetricBasedOnName(const std::string&);
    static Diversity* ChooseDiversityBasedOnIndex(const std::string&);
    ~DiversityMetricFactory() = default;


};



#endif //DIVERSITYMETRICFACTORY_H
