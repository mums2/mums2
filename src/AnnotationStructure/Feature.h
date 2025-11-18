//
// Created by Gregory Johnson on 11/13/25.
//

#ifndef MUMS2_FEATURE_H
#define MUMS2_FEATURE_H
#include <string>

#include "../Spectra/Spectra.h"

struct Feature {
    std::string ms1_id, ms2_id, formula;
    double mz, rt;
    Spectra spectra;
};
#endif //MUMS2_FEATURE_H