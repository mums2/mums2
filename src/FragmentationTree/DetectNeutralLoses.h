//
// Created by Gregory Johnson on 5/7/26.
//

#ifndef MUMS2_DETECTNEUTRALLOSES_H
#define MUMS2_DETECTNEUTRALLOSES_H
#include <vector>

/// Common List of Neutral loses are from the paper below:
/// Sebastian Böcker, Florian Rasche, Towards de novo identification of metabolites by analyzing tandem mass spectra,
/// Bioinformatics, Volume 24, Issue 16, August 2008, Pages i49–i55,
/// https://doi.org/10.1093/bioinformatics/btn270
///
class DetectNeutralLoses {
public:
    DetectNeutralLoses();
    double DetermineNeutralLoses(const std::vector<int>& elements) const;
private:
    const double score = std::log(2);
    const double scoreRadical = std::log(0.25);
    double GetDBE(const std::vector<int> &elements);
    std::vector<std::vector<int>> neutralLosesList;
    std::vector<bool> isRadical;
};

#endif //MUMS2_DETECTNEUTRALLOSES_H