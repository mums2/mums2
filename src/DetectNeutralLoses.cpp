//
// Created by Gregory Johnson on 5/7/26.
//

#include "FragmentationTree/DetectNeutralLoses.h"
#include <algorithm>
DetectNeutralLoses::DetectNeutralLoses() {
    constexpr size_t size = 17;
    neutralLosesList = std::vector<std::vector<int>>(size);
    const std::vector<int> Methyl{1, 3, 0, 0, 0, 0}; // CH3
    const std::vector<int> Methane{1, 4, 0, 0, 0, 0}; // CH4
    const std::vector<int> Oxy{1, 0, 0, 0, 0, 0}; // O
    const std::vector<int> Hydroxyl{0, 2, 0, 1, 0, 0}; // H20
    const std::vector<int> CarbonMonoxide{1, 0, 0, 1, 0, 0}; // CO
    const std::vector<int> Nitrogen{0, 0, 2, 0, 0, 0}; // N2
    const std::vector<int> Ammonia{0, 3, 1, 0, 0, 0}; // NH3
    const std::vector<int> Ethyl{2, 4, 0, 0, 0, 0}; // C2H4
    const std::vector<int> Formaldehyde{1, 2, 0, 1, 0, 0}; //CH2O
    const std::vector<int> Isobutene{4, 8, 0, 0, 0, 0}; // C4H8
    const std::vector<int> Isopentene{5, 8, 0, 0, 0, 0}; // C5H8
    const std::vector<int> FormicAcid{1, 2, 0, 2, 0, 0}; // CH2O2
    const std::vector<int> MalonicAcid{3, 2, 0, 3, 0, 0}; // C3H2O3
    const std::vector<int> Xylose{5, 8, 0, 4, 0, 0}; // C5H8O4
    const std::vector<int> Rhamnose{6, 10, 0, 4, 0, 0}; // C6H10O4
    const std::vector<int> Hexose{6, 10, 0, 5, 0, 0}; // C6H10O5
    const std::vector<int> GlucuronicAcid{6, 8, 0, 6, 0, 0}; // C6H8O6

    neutralLosesList[0] = Methyl;
    neutralLosesList[1] = Methane;
    neutralLosesList[2] = Oxy;
    neutralLosesList[3] = Hydroxyl;
    neutralLosesList[4] = CarbonMonoxide;
    neutralLosesList[5] = Nitrogen;
    neutralLosesList[6] = Ammonia;
    neutralLosesList[7] = Ethyl;
    neutralLosesList[8] = Formaldehyde;
    neutralLosesList[9] = Isobutene;
    neutralLosesList[10] = Isopentene;
    neutralLosesList[11] = FormicAcid;
    neutralLosesList[12] = MalonicAcid;
    neutralLosesList[13] = Xylose;
    neutralLosesList[14] = Rhamnose;
    neutralLosesList[15] = Hexose;
    neutralLosesList[16] = GlucuronicAcid;
    isRadical = std::vector<bool>(size, false);
    for (size_t i = 0; i < size; ++i) {
        const double rdbe = GetDBE(neutralLosesList[i]);
        double fraction = 0;
        std::modf(rdbe, &fraction);
        // Not radical
        if (fraction != 0) isRadical[i] = true;
    }
}

double DetectNeutralLoses::DetermineNeutralLoses(const std::vector<int> &elements) const {
    size_t detectedIndex = -1;
    for (size_t i = 0; i < neutralLosesList.size(); ++i) {
        if (neutralLosesList[i] == elements) {
            detectedIndex = i;
            break;
        }
    }
    if (detectedIndex == -1) return 0;
    if (isRadical[detectedIndex]) return scoreRadical;
    return score;
}

double DetectNeutralLoses::GetDBE(const std::vector<int>& elements) {
    // C H N O P S
    // 0th index is Carbon
    // 1st index is Hydrogen
    // 2nd Index is Nitrogen
    // 3rd Index is Oxygen
    // 4th Index is Phosphorus
    // 5th Index is Sulfur
    return elements[0] - (elements[1] * 0.5) + ((elements[2] + elements[4] ) * 0.5);
}
