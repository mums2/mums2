//
// Created by gregj on 4/26/2026.
//

#ifndef DECOMPRESULT_H
#define DECOMPRESULT_H
#include <vector>
#include <string>
struct DecompResult {
    std::vector<std::string> formula;
    std::vector<double> score;
    std::vector<double> exactmass;
    std::vector<double> rde;
};
#endif //DECOMPRESULT_H
