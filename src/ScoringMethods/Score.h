//
// Created by gregj on 10/2/2024.
//

#ifndef SCORE_H
#define SCORE_H
#include <Rcpp.h>
#include <vector>

#include "../Spectra/Spectra.h"

class Score {
public:
    virtual ~Score() = default;
    virtual double CalculateScore(const Spectra& firstSpectra, const Spectra& secondSpectra) = 0;

};

#endif //SCORE_H
