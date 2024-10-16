#include <Rcpp.h>
#include "Spectra/Spectra.h"
#include "ScoringMethods/ScoringFactory.h"
#include <vector>
#include <string.h>

// [[Rcpp::export]]
double ScoreMs2(std::string name1, std::vector<double> mz1, std::vector<double> intensity1, double precursorMz1, std::string name2, std::vector<double> mz2, std::vector<double> intensity2, double precursorMz2, Rcpp::List parameters) {
    ScoringFactory factory(parameters);
    Spectra firstSpectra(name1, mz1, intensity1, precursorMz1);
    Spectra secondSpectra(name2, mz2, intensity2, precursorMz2);
    
    double score = factory.CalculateScore(firstSpectra, secondSpectra);

    return score;
}