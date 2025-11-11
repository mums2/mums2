//
// Created by Gregory Johnson on 4/3/25.
//

#ifndef READSPECTRA_H
#define READSPECTRA_H
#include <string>
#include <Rcpp.h>
#include <fstream>
#include "../Spectra/Spectra.h"


class ReadSpectra {
public:
    static Rcpp::List ReadMGF(const std::string &filePath);
    static Rcpp::List ReadMSP(const std::string &filePath);
    static Spectra ReadSpectraFile(const std::string& filePath);

};



#endif //READSPECTRA_H
