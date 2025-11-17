//
// Created by Gregory Johnson on 4/3/25.
//

#ifndef READSPECTRA_H
#define READSPECTRA_H
#include <string>
#include <Rcpp.h>
#include <fstream>

#include "../AnnotationStructure/AnnotationNode.h"


class ReadSpectra {
public:
    Rcpp::List ReadMGF(const std::string &filePath);
    Rcpp::List ReadMSP(const std::string &filePath);
    std::vector<AnnotationNodeData> ReadMSPSpectra(const std::string &filePath);
};



#endif //READSPECTRA_H
