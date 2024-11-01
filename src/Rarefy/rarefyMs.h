#ifndef RAREFYMS_H
#define RAREFYMS_H

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

class RarefyMs {
    public:
    DataFrame rarefyMs(IntegerVector feature, IntegerVector abund, int size, int threshold);
};


#endif  // RAREFYMS_H


