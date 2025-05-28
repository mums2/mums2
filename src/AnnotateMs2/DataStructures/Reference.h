#ifndef REFERENCE
#define REFERENCE

#include <string>
#include <vector>
#include "../../Spectra/Spectra.h"

class Reference final {
    public:
    Reference(const int index, const double precursorMz, const Rcpp::String& formula, const std::vector<double> &mz,
          const std::vector<double>& intensity):
    index(index),formula(formula), spectra("", mz, intensity, precursorMz) {};

    int getIndex() const {
        return index;
    }
    
    double getPrecursorMz() const {
        return spectra.precursorMz;
    }

    const Rcpp::String& GetFormula() const {return formula;}
    std::vector<double> getSpecMz() const {
        return spectra.mz;
    }
    
    std::vector<double> getSpecIntensity() const {
        return spectra.intensity;
    }
    Spectra GetSpectra() const {
        return spectra;
    }
    
    private:
    int index;
    // double precursorMz;
    Rcpp::String formula;
    std::vector<double> mz;
    std::vector<double> intensity;
    Spectra spectra;
};

#endif //REFERENCE