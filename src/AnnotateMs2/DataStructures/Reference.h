#ifndef REFERENCE
#define REFERENCE

#include <string>
#include <vector>
#include "../../Spectra/Spectra.h"

class Reference final {
    public:
    Reference(int index, double precursorMz, std::vector<double> mz,
          std::vector<double> intensity):
    index(index), spectra("", mz, intensity, precursorMz) {};
    
    int getIndex() {
        return index;
    }
    
    double getPrecursorMz() {
        return spectra.precursorMz;
    }
    
    std::vector<double> getSpecMz() {
        return spectra.mz;
    }
    
    std::vector<double> getSpecIntensity() {
        return spectra.intensity;
    }
    Spectra GetSpectra() const {
        return spectra;
    }
    
    private:
    int index;
    double precursorMz;
    std::vector<double> mz;
    std::vector<double> intensity;
    Spectra spectra;
};

#endif //REFERENCE