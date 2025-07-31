#ifndef SPECTRA
#define SPECTRA

#include <string>
#include <vector>
struct Spectra final {
    Spectra(std::string name, const  std::vector<double>& mz, const  std::vector<double>& intensity, const double precursorMz):name(std::move(name)), 
    mz(mz),intensity(intensity),precursorMz(precursorMz) {};
    Spectra() = default;
    std::string name;
    std::vector<double> mz;
    std::vector<double> intensity;
    double precursorMz;
};

#endif //SPECTRA