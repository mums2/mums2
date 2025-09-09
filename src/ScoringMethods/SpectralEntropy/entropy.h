

#ifndef ENTROPY
#define ENTROPY
#define SPEC_TYPE
#include "../Score.h"
typedef double float_spec;

#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

extern "C"
{
    #include "CleanSpectrum.h"
    #include "SpectralEntropy.h"
}

class Entropy final : public Score {
public:
    std::vector<double> convert_matrix_to_vector(const std::vector<double>& mz,
                                                const std::vector<double>& intensity);
    double CalculateEntropySimilarity(const std::vector<double>& listOneMz,
                                        const std::vector<double>& listOneInt,
                                        const std::vector<double>& listTwoMz,
                                        const std::vector<double>& listTwoInt);
    Entropy() = default;
    explicit Entropy(const Rcpp::List&);

    double CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) override;

private:
    double ms2_tolerance_in_da_param{};
    double ms2_tolerance_in_ppm_param{};
    bool clean_spectra_param{};
    double min_mz_param{};
    double max_mz_param{};
    double noise_threshold_param{};
    int max_peak_num_param{};
    bool weighted_param{};
};

#endif // ENTROPY
