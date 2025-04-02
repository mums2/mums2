

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
    Rcpp::NumericVector convert_matrix_to_vector(const Rcpp::NumericMatrix peaks);
    Rcpp::NumericMatrix convert_vector_to_matrix(const Rcpp::NumericVector peaks, int nrow);
    double r_calculate_spectral_entropy(const Rcpp::NumericMatrix peaks);
    Rcpp::NumericMatrix r_clean_spectrum(const Rcpp::NumericMatrix peaks,
                                        float min_mz, float max_mz,
                                        float noise_threshold,
                                        float min_ms2_difference_in_da, float min_ms2_difference_in_ppm,
                                        int max_peak_num,
                                        bool normalize_intensity);
    double r_calculate_unweighted_entropy_similarity(const Rcpp::NumericMatrix peaks_a,
                                                  const Rcpp::NumericMatrix peaks_b,
                                                  float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
                                                  bool clean_spectra,
                                                  float min_mz, float max_mz,
                                                  float noise_threshold,
                                                  int max_peak_num);
    double r_calculate_entropy_similarity(const Rcpp::NumericMatrix peaks_a,
                                       const Rcpp::NumericMatrix peaks_b,
                                       float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
                                       bool clean_spectra,
                                       float min_mz, float max_mz,
                                       float noise_threshold,
                                       int max_peak_num);
    double CalculateEntropySimilarity(const std::vector<double>& listOneMz,
                                        const std::vector<double>& listOneInt,
                                        const std::vector<double>& listTwoMz,
                                        const std::vector<double>& listTwoInt);
    Entropy() = default;
    explicit Entropy(const Rcpp::List&);

    double CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) override;

private:
    float ms2_tolerance_in_da_param{};
    float ms2_tolerance_in_ppm_param{};
    bool clean_spectra_param{};
    float min_mz_param{};
    float max_mz_param{};
    float noise_threshold_param{};
    int max_peak_num_param{};
    bool weighted_param{};
};

#endif // ENTROPY
