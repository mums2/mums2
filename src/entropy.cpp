//
// Created by gregj on 10/2/2024.
//
#include "ScoringMethods/SpectralEntropy/entropy.h"


Entropy::Entropy(const Rcpp::List& parameters):ms2_tolerance_in_da_param(parameters["ms2_tolerance_in_da"]),
    ms2_tolerance_in_ppm_param(parameters["ms2_tolerance_in_ppm"]),clean_spectra_param(parameters["clean_spectra"]),
    min_mz_param(parameters["min_mz"]), max_mz_param(parameters["max_mz"]), noise_threshold_param(parameters["noise_threshold"]),
    max_peak_num_param(parameters["max_peak_num"]), weighted_param(parameters["weighted"]){}

double Entropy::CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) {
    return CalculateEntropySimilarity(firstSpectra.mz, firstSpectra.intensity,
        secondSpectra.mz, secondSpectra.intensity);
}



std::vector<double> Entropy::convert_matrix_to_vector(const std::vector<double>& mz,
    const std::vector<double>& intensity) {
    const size_t size = mz.size() * 2;
    std::vector<double> mzIntensityVector(size);
    size_t count = 0;
    for (size_t i = 0; i < size; i+=2) {
        mzIntensityVector[i] = mz[count];
        mzIntensityVector[i + 1] = intensity[count++];
    }
    return mzIntensityVector;
}


double Entropy::CalculateEntropySimilarity(const std::vector<double>& listOneMz,
                                    const std::vector<double>& listOneInt,
                                    const std::vector<double>& listTwoMz,
                                    const std::vector<double>& listTwoInt) {


    std::vector<double> peaks_a_vec = convert_matrix_to_vector(listOneMz, listOneInt);
    const int peaks_a_len = peaks_a_vec.size() / 2;
    double* peaks_a_ptr = peaks_a_vec.data();
    std::vector<double> peaks_b_vec = convert_matrix_to_vector(listTwoMz, listTwoInt);
    const int peaks_b_len = peaks_b_vec.size() / 2;
    double* peaks_b_ptr = peaks_b_vec.data();
    if (weighted_param)
        return calculate_entropy_similarity(
        peaks_a_ptr, peaks_a_len,
        peaks_b_ptr, peaks_b_len,
        ms2_tolerance_in_da_param, ms2_tolerance_in_ppm_param,
        clean_spectra_param,
        min_mz_param, max_mz_param,
        noise_threshold_param,
        max_peak_num_param);

    return calculate_unweighted_entropy_similarity( peaks_a_ptr, peaks_a_len,
        peaks_b_ptr, peaks_b_len,
        ms2_tolerance_in_da_param, ms2_tolerance_in_ppm_param,
        clean_spectra_param,
        min_mz_param, max_mz_param,
        noise_threshold_param,
        max_peak_num_param);
}