#include "Distance/Distance.h"


void Distance::CreateSpectraList(Rcpp::List data) {
    std::vector<double> pmz =  Rcpp::as<std::vector<double>>(data["pmz"]);
    std::vector<std::string> name =  Rcpp::as<std::vector<std::string>>(data["id"]);
    Rcpp::List mzInts = data["spectra"];
    const size_t n = pmz.size();
    spectraList = std::vector<Spectra>(n);
    for(size_t i = 0; i < n; i++){
        Rcpp::DataFrame specDataFrame = Rcpp::wrap(mzInts[i]);
        const Spectra spec(name[i], specDataFrame["mz"], specDataFrame["intensity"], pmz[i]);
        spectraList[i] = spec;
    }
}

void Distance::CalculateDistances(const double prec_threshold, const double cutoff,
    const ScoringFactory& scoreMethod, const int minPeaks) {
    const auto size = static_cast<int>(spectraList.size());
    for(int i = 0; i < size; i++) {
        Spectra firstSpectra = spectraList[i];
        for(int j = i + 1; j < size; j++){
            Spectra secondSpectra = spectraList[j];
            
            if (std::abs(firstSpectra.precursorMz - secondSpectra.precursorMz) > prec_threshold) {
                continue;
            }
            const double score = scoreMethod.CalculateScore(firstSpectra, secondSpectra, minPeaks);

            if ((1 - score) > cutoff) {
                continue;
            }
            
            SparseValue dist(i, j, (1 - score));
            sparseMatrix.push(dist);
        }
    }
} 

Rcpp::DataFrame Distance::ExtractMatrix() {
    const int len = static_cast<int>(sparseMatrix.size());
    Rcpp::IntegerVector iIndex(len);
    Rcpp::IntegerVector jIndex(len);
    Rcpp::NumericVector dist(len);
    
    for (int i = 0; i < len; i++) {
        SparseValue value = sparseMatrix.front();
        iIndex[i] = value.i + 1;
        jIndex[i] = value.j + 1;
        dist[i] = value.distance;
        sparseMatrix.pop();
    }

    Rcpp::DataFrame m = Rcpp::DataFrame::create(Rcpp::Named("i") = iIndex,
                                                Rcpp::Named("j") = jIndex,
                                                Rcpp::Named("dist") = dist);

    return m;
    
}
