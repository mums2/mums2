#include "Distance/Distance.h"


void Distance::CreateSpectraList(Rcpp::List data) {
    std::vector<double> pmz =  Rcpp::as<std::vector<double>>(data["pmz"]);
    std::vector<std::string> name =  Rcpp::as<std::vector<std::string>>(data["id"]);
    Rcpp::List mzInts = data["spectra"];
    
    const int n = pmz.size();
    for(int i = 0; i < n; i++){
        Rcpp::DataFrame specDataFrame = Rcpp::wrap(mzInts[i]);
        Spectra spec(name[i], specDataFrame["mz"], specDataFrame["intensity"], pmz[i]);
        spectraList.emplace_back(spec);
    }
}

void Distance::CalculateDistances(const double prec_threshold, double cutoff, const ScoringFactory& scoreMethod) {
    const auto size = static_cast<int>(spectraList.size());
    for(int i = 0; i < size; i++) {
        Spectra firstSpectra = spectraList[i];
        for(int j = i + 1; j < size; j++){
            Spectra secondSpectra = spectraList[j];
            
            if (std::abs(firstSpectra.precursorMz - secondSpectra.precursorMz) > prec_threshold) {
                continue;
            }
            const double score = scoreMethod.CalculateScore(firstSpectra, secondSpectra);

            if ((1 - score) > cutoff) {
                continue;
            }
            
            SparseValue dist(i, j, (1 - score));
            sparseMatrix.emplace_back(dist);
        }
    }
} 

Rcpp::DataFrame Distance::ExtractMatrix() {
    int len = sparseMatrix.size();
    Rcpp::IntegerVector i(len);
    Rcpp::IntegerVector j(len);
    Rcpp::NumericVector dist(len);
    
    for (int p = 0; p < len; p++) {
        i[p] = sparseMatrix[p].getI() + 1;
        j[p] = sparseMatrix[p].getJ() + 1;
        dist[p] = sparseMatrix[p].getDistance();
    }

    Rcpp::DataFrame m = Rcpp::DataFrame::create(Rcpp::Named("i") = i,
                                                Rcpp::Named("j") = j,
                                                Rcpp::Named("dist") = dist);

    return m;
    
}
