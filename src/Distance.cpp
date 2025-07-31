#include "Distance/Distance.h"
#include <RcppThread.h>
#include <mutex>

#include "CustomProgressBar/CliProgressBar.h"

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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
void Distance::CalculateDistances(const double prec_threshold, const double cutoff,
    const ScoringFactory& scoreMethod, const int minPeaks, const int numberOfThreads) {
    const auto size = static_cast<int>(spectraList.size());
    CliProgressBar p;
    for(int i = 0; i < size; i++) {
        Spectra firstSpectra = spectraList[i];
        std::mutex mutex;
        RcppThread::parallelFor(i + 1, size, [this, &firstSpectra, &prec_threshold,
            &minPeaks, &scoreMethod, &i, &cutoff, &mutex](int j) {
            const Spectra secondSpectra = spectraList[j];
            bool hasScored = false;
            double score = -1;
              if (std::abs(firstSpectra.precursorMz - secondSpectra.precursorMz) < prec_threshold) {
                  score = 1 - scoreMethod.CalculateScore(firstSpectra, secondSpectra, minPeaks);
                  hasScored = true;
              }
              if (hasScored && score < cutoff) {
                  mutex.lock();
                  const SparseValue dist(i, j, score);
                  sparseMatrix.push(dist);
                  mutex.unlock();
              }

        }, numberOfThreads);
        p.update(static_cast<float>(i)/static_cast<float>(size));
    }
    p.end_display();
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
