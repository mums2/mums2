//
// Created by Gregory Johnson on 11/19/25.
//

#include <Rcpp.h>
#include "CustomProgressBar/CliProgressBar.h"

// [[Rcpp::export]]
SEXP CreateProgressBarObject() {
    auto* progressBar = new CliProgressBar();
    return Rcpp::XPtr<CliProgressBar>(progressBar);
}

// [[Rcpp::export]]
void IncrementProgressBar(SEXP& progressBar, const float progress) {
    const Rcpp::XPtr<CliProgressBar> cliProgressBar(progressBar);
    cliProgressBar.get()->update(progress);
}

// [[Rcpp::export]]
void DestroyProgressBar(SEXP& progressBar) {
    const Rcpp::XPtr<CliProgressBar> cliProgressBar(progressBar);
    cliProgressBar.get()->end_display();
}

