//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <cstdint>
#include <algorithm>
#include <RcppThread.h>
#include <mutex>
#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "CustomProgressBar/CliProgressBar.h"
#include "CustomProgressBar/ETAProgressBar.h"
#include <regex>
#include "DataStructures/CommunityMatrix.h"
#include "DataStructures/CppMatrix.h"
#include "DiversityMetrics/Diversity.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityMetricFactory.h"
#include "FragmentationTree/FragmentationTree.h"
#include "FragmentationTree/GreedyHeuristic.h"
#include "HMDB/HumanMetabolomicsDB.h"
#include "Math/ParallelRandomNumberSitmo.h"
#include "Math/VectorMath.h"
#include "ScoringMethods/ScoringFactory.h"
#include "ScoringMethods/GNPS/GNPSScoringDynamicPriorityQueue.h"
#include "ScoringMethods/SpectralEntropy/entropy.h"
#include "Spectra/ReadSpectra.h"

CppMatrix CalculateDiversity(const CppMatrix& abundances, const std::string& diversityIndex) {
    std::string index = diversityIndex;
    std::transform(index.begin(), index.end(), index.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    CppMatrix results = diversity->CalculateDiversity(abundances, index);
    delete diversity;
    return results;
}

// [[Rcpp::export]]
SEXP CreateCommunityMatrix(Rcpp::NumericMatrix communityMatrix) {
    auto* matrix = new CommunityMatrix(communityMatrix);
    matrix->InitializeMatrix();
    return Rcpp::XPtr<CommunityMatrix>(matrix);
}

// [[Rcpp::export]]
SEXP GetCommunityMatrix(SEXP communityMatrix) {
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    return matrix.get()->GetCommunityMatrix();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const SEXP& communityMatrix, const uint64_t size,
    const uint64_t threshold, const int numberOfThreads, const int seed = 123) {

    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();
    const Rcpp::CharacterVector& rowNames = matrix.get()->GetRowNames();
    const Rcpp::CharacterVector& columnNames = matrix.get()->GetColumnNames();
    std::vector<std::string> names = Rcpp::as<std::vector<std::string> >(rowNames);
    const std::vector<std::vector<uint64_t>>& abundanceRanges = matrix.get()->GetAbundanceRanges();
    const std::vector<std::vector<uint64_t>>& communityAbundances = matrix.get()->GetCommunityAbundances();
    const std::vector<std::vector<uint64_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    const std::vector<uint64_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    std::vector<ParallelRandomNumberSitmo> rngEngines(row);
    for (int i = 0; i < row; ++i) {
        rngEngines[i] = ParallelRandomNumberSitmo(seed + i);
    }

    std::mutex mutex;
    RcppThread::parallelFor(0, row, [&communityAbundances, &eligibleIndexes,
    &abundanceRanges, &resultantMatrix, &sums, &rngEngines, &mutex, &size, &threshold](int i) {
        const std::vector<uint64_t> result = Rarefaction::Rarefy(communityAbundances[i], eligibleIndexes[i],
          abundanceRanges[i],rngEngines[i], size, sums[i], threshold);
        mutex.lock();
        for(const auto& index : eligibleIndexes[i]) {
          resultantMatrix(i, index) = result.at(index);
        }
        mutex.unlock();
    }, numberOfThreads);
    Rcpp::rownames(resultantMatrix) = rowNames;
    Rcpp::colnames(resultantMatrix) = columnNames;
    return resultantMatrix;
}

CppMatrix RarefactionCalculationParallelized(const std::vector<std::vector<uint64_t>>& communityAbundances,
    const std::vector<std::vector<uint64_t>>& eligibleIndexes,
    const std::vector<std::vector<uint64_t>>& abundanceRanges,
    const std::vector<uint64_t>& sums,
    ParallelRandomNumberSitmo& rngEngine,
    const int rows,
    const int columns,
    const uint64_t size,
    const uint64_t threshold) {

    std::vector<double> resultantMatrix(rows * columns);
    for (int i = 0; i < rows; i++) {
        const std::vector<uint64_t> data = Rarefaction::Rarefy(communityAbundances[i], eligibleIndexes[i],
        abundanceRanges[i],rngEngine, size, sums[i], threshold);
        size_t counter = 0;
        for (int j = columns * i; j < columns * i + columns; j++) {
            resultantMatrix[j] = static_cast<double>(data[counter++]);
        }
    }
    return CppMatrix(resultantMatrix, rows, columns);
}
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const SEXP& communityMatrix, const std::string& index,
    const uint64_t size, const uint64_t threshold, const bool subsample, const int numberOfThreads,
    const int iterations = 1000, const int seed = 123) {
    CliProgressBar p;
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);
    const Rcpp::CharacterVector samples = communityObject.get()->GetSampleNames();
    const size_t matrixRowSize = samples.size();
    int row = samples.size();
    if (index == "simpson" || index == "shannon")
        row = 1;
    CppMatrix diversityMatrix(std::vector<double>(row * samples.size(), 0), row, samples.size());
    const std::vector<std::vector<uint64_t>>& abundanceRanges = communityObject.get()->GetAbundanceRanges();
    const std::vector<std::vector<uint64_t>>& communityAbundances = communityObject.get()->GetCommunityAbundances();
    const std::vector<std::vector<uint64_t>>& eligibleIndexes = communityObject.get()->GetColumnEligibleIndexes();
    const std::vector<uint64_t>& sums = communityObject.get()->GetSums();
    std::vector<ParallelRandomNumberSitmo> rngEngines(iterations);
    for (int i = 0; i < iterations; ++i) {
        rngEngines[i] = ParallelRandomNumberSitmo(seed + i);
    }
    std::mutex mutex;
    int columnSize = communityObject.get()->GetColumn();
    int currentProgress = 0;
    const CppMatrix& matrix = communityObject.get()->GetCppMatrixOfAbundances();
    RcppThread::parallelFor(0, iterations, [&communityAbundances, &eligibleIndexes, &abundanceRanges,
        &diversityMatrix, &rngEngines, &matrix, &sums, &matrixRowSize, &columnSize, &size, &threshold,
        &subsample, &index, &iterations, &mutex, &p, &currentProgress](int i) {
        CppMatrix rarefyMatrix;

        if (subsample) {
            rarefyMatrix = RarefactionCalculationParallelized(communityAbundances, eligibleIndexes,
                 abundanceRanges, sums, rngEngines[i], matrixRowSize, columnSize, size, threshold);
        }
        mutex.lock();
        if (subsample)
            diversityMatrix += CalculateDiversity(rarefyMatrix, index);
        else
            diversityMatrix += CalculateDiversity(matrix, index);
        p.update(static_cast<float>(currentProgress++)/static_cast<float>(iterations));
        mutex.unlock();
    }, numberOfThreads);
    diversityMatrix/=iterations;
    p.end_display();

    Rcpp::NumericMatrix resultantMatrix = diversityMatrix.ToRcppMatrix();
    Rcpp::colnames(resultantMatrix) = samples;
    if(diversityMatrix.GetRowSize() <= 1) return resultantMatrix; // alpha diversity
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::List ReadMgf(const std::string& path) {
    return ReadSpectra::ReadMGF(path);
}

// [[Rcpp::export]]
Rcpp::List ReadMsp(const std::string& path) {
    return ReadSpectra::ReadMSP(path);
}

// [[Rcpp::export]]
Rcpp::NumericVector CompareMS2Ms1(const Rcpp::NumericVector& mz2, const Rcpp::NumericVector& mz1,
    const Rcpp::NumericVector& rt2, const Rcpp::NumericVector& rt1, const double mzThreshold,
    const double rtThreshold) {
    const auto currentSize = static_cast<size_t>(mz1.size());
    Rcpp::NumericVector resultsIndexes(currentSize, -1); // -1 means no match
    for (size_t i = 0; i < currentSize; i++) {
        double currentMz1 = mz1[i];
        double currentRt1 = rt1[i];
        Rcpp::NumericVector mzError = Rcpp::abs(currentMz1 - mz2) * 1e6 / currentMz1;
        Rcpp::NumericVector rtError = Rcpp::abs(currentRt1 - rt2);
        double bestDotProduct = 0;
        for (int j = 0; j < mzError.size(); j++) { // Pick score with the closest dotProduct value
            if (mzError[j] > mzThreshold || rtError[j] > rtThreshold) continue; // Over the threshold
            
            // Otherwise
            // Check if the similarity score (the dot product) is closer than the last one
            // If so replace
            double dotProduct = VectorMath::CosineScore({mz1[i], rt1[i]}, {mz2[j], rt2[j]});
            if (dotProduct < bestDotProduct) continue;
            resultsIndexes[i] = j + 1; // To match with R indexes add 1
            bestDotProduct = dotProduct;
            if (bestDotProduct >= 1) break; // If they are equal, we should break the loop and move on
        }
        
    }
    return resultsIndexes;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::export]]
std::string ComputeFragmentationTree(const Rcpp::List& molecularFormulas,
    const double parentMass, const int numberOfThreads) {
    const int size = molecularFormulas.size();
    FragmentationTree tree(molecularFormulas, parentMass);
    RcppThread::parallelFor(0, size, [&tree](int i) {
        tree.AddMolecularFormulaToGraph(i);
    }, numberOfThreads);
    return GreedyHeuristic::CalculateHeuristic(tree);
}


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

// [[Rcpp::export]]
std::vector<double> TestSimilarity(const std::vector<double>& mzOne,  std::vector<double>& intOne,
    const std::vector<double>& mzTwo,  std::vector<double>& intTwo, Rcpp::List scoringParamsSpectral,
    Rcpp::List ScoringParamsCosine, const double shift) {
    const ScoringFactory factory(ScoringParamsCosine);
    GNPSScoringDynamicPriorityQueue gnps;
    Entropy entropy;
    double res = entropy.CalculateEntropySimilarity(mzOne, intOne, mzTwo, intTwo);
    const auto res2 = gnps.ScoreRData(mzOne, intOne, mzTwo, intTwo, 0.5, shift);
    Rcpp::Rcout << "SpectralEntropy: " << res << std::endl;
    Rcpp::Rcout << "GNPS: " << res2[0] << std::endl;
    return(gnps.ScoreRData(mzOne, intOne, mzTwo, intTwo, 0.5, shift));
}

// [[Rcpp::export]]
void DistanceToPhylipFile(const Rcpp::NumericMatrix& mat, const std::vector<std::string>& names, const std::string& fileName) {
    const int size = mat.rows();
    std::ofstream output(fileName);
    // output.open(fileName);
    if (!output.is_open()) {
        Rcpp::stop("file not open");
    }
    output << size << std::endl;
    const size_t  nameSize = names.size();
    for (size_t i = 0; i < nameSize; i++) {
        output << names[i] << "\t";
        for (size_t j = 0; j < nameSize; j++) {
            output << mat(i, j) << "\t";
        }
        output << std::endl;
    }
    output.close();
}

// [[Rcpp::export]]
Rcpp::List DistanceDataFrameToMatrix(const Rcpp::DataFrame& distanceDataFrame) {
    const int size = distanceDataFrame.rows();
    const std::vector<int>& iIndexes = distanceDataFrame["i"];
    const std::vector<int>& jIndexes = distanceDataFrame["j"];
    const std::vector<double>& distances = distanceDataFrame["dist"];
    Rcpp::NumericMatrix result(size, size);
    std::unordered_map<int, int> nameSwap;
    int currentIndex = 0;

    for (int i = 0; i < size; i++) {
        if (nameSwap.find(iIndexes[i]) == nameSwap.end()) {
            nameSwap[iIndexes[i]] = currentIndex;
            currentIndex++;
        }
        if (nameSwap.find(jIndexes[i]) == nameSwap.end()) {
            nameSwap[jIndexes[i]] = currentIndex;
            currentIndex++;
        }
        result(nameSwap[iIndexes[i]], nameSwap[jIndexes[i]]) = distances[i];
        result(nameSwap[jIndexes[i]], nameSwap[iIndexes[i]]) = distances[i];
    }
    const int nameSize = nameSwap.size();
    Rcpp::IntegerVector iNames(nameSize);
    Rcpp::IntegerVector jNames(nameSize);
    int count = 0;
    for (const auto& keyPairs: nameSwap) {
        iNames[count] = keyPairs.first;
        jNames[count++] = keyPairs.second;
    }
    return Rcpp::List::create(Rcpp::Named("i") = iNames,
                              Rcpp::Named("j") = jNames,
                              Rcpp::Named("result") = result);
    // return result;
}
// [[Rcpp::export]]
void ToColumnFile(const Rcpp::DataFrame& matrix, const std::vector<std::string>& names,
    const std::string& fileName) {
    std::ofstream output(fileName);
    if (!output.is_open()) {
        Rcpp::stop("File not open");
    }
    const std::vector<int>& currentIndexes = matrix["i"];
    const std::vector<int>& otherIndexes = matrix["j"];
    const std::vector<double>& distances = matrix["dist"];

   for (size_t i = 0; i < currentIndexes.size(); i++) {
       const std::string& name = names[currentIndexes[i] - 1];
       const std::string& otherNames = names[otherIndexes[i] - 1];
       output << name << "\t" << otherNames << "\t" << distances[i] << std::endl;
   }
}

void ToPhylip(const Rcpp::DataFrame& matrix, const std::vector<std::string>& names,
    const std::string& fileName) {
    std::ofstream output(fileName);
    if (!output.is_open()) {
        Rcpp::stop("File not open");
    }
    output << names.size() << std::endl;
    const std::vector<int>& currentIndexes = matrix["i"];
    const std::vector<int>& otherIndexes = matrix["j"];
    const std::vector<double>& distances = matrix["dist"];

    for (size_t i = 0; i < currentIndexes.size(); i++) {
        const std::string& name = names[currentIndexes[i] - 1];
        const std::string& otherNames = names[otherIndexes[i] - 1];
        output << name << "\t" << otherNames << "\t" << distances[i] << std::endl;
    }
}

// [[Rcpp::export]]
void AddSpectrumToHMDBData(Rcpp::List& hmdbData, const std::vector<std::string>& metaboliteNames,
    const std::vector<std::string>& fileNames) {
    std::unordered_map<std::string, std::vector<int>> hmdbDataMap;
    Rcpp::List ls = hmdbData[0];
    Rcpp::List info = ls["info"];
    std::vector<std::string> keys = Rcpp::as<std::vector<std::string>>(info["key"]);
    for (int i = 0; i < hmdbData.size(); i++) {
        Rcpp::List ls = hmdbData[i];
        Rcpp::List info = ls["info"];
        const std::string name = Rcpp::as<std::string>(info["spectra_name"]);
        hmdbDataMap[name] = {};
    }

    for (size_t i = 0; i < metaboliteNames.size(); i++) {
        if (hmdbDataMap.find(metaboliteNames[i]) == hmdbDataMap.end()) continue;
        hmdbDataMap[metaboliteNames[i]].emplace_back(i);
    }
    for (int i = 0; i < hmdbData.size(); i++) {
        Rcpp::List ls = hmdbData[i];
        Rcpp::List info = ls["info"];
        Rcpp::List spectra = ls["spec"];
        Rcpp::List spec = hmdbData[i];
        const std::string name = Rcpp::as<std::string>(info["spectra_name"]);
        for (const auto& index : hmdbDataMap[name]) {
            std::ifstream input(fileNames[index]);
            if (!input.is_open()) {
                Rcpp::warning("File not found: " + fileNames[index]);
                continue;
            }
            double mz, intensity;
            std::list<double> mzValues;
            std::list<double> intensityValues;
            while (input >> mz >> intensity) {
                mzValues.push_back(mz);
                intensityValues.push_back(intensity);
            }
            spectra["mz"] = wrap(mzValues);
            spectra["inÍtensity"] = wrap(intensityValues);
        }
    }

}

// [[Rcpp::export]]
SEXP CreateHumanMetabolomicsDB() {
    auto* hmdb = new HumanMetabolomicsDB();
    return Rcpp::XPtr<HumanMetabolomicsDB>(hmdb);
}

// [[Rcpp::export]]
void AddHumanMetabolomicNode(SEXP& hmdbPtr, const std::vector<std::string>& names,
    const std::vector<std::string>& values) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->AddHumanMetabolomicNode(HumanMetabolomicsDBNode(names,values));
}

// [[Rcpp::export]]
void PrintHMDBNames(const SEXP& hmdbPtr) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->PrintHumanMetabolomicsDB();
}
// [[Rcpp::export]]
void AddSpectra(SEXP& hmdbPtr, const std::vector<std::string>& fileNames,
    const std::vector<std::string>& databaseNames) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    for (size_t i = 0; i < fileNames.size(); i++) {
       hmdbPointer.get()->AddSpectraFiles(fileNames[i], databaseNames[i]);
    }
}

// [[Rcpp::export]]
void ProcessMs2Files(SEXP& hmdbPtr) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->ProcessSpectraFiles();
}

// [[Rcpp::export]]
Rcpp::List GetList(SEXP& hmdbPtr) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    return hmdbPointer.get()->ConstructDataBase();
}

// [[Rcpp::export]]
void ReadSpectraFile(const std::string& filePath) {
    ReadSpectra::ReadSpectraFile(filePath);
}