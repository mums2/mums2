
#include <vector>
#include <string>
#include <Rcpp.h>
#include <RcppThread.h>
#include <mutex>

#include "CustomProgressBar/CliProgressBar.h"
#include "Decomposition/DecomposeMass.h"
#include "Decomposition/DecomposeMassInputData.h"
#include "FragmentationTree/FragmentationTree.h"
#include "FragmentationTree/GreedyHeuristic.h"


std::vector<DecompositionMassInputData> CreateInputData(const Rcpp::NumericVector& mzData, const Rcpp::List& masses) {
	const int size = static_cast<int>(masses.size());
	std::vector<DecompositionMassInputData> decompInput(size);
	for (int i = 0; i < size; i++) {
		const Rcpp::NumericVector& mass = masses[i];
		const int massSize = mass.size();
		const double parentMass = mzData[i];
		decompInput[i].masses.reserve(massSize);
		decompInput[i].parentMass = parentMass;
		for (int j = massSize - 1; j >= 0; j--) {
			if (parentMass < mass[j]) continue;
			decompInput[i].masses.emplace_back(mass[j]);
		}
	}
	return decompInput;
}

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppThread)]]
//[[Rcpp::export]]
std::vector<std::string> DeNovoMolecularFormulaPrediction(const Rcpp::NumericVector& mzData, const Rcpp::List& masses,
	const double ppm, const int numberOfThreads = 1) {
	const int size = static_cast<int>(masses.size());
	const std::vector<DecompositionMassInputData> inputData = CreateInputData(mzData, masses);
	std::mutex mutex;

	std::vector<std::vector<DecompResult>> allNodes(size);
	DecomposeMass decomposeMass;
	CliProgressBar progressBar;
	int counter = 0;
	Rcpp::Rcout << "Calculating decomposition masses..." << std::endl;
	RcppThread::parallelFor(0, size, [&inputData, &decomposeMass,
		&allNodes, &ppm, &mutex, &progressBar, &counter, &size](int i) {
			const std::vector<DecompResult> result =
			decomposeMass.GenerateMolecularFormulas(inputData[i], 1, ppm);
		mutex.lock();
		allNodes[i] = result;
		counter++;
		progressBar.update(static_cast<float>(counter)/static_cast<float>(size));

		mutex.unlock();
	}, numberOfThreads);
	progressBar.end_display();
	DetectNeutralLoses neutralLoses;
	CliProgressBar progressBar2;
	std::vector<std::string> resultantFormulas(size);
	counter = 0;
	Rcpp::Rcout << "Calculating fragmentation trees..." << std::endl;
	for (size_t i = 0; i < size; i++) {
		FragmentationTree tree(allNodes[i], inputData[i].parentMass);
		const std::vector<int>& colorRanges = tree.GetColorRanges();
		if (colorRanges.at(1) <= 0) continue; // Mean there are no decompositions!
		const int colorRangesSize = static_cast<int>(colorRanges.size());
		for (int j = colorRangesSize - 1 ; j >= 2; j--) {
			const int endingIndex = colorRanges.at(j);
			const int startingIndex =  colorRanges.at(j - 1);
			// Start One color Up
			const int fragmentEndPosition = colorRanges.at(j - 1);
			const int fragmentStartPosition = colorRanges.at(j - 2);
			RcppThread::parallelFor(fragmentStartPosition, fragmentEndPosition, [&tree, &neutralLoses, &startingIndex,
				&endingIndex](int j) {
				tree.AddMolecularFormulaToGraph(j, startingIndex, endingIndex, neutralLoses);
			}, numberOfThreads);
		}
		resultantFormulas[i] = GreedyHeuristic::CalculateHeuristic(tree);
		counter++;
		progressBar2.update(static_cast<float>(counter)/static_cast<float>(size));
	}
	progressBar2.end_display();
	return resultantFormulas;
}
