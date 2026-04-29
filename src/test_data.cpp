
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
		decompInput[i].masses = std::vector<double>(massSize, 0);
		decompInput[i].parentMass = mzData[i];
		for (int j = 0; j < massSize; j++) {
			decompInput[i].masses[j] = mass[j];
		}
	}
	return decompInput;
}




// std::vector<std::string> DecomposeMasses(const Rcpp::NumericVector& mzData, const Rcpp::List& masses,
// 	const double ppm, const int numberOfThreads = 1) {
// 	const int size = static_cast<int>(masses.size());
// 	const std::vector<DecompositionMassInputData> inputData = CreateInputData(mzData, masses);
// 	std::mutex mutex;
// 	std::vector<std::string> resultantFormulas;
// 	DecomposeMass decomposeMass;
// 	CliProgressBar progressBar;
// 	float counter = 0;
// 	RcppThread::parallelFor(0, size, [&inputData, &decomposeMass, &resultantFormulas,
// 		&ppm, &mutex, &progressBar, &counter, &size](int i) {
// 		const std::vector<FragmentationNode> result =
// 			decomposeMass.GenerateMolecularFormulas(inputData[i], 1, ppm);
// 		FragmentationTree tree(result, inputData[i].parentMass);
// 		const std::vector<int>& availableIndexes = tree.GetColorZeroFormulas();
// 		for (const int& availableIndex : availableIndexes) {
// 			tree.AddMolecularFormulaToGraph(availableIndex);
// 		}
// 		std::string resultantFormula = GreedyHeuristic::CalculateHeuristic(tree);
// 		mutex.lock();
// 		resultantFormulas.emplace_back(resultantFormula);
// 		progressBar.update(++counter/static_cast<float>(size));
// 		mutex.unlock();
// 	}, numberOfThreads);
// 	return resultantFormulas;
// }
//

//[[Rcpp::export]]
std::vector<std::string> DecomposeMassesOther(const Rcpp::NumericVector& mzData, const Rcpp::List& masses,
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

	CliProgressBar progressBar2;
	std::vector<std::string> resultantFormulas(size);
	counter = 0;
	Rcpp::Rcout << "Calculating fragmentation trees..." << std::endl;
	for (size_t i = 0; i < size; i++) {
		FragmentationTree tree(allNodes[i], inputData[i].parentMass);
		const std::vector<int>& availableIndexes = tree.GetColorZeroFormulas();
		const int currentSize = static_cast<int>(availableIndexes.size());
		if (currentSize <= 0) continue; // Mean there are no decompositions!
		RcppThread::parallelFor(0, currentSize, [&tree, &availableIndexes](int j) {
			tree.AddMolecularFormulaToGraph(availableIndexes[j]);
		}, numberOfThreads);
		resultantFormulas[i] = GreedyHeuristic::CalculateHeuristic(tree);
		counter++;
		progressBar2.update(static_cast<float>(counter)/static_cast<float>(size));
	}
	progressBar2.end_display();
	return resultantFormulas;
}




// void DecomposeMasses1(const double mass,
// 	const double ppm, const int numberOfThreads = 1) {
// 	constexpr DecomposeMass decomposeMass;
// 		const auto result =
// 			decomposeMass.DecomposeMassFormulas(mass, 1, ppm);
// 	const auto score = decomposeMass.GenerateResults(result, 0, 0);
// 	FragmentationTree tree(allNodes[i], inputData[i].parentMass);
// 	const std::vector<int>& availableIndexes = tree.GetColorZeroFormulas();
// 		tree.AddMolecularFormulaToGraph(availableIndexes[j]);
//
// }

//[[Rcpp::export]]
std::string DecomposeMasses3(const Rcpp::NumericVector& mzData, const std::vector<double>& masses,
	const double ppm, const int numberOfThreads = 1) {
	const int size = static_cast<int>(masses.size());
	// const std::vector<DecompositionMassInputData> inputData = CreateInputData(mzData, masses);
	std::mutex mutex;
	std::vector<std::vector<FragmentationNode>> decompResults;
	DecomposeMass decomposeMass;
	CliProgressBar p;
	RcppThread::ProgressBar bar(size, 5);
	int counter = 0;
	RcppThread::parallelFor(0, size, [&masses, &decomposeMass, &decompResults,
		&ppm, &p, &size, &counter, &mutex](int i) {
		auto result = decomposeMass.DecomposeMassFormulas(masses[i], 1, ppm);

		// decompResults.emplace_back(result);
		// mutex.unlock();
		mutex.lock();
		counter++;
		p.update(static_cast<float>(counter)/static_cast<float>(size));
		mutex.unlock();
	}, numberOfThreads);
	p.end_display();
	return "";
}
