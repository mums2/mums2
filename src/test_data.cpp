
#include <vector>
#include <string>
#include <Rcpp.h>
#include <RcppThread.h>
#include <mutex>
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



//[[Rcpp::export]]
std::vector<std::string> DecomposeMasses(const Rcpp::NumericVector& mzData, const Rcpp::List& masses,
	const double ppm, const int numberOfThreads = 1) {
	const int size = static_cast<int>(masses.size());
	const std::vector<DecompositionMassInputData> inputData = CreateInputData(mzData, masses);
	std::mutex mutex;
	// std::vector<std::vector<FragmentationNode>> decompResults;
	std::vector<std::string> resultantFormulas;
	DecomposeMass decomposeMass;
	RcppThread::ProgressBar bar(size, 1);
	RcppThread::parallelFor(0, size, [&inputData, &decomposeMass, &resultantFormulas,
		&ppm, &mutex, &bar](int i) {
		const std::vector<FragmentationNode> result =
			decomposeMass.GenerateMolecularFormulas(inputData[i], 1, ppm);
		FragmentationTree tree(result, inputData[i].parentMass);
		const std::vector<int>& availableIndexes = tree.GetColorZeroFormulas();
		for (const int& availableIndex : availableIndexes) {
			tree.AddMolecularFormulaToGraph(availableIndex);
		}
		std::string resultantFormula = GreedyHeuristic::CalculateHeuristic(tree);
		mutex.lock();
		resultantFormulas.emplace_back(resultantFormula);
		mutex.unlock();
		++bar;
	}, numberOfThreads);
	return resultantFormulas;
}

//[[Rcpp::export]]
std::string DecomposeMasses3(const Rcpp::NumericVector& mzData, const Rcpp::List& masses,
	const double ppm, const int numberOfThreads = 1) {
	const int size = static_cast<int>(masses.size());
	const std::vector<DecompositionMassInputData> inputData = CreateInputData(mzData, masses);
	std::mutex mutex;
	std::vector<std::vector<FragmentationNode>> decompResults;
	DecomposeMass decomposeMass;
	// RcppThread::ProgressBar bar(size, 1);
	// RcppThread::parallelFor(0, size, [&inputData, &decomposeMass, &decompResults,
	// 	&ppm, &mutex, &bar](int i) {
	// 	std::vector<FragmentationNode> result = decomposeMass.GenerateMolecularFormulas(inputData[i], 1, ppm);
	// 	mutex.lock();
	// 	decompResults.emplace_back(result);
	// 	mutex.unlock();
	// 	++bar;
	// }, numberOfThreads);
	const std::vector<FragmentationNode> result =
		decomposeMass.GenerateMolecularFormulas(inputData[5], 1, ppm);
	FragmentationTree tree(result, inputData[5].parentMass);
	const std::vector<int>& availableIndexes = tree.GetColorZeroFormulas();
	for (const int& availableIndex : availableIndexes) {
		tree.AddMolecularFormulaToGraph(availableIndex);
	}
	return GreedyHeuristic::CalculateHeuristic(tree);
}

