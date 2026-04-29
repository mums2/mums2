//
// Created by gregj on 4/26/2026.
//

#include "Decomposition/DecomposeMass.h"

#include <algorithm>
#include <complex>
#include <cmath>
#include <numeric>




char DecomposeMass::GetParity(const ims::ComposedElement& molecule, const int charge) const {
	const bool massEven =  static_cast<int>(molecule.getMass()) % 2 == 0 ? true : false;
	const bool nitrogenEven = static_cast<int>(molecule.getElementAbundance("N")) % 2 == 0 ? true : false;
	const bool chargeEven = abs(charge) % 2 == 0 ? true : false;

	return (massEven ^ nitrogenEven ^ chargeEven) == true ? 'e' : 'o' ;
}


bool DecomposeMass::IsValidMyNitrogenRule(const ims::ComposedElement& molecule, const int z) const {
	const bool massOdd =  static_cast<int>(molecule.getNominalMass()) % 2 == 1;
	const bool massEven = !massOdd;

	const bool nitrogenOdd = static_cast<int>(molecule.getElementAbundance("N")) % 2 == 1;
	const bool nitrogenEven = !nitrogenOdd;

	const bool zOdd = z % 2 == 1;
	const bool zEven = !zOdd;

	return (  zEven & massEven & nitrogenEven )
	  |    (  zEven & massOdd  & nitrogenOdd  )
	  |    (  zOdd  & massEven & nitrogenOdd  )
	  |    (  zOdd  & massOdd  & nitrogenEven );

}

float DecomposeMass::GetDBE(const ims::ComposedElement& molecule, int z) const {
	// {{{

	return (1 + static_cast<int>(molecule.getElementAbundance("C"))
		+ static_cast<int>(molecule.getElementAbundance("Si"))
		- 0.5 * (static_cast<int>(molecule.getElementAbundance("H"))
			 +static_cast<int>(molecule.getElementAbundance("F"))
			 +static_cast<int>(molecule.getElementAbundance("Cl"))
			 +static_cast<int>(molecule.getElementAbundance("Br"))
			 +static_cast<int>(molecule.getElementAbundance("I")))
		+ 0.5 * (static_cast<int>(molecule.getElementAbundance("N"))
			 + static_cast<int>(molecule.getElementAbundance("P"))));
}


bool DecomposeMass::IsWithinElementRange(const ims::ComposedElement& molecule, const ims::ComposedElement& minElements,
	const ims::ComposedElement& maxElements) const {

	//iterating over Elements present in minElements
	for (auto it = minElements.getElements().begin();
		 it != minElements.getElements().end(); ++it) {

		const int minCount = static_cast<int>(minElements.getElementAbundance((it->first).getName()));

		const int count = static_cast<int>(molecule.getElementAbundance((it->first).getName()));

		if (count < minCount) {
			return false;
		}

	}

	//must iterate over Elements present in molecule and ignore cases of Elements not defined in maxElements
	for (auto it = maxElements.getElements().begin();
		 it != maxElements.getElements().end(); ++it) {

		const int maxCount = static_cast<int>(maxElements.getElementAbundance((it->first).getName()));

		const int count = static_cast<int>(molecule.getElementAbundance((it->first).getName()));

		// TODO: Fails e.g. for "C2N0"
		// JL: I fixed this error by changing "maxcount > 0" to "maxcount >= 0"; however this operator is only available for C++ >20 so it might not compile on old systems
		if (maxCount >= 0 && count > maxCount) {
			return false;
		}
	}

	return true;
}

void DecomposeMass::InitializeCHNOPS(ims::Alphabet& chnops, const int maxisotopes) const {

    ims::IsotopeDistribution::SIZE = maxisotopes;
	ims::IsotopeDistribution::ABUNDANCES_SUM_ERROR = 0.00001;

// Hydrogen

	ims::IsotopeDistribution::nominal_mass_type massH = 1;
	ims::IsotopeDistribution::peaks_container peaksH;
	peaksH.emplace_back(0.00782503, 0.99988500);
	peaksH.emplace_back(0.01410178, 0.00011500);
	peaksH.emplace_back(0.01604928, 0.00000000);

	ims::IsotopeDistribution distributionH(peaksH, massH);

// Oxygen
	ims::IsotopeDistribution::nominal_mass_type massO = 16;
	ims::IsotopeDistribution::peaks_container peaksO;
	peaksO.emplace_back(-0.00508538, 0.99757000);
	peaksO.emplace_back(-0.00086824, 0.00038000);
	peaksO.emplace_back(-0.00084039, 0.00205000);

	ims::IsotopeDistribution distributionO(peaksO, massO);

// Carbonate

	ims::IsotopeDistribution::nominal_mass_type massC = 12;
	ims::IsotopeDistribution::peaks_container peaksC;
	peaksC.emplace_back(0.0, 0.98930000);
	peaksC.emplace_back(0.00335484, 0.01070000);
	peaksC.emplace_back(0.00324199, 0.00000000);
	ims::IsotopeDistribution distributionC(peaksC, massC);

// Nitrogen

	ims::IsotopeDistribution::nominal_mass_type massN = 14;
	ims::IsotopeDistribution::peaks_container peaksN;
	peaksN.emplace_back(0.00307400, 0.99636000);
	peaksN.emplace_back(0.00010890, 0.00364000);

	ims::IsotopeDistribution distributionN(peaksN, massN);

// Sulfur

	ims::IsotopeDistribution::nominal_mass_type massS = 32;
	ims::IsotopeDistribution::peaks_container peaksS;
	peaksS.emplace_back(-0.02792883, 0.94990000);
	peaksS.emplace_back(-0.02854109, 0.00750000);
	peaksS.emplace_back(-0.03213300, 0.04250000);
	peaksS.emplace_back(-0.00000000, 0.00000000);
	peaksS.emplace_back(-0.03291929, 0.00010000);


	ims::IsotopeDistribution distributionS(peaksS, massS);

// Phosphor
	ims::IsotopeDistribution::nominal_mass_type massP = 31;
	ims::IsotopeDistribution::peaks_container peaksP;
	peaksP.emplace_back(-0.02623800, 1.0);

	ims::IsotopeDistribution distributionP(peaksP, massP);

	ims::Alphabet::element_type H("H", distributionH);
	ims::Alphabet::element_type C("C", distributionC);
	ims::Alphabet::element_type N("N", distributionN);
	ims::Alphabet::element_type O("O", distributionO);
	ims::Alphabet::element_type P("P", distributionP);
	ims::Alphabet::element_type S("S", distributionS);

	chnops.push_back(H);
	chnops.push_back(C);
	chnops.push_back(N);
	chnops.push_back(O);
	chnops.push_back(P);
	chnops.push_back(S);
}




DecompResult DecomposeMass::GenerateResults(const std::multimap<double, ims::ComposedElement,
                                                          std::greater<double>>& scores, const int z,
                                                          const int color) const {

	// Build result set to be returned as a list to R.
	std::vector<std::string> formula;
	std::vector<double> score;
	std::vector<double> exactmass;
	// std::vector<int> charge;
	//std::vector<FragmentationNode> nodes;
	 formula.reserve(scores.size());
	 score.reserve(scores.size());
	 exactmass.reserve(scores.size());
	// charge.reserve(scores.size());


	// Chemical rules
	// std::vector<std::string> parity(scores.size());
	// std::vector<bool> valid(scores.size());
	// std::vector<double> DBE(scores.size());

	// std::vector<std::string> colNames(2);
	// colNames[0] = "mass";
	// colNames[1] = "intensity";

	// outputs molecules & their scores.

	for (const auto& scoreResult : scores) {
		if (!IsValidMyNitrogenRule(scoreResult.second, z)) {
			continue;
		}
		score.emplace_back(scoreResult.first);
		formula.emplace_back(scoreResult.second.getSequence());
		exactmass.emplace_back(scoreResult.second.getMass());
		// charge.emplace_back(z);

		// Chemical rules
		// parity[i] = GetParity(scoreResult.second, z);
		//
		//
		// DBE[i] = GetDBE(scoreResult.second, z);
		// nodes.emplace_back(color, -1,
		// scoreResult.first, 0, MolecularFormula(scoreResult.second.getSequence(),
		// 	scoreResult.second.getMass()));
	}
	if (formula.empty() && !scores.empty()) {
		// this means all the formulas were invalid
		const auto it = scores.cbegin();
		score.emplace_back(it->first);
		formula.emplace_back(it->second.getSequence());
		exactmass.emplace_back(it->second.getMass());
		// nodes.emplace_back(color, -1,
		// it->first, 0, MolecularFormula(it->second.getSequence(),
		// 	it->second.getMass()));
	}
	return DecompResult{formula, score, exactmass, std::vector<int>(formula.size(), color)};
}





std::multimap<double, ims::ComposedElement,
	std::greater<double> > DecomposeMass::DecomposeMassFormulas(double mass, double intensity, double ppm) const {
	std::vector<double> abundances {intensity};
	std::vector<double> masses {mass};
	constexpr int maxIsotopes = 10;
	constexpr double precision = 1.0e-05;
	double error = ppm + 0.0001 / masses[0] * 1000000.0;
	error *= masses[0] * 1.0e-06;
	ims::Alphabet alphabet;
	std::vector<std::string> elements_order;
	InitializeCHNOPS(alphabet, maxIsotopes);
	elements_order.emplace_back("C");
	elements_order.emplace_back("H");
	elements_order.emplace_back("N");
	elements_order.emplace_back("O");
	elements_order.emplace_back("P");
	elements_order.emplace_back("S");

	ims::Weights weights(alphabet.getMasses(), precision);
	weights.divideByGCD();

	ims::RealMassDecomposer decomposer(weights);

	double abundances_sum = 0.0;
	for (const double& abundance : abundances) {
		abundances_sum += abundance;
	}
	for (double& abundance : abundances) {
		abundance /= abundances_sum;
	}

	std::vector<double> peaklist_masses;
	ims::DistributionProbabilityScorer::abundances_container peaklist_abundances;
	for (std::vector<double>::size_type mi = 0; mi < masses.size() && mi < abundances.size(); ++mi) {
		peaklist_masses.push_back(masses[mi]);
		peaklist_abundances.push_back(abundances[mi]);
	}

	ims::DistributionProbabilityScorer scorer(peaklist_masses, peaklist_abundances);

	std::vector<std::pair<ims::ComposedElement, double>> nonnormalized_scores;

	std::multimap<double, ims::ComposedElement,
	std::greater<double> > scores;
	std::vector<std::vector<unsigned>> decompositions =
		decomposer.getDecompositions(masses[0], error);

	double accumulated_score = 0.0;
	ims::ComposedElement minElements("C0", alphabet);
	ims::ComposedElement maxElements("C999999", alphabet);



		for (const auto & decomposition : decompositions) {

		// creates a candidate molecule out of elemental composition and a set of elements
		ims::ComposedElement candidate_molecule(decomposition, alphabet);

		// Check minimum/maximum element counts
		if (!IsWithinElementRange(candidate_molecule, minElements, maxElements)) {
			continue;
		}

		// updates molecules isotope distribution (since its not calculated upon creation:
		// it would be time consuming before applying chemical filter)
		candidate_molecule.updateIsotopeDistribution();
		// updates molecules sequence in a order of elements(atoms) one would like it
		// to appear
		candidate_molecule.updateSequence(&elements_order);

		// gets a theoretical isotope distribution of the candidate molecule
		ims::IsotopeDistribution candidate_molecule_distribution =
				candidate_molecule.getIsotopeDistribution();

		// gets a sequence of the candidate molecule
		std::string candidate_molecule_sequence =
				candidate_molecule.getSequence();

		// extracts masses and abundances from isotope distribution of the candidate molecule
		std::vector<double> candidate_masses = candidate_molecule_distribution.getMasses();
		std::vector<double> candidate_abundances = candidate_molecule_distribution.getAbundances();

		// normalizes candidate abundances if the size of the measured peaklist is less than
		// the size of theoretical isotope distribution. This is always the case since our
		// theoretical distributions are limited to by default 10 peaks and measured peaklists contain
		// less than 10 peaks

		long long size =
			static_cast<long long>(std::min(peaklist_abundances.size(), candidate_abundances.size()));

		if (size < static_cast<long long>(candidate_abundances.size())) {
			// normalizes the isotope distribution abundances with respect to the number of elements in peaklist
			double sum = std::accumulate(candidate_abundances.begin(),
							candidate_abundances.begin() + size,
							0.0);
			if (std::fabs(sum - 1) > ims::IsotopeDistribution::ABUNDANCES_SUM_ERROR) {
				double scale = 1.0/sum;
				// for (const auto& abund  : candidate_abundances) {
				// 	RcppThread::Rcout << "Candidate Abundances: " << abund << std::endl;
				// }

				std::transform(candidate_abundances.begin(),			// begin of source range
					candidate_abundances.begin() + size,		// end of source range
					candidate_abundances.begin(), 			// destination
					[scale](const double a) { return a * scale; }); // operation (*scale)
			}

		}

		// calculates a score
		double score = scorer.score(candidate_masses, candidate_abundances);

		// stores the sequence with non-normalized score

		nonnormalized_scores.emplace_back(candidate_molecule, score);

		// accumulates scores
		accumulated_score += score;

	}

	for (const auto& nonnormalized_score : nonnormalized_scores) {
		double normalized_score = nonnormalized_score.second;
		if (accumulated_score > 0.0) {
			normalized_score /= accumulated_score;
		}
		// stores the sequence with the score
		scores.emplace(normalized_score, nonnormalized_score.first);
	}
	return scores;
}

std::vector<DecompResult> DecomposeMass::GenerateMolecularFormulas(const DecompositionMassInputData& inputData,
	const double intensity, const double ppm) const {
	const int size = static_cast<int>(inputData.masses.size() + 1);

	std::vector<DecompResult> results(size);
				results[0] = GenerateResults(DecomposeMassFormulas(inputData.parentMass, intensity, ppm),
					0, 0);

	if (results.empty()) return {};

	for (int i = 1; i < size; ++i) {
		const double currentMass = inputData.masses[i - 1];
		if (currentMass > inputData.parentMass) continue;

		results[i] = GenerateResults(DecomposeMassFormulas(currentMass, intensity, ppm),
					0, i);
	}
	return results;
}

SEXP DecomposeMass::DecompToRObject(const DecompResult& decompResult) {
	return Rcpp::List::create(Rcpp::Named("formulas") = decompResult.formula,
	                   Rcpp::Named("scores") = decompResult.score,
	                   Rcpp::Named("exactmass") = decompResult.exactmass,
	                   Rcpp::Named("color") = decompResult.color);
}
