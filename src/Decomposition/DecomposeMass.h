//
// Created by gregj on 4/26/2026.
//

#ifndef DECOMPOSEMASS_H
#define DECOMPOSEMASS_H
#include "../RdisopHeaderFiles/alphabet.h"
#include "../RdisopHeaderFiles/distributionprobabilityscorer.h"
#include "../RdisopHeaderFiles/realmassdecomposer.h"
#include "../RdisopHeaderFiles/weights.h"
#include "../RdisopHeaderFiles/composedelement.h"
#include "DecompResult.h"
#include "DecomposeMassInputData.h"
#include "../FragmentationTree/FragmentationNode.h"
#include <Rcpp.h>



class DecomposeMass {
public:
    std::multimap<double, ims::ComposedElement,
	std::greater<double>> DecomposeMassFormulas(double mass, double intensity, double ppm = 2) const;
    std::vector<std::multimap<double, ims::ComposedElement,
    std::greater<double> >> GenerateMolecularFormulas(const DecompositionMassInputData& inputData,
        double intensity = 1, double ppm = 2) const;
    SEXP DecompToRObject(const DecompResult& decompResult);
private:
    void InitializeCHNOPS(ims::Alphabet& chnops, int maxisotopes) const;
    char GetParity(const ims::ComposedElement& molecule, int charge=0) const;
    bool IsValidMyNitrogenRule(const ims::ComposedElement& molecule, int z) const;
    float GetDBE(const ims::ComposedElement& molecule, int z) const;
    bool IsWithinElementRange(const ims::ComposedElement& molecule, const ims::ComposedElement& minElements,
    const ims::ComposedElement& maxElements) const;
public:
    std::vector<FragmentationNode> GenerateResults(const std::multimap<double, ims::ComposedElement,
                                                          std::greater<double>>& scores, int z, int color) const;


};



#endif //DECOMPOSEMASS_H
