#ifndef QUERY
#define QUERY

#include <string>
#include <utility>
#include <vector>
#include "../../Spectra/Spectra.h"

class Query {
    public:
    Query(std::string variableId, std::string ms2Id, const double ms2Mz,
          const double ms2Rt, const Rcpp::String& formula, const std::vector<double> &mz,
          const std::vector<double> &intensity):
    formula(formula), variableId(std::move(variableId)), ms2Id(std::move(ms2Id)),
    ms2Rt(ms2Rt), spectra(variableId, mz, intensity, ms2Mz){};
    

    
    void AddToScores(const double score) {
        scores.emplace_back(score);
    }

    void AddToRefComps(const int refComp) {
        refComps.emplace_back(refComp);
    }

    // Getters

    std::vector<int> GetRefComps() const {return refComps;}
    std::vector<double> GetScores() const {return scores;}

    Spectra GetSpectra() const {
        return spectra;
    }
    
    std::string getVariableId() {
        return variableId;
    }
    std::string getMs2Id() {
        return ms2Id;
    }
    const Rcpp::String& GetFormula() const {return formula;}
    double getMs2Pmz() const {
        return spectra.precursorMz;
    }
    
    double getMs2Rt() const {
        return ms2Rt;
    }
    
    std::vector<double> getSpecMz() const {
        return spectra.mz;
    }
    
    std::vector<double> getSpecIntensity()const {
        return spectra.intensity;
    }
    
    private:
    Rcpp::String formula;
    std::string variableId;
    std::string ms2Id;
    double ms2Rt;
    Spectra spectra;
    std::vector<int> refComps;
    std::vector<double> scores;
};

#endif //QUERY