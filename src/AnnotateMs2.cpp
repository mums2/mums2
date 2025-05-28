#include "AnnotateMs2/AnnotateMs2.h"
#include "Utils/Utils.h"
#include <Rcpp.h>

#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"


void AnnotateMs2::createQueryList(const std::vector<std::string>& variableId, const std::vector<std::string> &ms2Id,
                                  const std::vector<float> &ms2Mz, const std::vector<float> &ms2Rt, const Rcpp::StringVector& formulas,
                                  const Rcpp::List& ms2Spectra) {
              
              const size_t n = ms2Id.size();
              for (size_t i = 0; i < n; i++) {
                  Rcpp::DataFrame specDf = Rcpp::wrap(ms2Spectra[i]);

                  Query query(variableId[i], ms2Id[i], ms2Mz[i], ms2Rt[i],
                      formulas[i], specDf["mz"], specDf["intensity"]);

                  queryList.emplace_back(query);
              }
              
              Rcpp::Rcout << "added " << queryList.size() << " query spectra for annotation." << std::endl;

          }

void AnnotateMs2::createRefList(Rcpp::List reference) {

    int n = reference.length();

    for (int i = 0; i < n; i++) {
        Rcpp::List ref = Rcpp::wrap(reference[i]);
        Rcpp::DataFrame info = Rcpp::wrap(ref[0]);
        Rcpp::DataFrame spectra = Rcpp::wrap(ref[1]);

        Utils util;
        Rcpp::StringVector infoValues = info["value"];
        std::vector<int> pmzPos = util.my_grep(info["key"], "precursormz");
        std::vector<int> formulaIndex = util.my_grep(info["key"], "formula");

        std::string formula;
        if (formulaIndex[0] != -1) {
            Rcpp::String val = infoValues[formulaIndex[0]];
            if (val != "NA" && val != "NULL")
                formula = Rcpp::as<std::string>(infoValues[formulaIndex[0]]);
        }
        // need to add stop if precursor is detected more than once - should not happen.
        // if (pmzPos.size() > 1) {
            
        // }
        

        double pmz = -1;
        const int pmzIndex = pmzPos[0];
        const std::string& value = Rcpp::as<std::string>(infoValues[pmzIndex]);
        if (pmzIndex != -1 && !value.empty() && value != "NA") { // Meaning the value was found
            pmz = std::stod(value);
        }
        Reference referenceData(i, pmz, formula, spectra["mz"], spectra["intensity"]);
        referenceList.emplace_back(referenceData);
    }
    
    Rcpp::Rcout << "added " << referenceList.size() << " references for annotation." << std::endl;
}

 Rcpp::DataFrame AnnotateMs2::getMatches(double threshold, const ScoringFactory& factory, const double minScore ,
     const double chemicalMinScore) {
    const int nQuery = queryList.size();
    const int nRef = referenceList.size();
    
    // int nMatches = 0;
    // for matches
    Rcpp::CharacterVector query_ms1_id;
    Rcpp::CharacterVector query_ms2_id;
    Rcpp::CharacterVector query_formula;
    Rcpp::NumericVector query_mz;
    Rcpp::NumericVector query_rt;
    Rcpp::IntegerVector ref_idx;
    Rcpp::NumericVector formulaSimilarity;
    Rcpp::NumericVector score;
    
    for (int i = 0; i < nQuery; i++) {
        Spectra currentQuerySpectra = queryList[i].GetSpectra();
        const double qPmz = currentQuerySpectra.precursorMz;
        
        for (int j = 0; j < nRef; j++) {
            if (abs(qPmz - referenceList[j].getPrecursorMz()) > threshold) {
                continue;
            }
            
            const int refIdx = referenceList[j].getIndex();
            Spectra referenceSpectra = referenceList[refIdx].GetSpectra();
            const double compScore = factory.CalculateScore(currentQuerySpectra,
                referenceSpectra, minPeaks);



            if (compScore >= minScore) {
                // nMatches += 1;
                const auto queryFormula = queryList[i].GetFormula().get_cstring();
                const auto other = referenceList[j].GetFormula().get_cstring();
                const double chemicalSimilarity = MolecularFormulaSimilarity::ComputeSimilarity(queryList[i].GetFormula(),
                referenceList[j].GetFormula());
                if (chemicalSimilarity < chemicalMinScore) continue;
                queryList[i].AddToRefComps(refIdx);
                queryList[i].AddToScores(compScore);
                formulaSimilarity.push_back(chemicalSimilarity);
                query_formula.push_back(queryList[i].GetFormula());
                query_ms1_id.push_back(queryList[i].getVariableId());
                query_ms2_id.push_back(queryList[i].getMs2Id());
                query_mz.push_back(queryList[i].getMs2Pmz());
                query_rt.push_back(queryList[i].getMs2Rt());
                ref_idx.push_back(refIdx + 1);
                score.push_back(compScore);
            }                     
        }
    }
    
    Rcpp::DataFrame m = Rcpp::DataFrame::create(Rcpp::Named("query_ms1_id") = query_ms1_id,
                                                Rcpp::Named("query_ms2_id") = query_ms2_id,
                                                Rcpp::Named("query_mz") = query_mz,
                                                Rcpp::Named("query_rt") = query_rt,
                                                Rcpp::Named("ref_idx") = ref_idx,
                                                Rcpp::Named("query_formula") = query_formula,
                                                Rcpp::Named("chemical_similarity") = formulaSimilarity,
                                                Rcpp::Named("score") = score);

    return m;
}
