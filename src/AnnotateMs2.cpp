#include "AnnotateMs2/AnnotateMs2.h"
#include "Utils/Utils.h"
#include <Rcpp.h>


void AnnotateMs2::createQueryList(std::vector<std::string> variableId, std::vector<std::string> ms2Id, 
          std::vector<float> ms2Mz, std::vector<float> ms2Rt, Rcpp::List ms2Spectra) {
              
              int n = ms2Id.size();
              for (int i = 0; i < n; i++) {
                  Rcpp::DataFrame specDf = Rcpp::wrap(ms2Spectra[i]);

                  Query query(variableId[i], ms2Id[i], ms2Mz[i], ms2Rt[i], specDf["mz"], specDf["intensity"]);
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
        std::vector<int> pmzPos = util.my_grep(info["key"], "precursormz");
        
        // need to add stop if precursor is detected more than once - should not happen.
        // if (pmzPos.size() > 1) {
            
        // }
        
        std::vector<std::string> infoValues = Rcpp::as<std::vector<std::string>>(info["value"]);
        double pmz = std::stod(infoValues[pmzPos[0]]);
        
        Reference reference(i, pmz, spectra["mz"], spectra["intensity"]);
        referenceList.emplace_back(reference); 
    }
    
    Rcpp::Rcout << "added " << referenceList.size() << " references for annotation." << std::endl;
}

 Rcpp::DataFrame AnnotateMs2::getMatches(double threshold, const ScoringFactory& factory, double minScore) {
    int nQuery = queryList.size();
    int nRef = referenceList.size();
    
    int nMatches = 0;
    // for matches
    Rcpp::CharacterVector query_ms1_id;
    Rcpp::CharacterVector query_ms2_id;
    Rcpp::NumericVector query_mz;
    Rcpp::NumericVector query_rt;
    Rcpp::IntegerVector ref_idx;
    Rcpp::NumericVector score;
    
    for (int i = 0; i < nQuery; i++) {
        Spectra currentQuerySpectra = queryList[i].GetSpectra();
        double qPmz = currentQuerySpectra.precursorMz;
        
        for (int j = 0; j < nRef; j++) {
            if (abs(qPmz - referenceList[j].getPrecursorMz()) > threshold) {
                continue;
            }
            
            int refIdx = referenceList[j].getIndex();
            Spectra referenceSpectra = referenceList[refIdx].GetSpectra();
            double compScore = factory.CalculateScore(currentQuerySpectra, referenceSpectra);
            
            if (compScore >= minScore) {
                nMatches += 1;
                queryList[i].AddToRefComps(refIdx);
                queryList[i].AddToScores(compScore);
                
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
                                                Rcpp::Named("score") = score);

    return m;
}
