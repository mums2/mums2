#ifndef ANNOTATEMS2
#define ANNOTATEMS2

#include <Rcpp.h>
#include <string>
#include <vector>
#include "DataStructures/Query.h"
#include "DataStructures/Reference.h"
#include "../ScoringMethods/ScoringFactory.h"

class AnnotateMs2 {
public:
    AnnotateMs2(const size_t minPeaks):minPeaks(minPeaks) {}
    ~AnnotateMs2() = default;
    void createQueryList(std::vector<std::string> variableId, std::vector<std::string> ms2Id,
      std::vector<float> ms2Mz, std::vector<float> ms2Rt, Rcpp::List ms2Spectra);
    void createRefList(Rcpp::List reference);
    Rcpp::DataFrame getMatches(double threshold, const ScoringFactory& factory, double minScore);



private:
    std::vector<Query> queryList;
    std::vector<Reference> referenceList;
    size_t minPeaks = 0;
};

#endif //ANNOTATEMS2