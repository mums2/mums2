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
    explicit AnnotateMs2(const size_t minPeaks):minPeaks(minPeaks) {}
    ~AnnotateMs2() = default;
    void createQueryList(const std::vector<std::string>& variableId, const std::vector<std::string> &ms2Id,
          const std::vector<float> &ms2Mz, const std::vector<float> &ms2Rt, const Rcpp::StringVector& formulas,
          const Rcpp::List& ms2Spectra);
    void createRefList(Rcpp::List reference);
    Rcpp::DataFrame getMatches(double threshold, const ScoringFactory& factory, double minScore,
        double chemicalMinScore);



private:
    std::vector<Query> queryList;
    std::vector<Reference> referenceList;
    size_t minPeaks = 0;
};

#endif //ANNOTATEMS2