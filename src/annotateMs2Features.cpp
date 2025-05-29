#include "AnnotateMs2/AnnotateMs2.h"
#include "Utils/Utils.h"
#include <Rcpp.h>

#include "AnnotateMs2/DataStructures/AnnotationNode.h"
#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"
#include "ScoringMethods/ScoringFactory.h"

// [[Rcpp::export]]
Rcpp::DataFrame AnnotateMs2Features(const std::vector<std::string>& variableId, const std::vector<std::string>& ms2Id, 
          const std::vector<float>& ms2Mz, const std::vector<float>& ms2Rt, const Rcpp::StringVector& formulas,
          const Rcpp::List& ms2Spectra, const Rcpp::List& reference, const Rcpp::List& parameters,
          const double precursorThreshold, const double minScore, const double chemicalMinScore,
          const size_t minPeaks) {

    const ScoringFactory factory(parameters);
    AnnotateMs2 annotateMs2(minPeaks);
    annotateMs2.createQueryList(variableId, ms2Id, ms2Mz, ms2Rt, formulas, ms2Spectra);
    annotateMs2.createRefList(reference);
    Rcpp::DataFrame matches = annotateMs2.getMatches(precursorThreshold, factory, minScore, chemicalMinScore);
    return matches;

    
}


// [[Rcpp::export]]
Rcpp::DataFrame AnnotateMs2Features2(const Rcpp::DataFrame& queryList, const Rcpp::List querySpectra,
    const Rcpp::List referenceList, const Rcpp::List& scoringParameters, const Rcpp::StringVector& formulas,
    const double precursorThreshold,const double minScoreThreshold, const double chemicalMinScore,
    const size_t minPeaks) {
    const ScoringFactory factory(scoringParameters);
    const auto querySpectraSize = static_cast<size_t>(querySpectra.size());
    const auto referenceListSize = static_cast<size_t>(referenceList.size());
    std::vector<AnnotationNode> querySpectaList(querySpectraSize);
    std::vector<AnnotationNode> referenceSpectraList(referenceListSize);
    for (size_t i = 0; i < querySpectraSize; ++i) {
        const Rcpp::List& spectra = querySpectra[i];
        const Rcpp::NumericVector& queryData = queryList["mz"];
        querySpectaList[i] = AnnotationNode(spectra["mz"], spectra["intensity"],
            formulas[i], queryData[i], i);
    }

    Utils utils;
    for (size_t i = 0; i < referenceListSize; ++i) {
        const Rcpp::List& ref = referenceList[i];
        const Rcpp::List& info = ref["info"];
        const Rcpp::List& spec = ref["spec"];
        Rcpp::StringVector values = info["value"];
        Rcpp::String formula;
        std::string mz = "-1";
        std::string a;
        const std::vector<int> indexes = utils.my_grep(info["key"], "precursormz");
        const std::vector<int> formulaIndex = utils.my_grep(info["key"], "formula");

        if (indexes[0] != -1) {
            std::string value = Rcpp::as<std::string>(values[indexes[0]]);
            if (value != "NA" && value != "NULL")
                mz = Rcpp::as<std::string>(values[indexes[0]]);
        }


        if (formulaIndex[0] != -1) {
            std::string value = Rcpp::as<std::string>(values[formulaIndex[0]]);
            if (value != "NA" && value != "NULL")
                formula = values[formulaIndex[0]];
        }

        referenceSpectraList[i] = AnnotationNode(spec["mz"],
            spec["intensity"], formula, std::stod(mz), i);
    }
    const Rcpp::StringVector& ms1Ids = queryList["ms1_compound_id"];
    const Rcpp::StringVector& ms2Ids = queryList["ms2_spectrum_id"];
    const Rcpp::StringVector& queryRts = queryList["rt"];

    std::list<std::string> query_ms1_id;
    std::list<std::string> query_ms2_id;
    std::list<std::string> query_formula;
    std::list<double> query_mz;
    std::list<std::string> query_rt;
    std::list<int> ref_idx;
    std::list<double> formulaSimilarity;
    std::list<double> score;
    for (const auto& query : querySpectaList) {
        double currentPrecursorMz = query.GetPrecursorMz();
        Spectra currentSpectra = query.GetSpectra();
        int index = query.GetIndex();
        for (const auto& ref : referenceSpectraList) {
            if ((std::abs(currentPrecursorMz - ref.GetPrecursorMz())) * 1e6 / currentPrecursorMz > precursorThreshold) continue;
            const double chemicalComparison = MolecularFormulaSimilarity::ComputeSimilarity(query.GetFormula(), ref.GetFormula());
            if (chemicalComparison < chemicalMinScore) continue;
            const double comparisonScore = factory.CalculateScore(currentSpectra, ref.GetSpectra(), minPeaks);
            if (comparisonScore < minScoreThreshold) continue;
            query_ms1_id.emplace_back(ms1Ids[index]);
            query_ms2_id.emplace_back(ms2Ids[index]);
            query_formula.emplace_back(formulas[index]);
            query_mz.emplace_back(query.GetPrecursorMz());
            query_rt.emplace_back(queryRts[index]);
            ref_idx.emplace_back(ref.GetIndex());
            formulaSimilarity.emplace_back(chemicalComparison);
            score.emplace_back(comparisonScore);
        }
    }


    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::Named("query_ms1_id") = query_ms1_id,
                                                Rcpp::Named("query_ms2_id") = query_ms2_id,
                                                Rcpp::Named("query_mz") = query_mz,
                                                Rcpp::Named("query_rt") = query_rt,
                                                Rcpp::Named("ref_idx") = ref_idx,
                                                Rcpp::Named("query_formula") = query_formula,
                                                Rcpp::Named("chemical_similarity") = formulaSimilarity,
                                                Rcpp::Named("score") = score);

    std::unordered_map<std::string, Rcpp::CharacterVector> uniqueColumns;
    int rowSize = df.nrow();
    int count = 0;
    for (const auto index : ref_idx) {
        const Rcpp::List& ref = referenceList[index];
        const Rcpp::List& info = ref["info"];
        Rcpp::StringVector values = info["value"];
        Rcpp::CharacterVector headers = info["key"];
        for (int j = 0; j < headers.size(); j++) {
            std::string header = Rcpp::as<std::string>(headers[j]);
            if (uniqueColumns.find(header) == uniqueColumns.end()) {
                Rcpp::CharacterVector newColumn(rowSize);
                uniqueColumns[header] = newColumn;
            }
            uniqueColumns[header][count] = values[j];
        }
        count++;
    }
    for (const auto& value : uniqueColumns) {
        df.push_back(value.second, value.first);
    }
    return df;
}