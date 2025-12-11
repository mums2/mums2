//
// Created by Gregory Johnson on 11/13/25.
//

#include "AnnotationStructure/Annotation.h"

Annotation::Annotation(const std::queue<AnnotatedNode>& nodes):
annotatedNodes(nodes){}

Rcpp::DataFrame Annotation::CreateAnnotationDataFrame(const std::vector<Feature>& features,
    const std::vector<AnnotationNode>& annotationNodes) {
    std::unordered_map<std::string, std::vector<std::string>> columnData;
    const size_t size = annotatedNodes.size();
    const std::vector<std::string> headers {"query_ms1_id", "query_ms2_id", "query_mz",
        "query_rt", "ref_idx", "query_formula", "chemical_similarity", "score"};
    for (const auto& header : headers) {
        columnData[header] = std::vector<std::string>(size);
    }

    int counter = 0;
    while (!annotatedNodes.empty()) {
        const AnnotatedNode& node = annotatedNodes.front();
        const Feature& feature = features[node.featureIndex];
        const AnnotationNode& data = annotationNodes[node.annotationNodeIndex];
        columnData["query_ms1_id"][counter] = feature.ms1_id;
        columnData["query_ms2_id"][counter] = feature.ms2_id;
        columnData["query_mz"][counter] = std::to_string(feature.mz);
        columnData["query_rt"][counter] = std::to_string(feature.rt);
        columnData["ref_idx"][counter] = std::to_string(data.referenceIndex);
        columnData["query_formula"][counter] = feature.formula;
        columnData["chemical_similarity"][counter] = std::to_string(node.formulaSimilarity);
        columnData["score"][counter] = std::to_string(node.score);
        for (const auto& keyValue : data.keyValues) {
            if (columnData.find(keyValue.key) == columnData.end()) {
                std::vector<std::string> values(size);
                values[counter] = keyValue.value;
                columnData[keyValue.key] = values;
                continue;
            }
            columnData[keyValue.key][counter] = keyValue.value;
        }
        annotatedNodes.pop();
        counter++;
    }
    Rcpp::DataFrame df = Rcpp::DataFrame::create();
    for (const auto& map : columnData) {
        df.push_back(map.second, map.first);
    }
    return df;
}
