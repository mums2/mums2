//
// Created by Gregory Johnson on 11/13/25.
//

#include "AnnotationStructure/Annotation.h"

Annotation::Annotation(const std::queue<AnnotatedNode>& nodes):
annotatedNodes(nodes){}

Rcpp::DataFrame Annotation::CreateAnnotationDataFrame() {
    std::unordered_map<std::string, std::vector<std::string>> columnData;
    const size_t size = annotatedNodes.size();
    int counter = 0;
    while (!annotatedNodes.empty()) {
        const AnnotatedNode& node = annotatedNodes.front();
        columnData["query_ms1_id"].emplace_back(node.feature.ms1_id);
        columnData["query_ms2_id"].emplace_back(node.feature.ms2_id);
        columnData["query_mz"].emplace_back(std::to_string(node.feature.mz));
        columnData["query_rt"].emplace_back(std::to_string(node.feature.rt));
        columnData["ref_idx"].emplace_back(std::to_string(node.node.referenceIndex));
        columnData["query_formula"].emplace_back(node.feature.formula);
        columnData["chemical_similarity"].emplace_back(std::to_string(node.formulaSimilarity));
        columnData["score"].emplace_back(std::to_string(node.score));
        const AnnotationNode& data = node.node;
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
