//
// Created by Gregory Johnson on 11/13/25.
//

#ifndef MUMS2_ANNOTATIONNODE_H
#define MUMS2_ANNOTATIONNODE_H
#include <string>
#include <list>
#include "../Distance/DataStructures/Spectra.h"

struct KeyValues {
    std::string key;
    std::string value;
};
struct AnnotationNodeData {
    double precursorMz;
    std::string name;
    std::string chemicalFormula;
    std::list<KeyValues> keyValues;
    std::list<std::string> keys;
    std::list<std::string> values;
    Spectra spectra;
    size_t referenceIndex;
};
#endif //MUMS2_ANNOTATIONNODE_H