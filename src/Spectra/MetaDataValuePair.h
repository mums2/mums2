//
// Created by gregj on 4/8/2025.
//

#ifndef METADATAVALUEPAIR_H
#define METADATAVALUEPAIR_H
#include <string>

struct MetaDataValuePair {
    MetaDataValuePair(std::string& key, std::string& value) : key(std::move(key)), value(std::move(value)) {}
    std::string key;
    std::string value;
};
#endif //METADATAVALUEPAIR_H
