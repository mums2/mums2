//
// Created by Gregory Johnson on 4/3/25.
//

#include "Spectra/ReadSpectra.h"

#include "Spectra/MetaDataValuePair.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

Rcpp::List ReadSpectra::ReadMGF(const std::string& filePath) {
    std::ifstream spectraData(filePath);
    std::string line;
    std::list<std::list<double>> mzContainer;
    std::list<std::list<double>> intensityContainer;
    std::list<double> mz;
    std::list<double> intensity;
    std::unordered_map<std::string, std::vector<std::string>> map;

    spectraData.unsetf(std::ios_base::skipws);
    // count the newlines with an algorithm specialized for counting:
    unsigned long line_count = std::count(
        std::istream_iterator<char>(spectraData),
        std::istream_iterator<char>(),
        '\n');
    Progress p(line_count, true);
    spectraData.clear();
    spectraData.seekg(0);
    spectraData.setf(std::ios_base::skipws);
    while(std::getline(spectraData,  line)) {
        p.increment();
        if(line.find("BEGIN IONS") != std::string::npos) continue;
        if(line.find("END IONS") != std::string::npos) { // ending

            mzContainer.emplace_back(mz);
            intensityContainer.emplace_back(intensity);
            mz.clear();
            intensity.clear();
            continue;
        }
        if(line.find("PEPMASS") != std::string::npos || line.find("RTINMINUTES") != std::string::npos ||
            line.find("RTINSECONDS") != std::string::npos) { // headers and/or metadata
            std::vector<std::string> values;
            std::string delimiter = "=";
            const auto pos = line.find(delimiter);
            const std::string valueName = line.substr(0, pos);
            const std::string value = line.substr(pos + delimiter.length(), line.size());
            map[valueName].emplace_back(value);
            continue;
        }
        if(line.find('=') == std::string::npos && line.find(' ') != std::string::npos) { // peak mz/intensities
            std::vector<std::string> values;
            std::string delimiter = " ";
            const auto pos = line.find(delimiter);
            const std::string mzValue = line.substr(0, pos);
            const std::string intensityValue = line.substr(pos + delimiter.length(), line.size());
            mz.emplace_back(std::stod(mzValue));
            intensity.emplace_back(std::stod(intensityValue));
        }

    }
    const int spectraPeaks = static_cast<int>(mzContainer.size());
    Rcpp::DataFrame dataFrame;
    for (const auto& value : map) {
        dataFrame.push_back(value.second, value.first);
    }
    Rcpp::List mzIntensityList(spectraPeaks);
    for (int i = 0; i < spectraPeaks; i++) {
        mzIntensityList[i] = Rcpp::List::create(Rcpp::Named("mz") = mzContainer.front(),
        Rcpp::Named("intensity") = intensityContainer.front());
        intensityContainer.pop_front();
        mzContainer.pop_front();
    }
    return Rcpp::List::create(Rcpp::Named("ms2_table") = dataFrame,
        Rcpp::Named("mzIntensityList") = mzIntensityList);
}


Rcpp::List ReadSpectra::ReadMSP(const std::string& filePath) {
    std::ifstream spectraData(filePath);
    std::string line;
    std::list<std::list<double>> mzContainer;
    std::list<std::list<double>> intensityContainer;
    std::list<double> mz;
    std::list<double> intensity;
    std::list<std::list<MetaDataValuePair>> metaDataKeyContainer;
    std::list<MetaDataValuePair> metaDataKeys;
    spectraData.unsetf(std::ios_base::skipws);
    // count the newlines with an algorithm specialized for counting:
    unsigned long line_count = std::count(
        std::istream_iterator<char>(spectraData),
        std::istream_iterator<char>(),
        '\n');
    Progress p(line_count, true);
    spectraData.clear();
    spectraData.seekg(0);
    spectraData.setf(std::ios_base::skipws);
    while(std::getline(spectraData,  line)) {
        p.increment();
        if(line.empty()) { // end of the current block
            mzContainer.emplace_back(mz);
            intensityContainer.emplace_back(intensity);
            metaDataKeyContainer.emplace_back(metaDataKeys);
            metaDataKeys.clear();
            mz.clear();
            intensity.clear();
            long long len = spectraData.tellg();
            getline(spectraData, line);
            if (!line.empty())
                spectraData.seekg(len ,std::ios_base::beg);
            continue;
        }
        if(line.find(':') != std::string::npos) { // headers and/or metadata
            std::vector<std::string> values;
            std::string delimiter = ": ";
            const auto pos = line.find(delimiter);
            std::string valueName = line.substr(0, pos);
            std::string value = line.substr(pos + delimiter.length(), line.size());
            metaDataKeys.emplace_back(valueName, value);
            continue;
        }

        if(line.find('\t') == std::string::npos && line.find(' ') == std::string::npos) continue;
        std::string delimiter = "\t";
        if (line.find(' ') != std::string::npos)
            delimiter = ' ';
        // peak mz/intensities
        std::vector<std::string> values;
        const auto pos = line.find(delimiter);
        const std::string mzValue = line.substr(0, pos);
        const std::string intensityValue = line.substr(pos + delimiter.length(), line.size());
        mz.emplace_back(std::stod(mzValue));
        intensity.emplace_back(std::stod(intensityValue));

    }
    const int spectraPeaks = static_cast<int>(mzContainer.size());
    Rcpp::DataFrame dataFrame;
    Rcpp::List mspList(spectraPeaks);

    for (int i = 0; i < spectraPeaks; i++) {
        Rcpp::DataFrame peaksDf = Rcpp::DataFrame::create(Rcpp::Named("mz") = mzContainer.front(),
            Rcpp::Named("intensity") = intensityContainer.front());

        size_t counter = 0;
        const std::list<MetaDataValuePair> keyPairs = metaDataKeyContainer.front();
        const size_t size = keyPairs.size();
        std::vector<std::string> keys(size);
        std::vector<std::string> values(size);
        for (const auto& value : metaDataKeyContainer.front()) {
            keys[counter] = value.key;
            values[counter++] = value.value;
        }
        Rcpp::DataFrame metaData = Rcpp::DataFrame::create(Rcpp::Named("key") = keys,
            Rcpp::Named("value") = values);
        mspList[i] = Rcpp::List::create(Rcpp::Named("info") = metaData,
            Rcpp::Named("spec") = peaksDf);
        intensityContainer.pop_front();
        mzContainer.pop_front();
        metaDataKeyContainer.pop_front();
    }
    return mspList;
}
