//
// Created by Gregory Johnson on 4/3/25.
//

#include "Spectra/ReadSpectra.h"

#include "Spectra/MetaDataValuePair.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "CustomProgressBar/CliProgressBar.h"
#include "CustomProgressBar/ETAProgressBar.h"

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
    float line_count = static_cast<float>(std::count(
         std::istream_iterator<char>(spectraData),
         std::istream_iterator<char>(),
         '\n'));
    CliProgressBar p;
    spectraData.close();
    spectraData.open(filePath);
    float currentLine = 0;
    while(std::getline(spectraData,  line)) {
        p.update(++currentLine/line_count);
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
    spectraData.close();
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
    p.end_display();
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
    float line_count = static_cast<float>(std::count(
        std::istream_iterator<char>(spectraData),
        std::istream_iterator<char>(),
        '\n'));
    CliProgressBar p;
    spectraData.close();
    spectraData.open(filePath);
    float currentLine = 0;
    auto isSpaces = [](unsigned char const c) { return std::isspace(c); };
    while(std::getline(spectraData,  line)) {
        p.update(++currentLine/line_count);
        bool whiteSpace = std::all_of(line.begin(), line.end(), isSpaces);
        if(line.empty() || whiteSpace) { // end of the current block
            mzContainer.emplace_back(mz);
            intensityContainer.emplace_back(intensity);
            metaDataKeyContainer.emplace_back(metaDataKeys);
            metaDataKeys.clear();
            mz.clear();
            intensity.clear();
            getline(spectraData, line);
            if (line.empty())
                continue;

        }
        if(line.find(':') != std::string::npos) { // headers and/or metadata
            std::vector<std::string> values;
            std::string delimiter = ": ";
            const auto pos = line.find(delimiter);
            std::string valueName = line.substr(0, pos);
            std::transform(valueName.begin(), valueName.end(), valueName.begin(), ::tolower);
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
        if (mzValue.empty() || intensityValue.empty()) continue;
        if (!std::isdigit(mzValue[0]) || !std::isdigit(intensityValue[0])) continue;
        mz.emplace_back(std::stod(mzValue));
        intensity.emplace_back(std::stod(intensityValue));


    }
    spectraData.close();
    const int spectraPeaks = static_cast<int>(mzContainer.size());
    Rcpp::DataFrame dataFrame;
    Rcpp::List mspList(spectraPeaks);
    for (int i = 0; i < spectraPeaks; i++) {
        Rcpp::List peaksDf = Rcpp::List::create(Rcpp::Named("mz") = mzContainer.front(),
            Rcpp::Named("intensity") = intensityContainer.front());

        const std::list<MetaDataValuePair> keyPairs = metaDataKeyContainer.front();
        std::unordered_map<std::string, std::vector<std::list<std::string>::iterator>> duplicateNames;
        std::list<std::string> keys;
        // std::list<std::string>::iterator k = keys.begin();
        std::list<std::string> values;
        for (const auto& value : metaDataKeyContainer.front()) {
            if (!duplicateNames[value.key].empty()) {
                // We have a duplicate, we now need to combine the values
                const auto firstFoundIndex = duplicateNames[value.key][0];
                firstFoundIndex->append(";" + value.value);
                continue;
            }
            keys.emplace_back(value.key);
            values.emplace_back(value.value);
            duplicateNames[value.key].emplace_back(--values.end());
        }
        Rcpp::List metaData = Rcpp::List::create(Rcpp::Named("key") = keys,
            Rcpp::Named("value") = values);
        mspList[i] = Rcpp::List::create(Rcpp::Named("info") = metaData,
            Rcpp::Named("spec") = peaksDf);
        intensityContainer.pop_front();
        mzContainer.pop_front();
        metaDataKeyContainer.pop_front();
    }
    p.end_display();
    return mspList;
}

Spectra ReadSpectra::ReadSpectraFile(const std::string &filePath) {
    std::ifstream spectraFile(filePath);
    if (!spectraFile.is_open()) {
        Rcpp::stop("Could not open spectra file.");
    }

    double mz, intensity;
    std::vector<double> mzValues;
    std::vector<double> intensityValues;
    mzValues.reserve(1000); // May need to read the lines
    intensityValues.reserve(1000);
    std::string line;
    while (std::getline(spectraFile, line)) {
        std::istringstream iss(line);
        iss >> mz >> intensity;
        mzValues.emplace_back(mz);
        intensityValues.emplace_back(intensity);
    }
    return Spectra {"", mzValues, intensityValues, 0};
}
