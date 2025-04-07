//
// Created by Gregory Johnson on 4/3/25.
//

#include "Spectra/ReadSpectra.h"

Rcpp::List ReadSpectra::ReadMGF(const std::string& filePath) {
    std::ifstream spectraData(filePath);
    std::string line;
    std::list<std::list<double>> mzContainer;
    std::list<std::list<double>> intensityContainer;
    std::list<double> mz;
    std::list<double> intensity;
    std::unordered_map<std::string, std::vector<std::string>> map;
    while(std::getline(spectraData,  line)) {
        if(line.find("BEGIN IONS") != std::string::npos) continue;
        if(line.find("END IONS") != std::string::npos) { // ending
            // std::vector<double> tmpMz;
            // tmpMz.insert(tmpMz.end(), std::make_move_iterator(mz.begin()), std::make_move_iterator(mz.end()));
            // std::vector<double> tmpInt;
            // tmpInt.insert(tmpInt.end(), std::make_move_iterator(intensity.begin()), std::make_move_iterator(intensity.end()));
            // spectraTables.emplace_front(MapUnordedMapToSpectraTable(map,tmpMz, tmpInt));
            mzContainer.emplace_back(mz);
            intensityContainer.emplace_back(intensity);
            mz.clear();
            intensity.clear();
            continue;
        }
        if(line.find('=') != std::string::npos) { // peak mz/intensities
            std::vector<std::string> values;
            std::string delimiter = "=";
            const auto pos = line.find(delimiter);
            const std::string valueName = line.substr(0, pos);
            const std::string value = line.substr(pos + delimiter.length(), line.size());
            map[valueName].emplace_back(value);
            continue;
        }
        if(line.find(' ') != std::string::npos) { // peak mz/intensities
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
        mzIntensityList[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = mzContainer.front(),
        Rcpp::Named("intensity") = intensityContainer.front());
        intensityContainer.pop_front();
        mzContainer.pop_front();
    }
    return Rcpp::List::create(Rcpp::Named("ms2_table") = dataFrame,
        Rcpp::Named("mzIntensityList") = mzIntensityList);
}
