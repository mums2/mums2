//
// Created by Gregory Johnson on 4/3/25.
//

#include "ReadSpectra.h"

bool ReadSpectra::ReadMGF(std::string filePath) {
    std::ifstream spectraData(filePath);
    std::string line;
    std::list<double> mz;
    std::list<double> intensity;
    SpectraTable spectraTable;
    std::list<SpectraTable> spectraTables;
    std::unordered_map<std::string, std::string> map;
    while(std::getline(spectraData,  line)) {
        if(line.find("BEGIN IONS") != std::string::npos) continue;
        if(line.find("END IONS") != std::string::npos) { // ending
            std::vector<double> tmpMz;
            tmpMz.insert(tmpMz.end(), std::make_move_iterator(mz.begin()), std::make_move_iterator(mz.end()));
            std::vector<double> tmpInt;
            tmpInt.insert(tmpInt.end(), std::make_move_iterator(intensity.begin()), std::make_move_iterator(intensity.end()));
            spectraTables.emplace_front(MapUnordedMapToSpectraTable(map,tmpMz, tmpInt));
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
            map[valueName] = value;
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
}
