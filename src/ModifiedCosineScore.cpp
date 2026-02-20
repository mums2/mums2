//
// Created by gregj on 2/21/2024.
//

#include "ScoringMethods/GNPS/ModifiedCosineScore.h"
#include <algorithm>
ModifiedCosineScore::ModifiedCosineScore(const Rcpp::List& parameter):tolerance(parameter["tolerance"]) {

}

std::vector<double> ModifiedCosineScore::ScoreRData(const std::vector<double>& listOneMZ, std::vector<double>& listOneInt,
    const std::vector<double>& listTwoMz, std::vector<double>& listTwoInt , const double tolerance, const double shift)
{
    SquareRootNormalize::Normalize(listOneInt);
    SquareRootNormalize::Normalize(listTwoInt);
    int sizeOfMap = 0;
    int shiftMapSize = 0;
    auto map = ConstructPeaks(listOneMZ, listTwoMz, tolerance,0, sizeOfMap);
    std::unordered_map<int, std::unordered_set<int>> shift_map;
    if(std::abs(shift) > tolerance)
    {
        shift_map = ConstructPeaks(listOneMZ, listTwoMz, tolerance,shift, shiftMapSize);
    }
    auto shiftQueue = ConstructPriorityQueue(shift_map, listOneInt, listTwoInt, shiftMapSize);
    auto queue = ConstructPriorityQueue(map, listOneInt, listTwoInt, sizeOfMap);
    queue.insert(queue.end(), shiftQueue.begin(), shiftQueue.end());
    std::make_heap(queue.begin(), queue.end(), CompareScores());
    int peakCount = 0;
    return {ScoreMatches(queue, listOneInt.size(), peakCount),
    static_cast<double>(peakCount)};
}

double ModifiedCosineScore::CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) {

    std::vector<double> intensityListFirst = firstSpectra.intensity;
    std::vector<double> intensityListSecond = secondSpectra.intensity;
    return ScoreRData(firstSpectra.mz, intensityListFirst, secondSpectra.mz,
        intensityListSecond, tolerance,
        (firstSpectra.precursorMz - secondSpectra.precursorMz))[0];
}

std::unordered_map<int, std::unordered_set<int>> ModifiedCosineScore::ConstructPeaks(const std::vector<double>& mzVector,
    const std::vector<double>& otherMzVector, const double tol, const double shift, int& totalSizeOfIndexMap)
{
    //Declare a vector and look at its size
    const double adjTolerance = tol + 0.000001;
    std::unordered_map<int, std::unordered_set<int>> indexMap;
    totalSizeOfIndexMap = 0;
    // Iterate over
    for(size_t i = 0; i < mzVector.size(); i++)
    {
        const double highBounds = mzVector[i] + adjTolerance - shift;
        const double lowBounds = mzVector[i]  - adjTolerance - shift;
        auto lowerIters = std::lower_bound(otherMzVector.begin(), otherMzVector.end(), lowBounds);
        while(lowerIters != otherMzVector.end() && *lowerIters < highBounds)
        {
            indexMap[i].insert(lowerIters - otherMzVector.begin());
            ++lowerIters;
            totalSizeOfIndexMap++;
        }
        //Should probably construct the Priorty queue here?
        //Although im sending over a good amount of memory
    }
    return indexMap;
}

std::vector<ScoreValues> ModifiedCosineScore::ConstructPriorityQueue(std::unordered_map<int, std::unordered_set<int>>& map,
    const std::vector<double>& intensitiesOne, const std::vector<double>& intensitiesTwo, int sizeOfHeapArray)
{
    //std::priority_queue<ScoreValues, std::vector<ScoreValues>, CompareScores> queue;
    std::vector<ScoreValues> score_values_vector(sizeOfHeapArray);
    int count = 0;
    for(auto& value : map) {
        const double intensityOneValue = intensitiesOne[value.first];
        for(const auto index : value.second) {
            const double intensityTwoValue = intensitiesTwo[index];
            const double score = intensityOneValue * intensityTwoValue;
            score_values_vector[count++] = {value.first, index,score};
        }
    }
    return score_values_vector;
}

double ModifiedCosineScore::ScoreMatches(
    std::vector<ScoreValues>& queue,const size_t countOfSpectraOne, int& numberOfPeakMatches) {
    std::unordered_set<int> usedPeakOne;
    std::unordered_set<int> usedPeakTwo;
    double totalScore = 0;
    while(!queue.empty()) {
        ScoreValues value = queue.front();
        if(usedPeakOne.find(value.indexOne) != usedPeakOne.end()||
            usedPeakTwo.find(value.indexTwo) != usedPeakTwo.end()) {
            std::pop_heap(queue.begin(), queue.end(), CompareScores());
            queue.pop_back();
            continue;
        }
        usedPeakOne.insert(value.indexOne);
        usedPeakTwo.insert(value.indexTwo);
        totalScore += value.score;
        numberOfPeakMatches ++;
        if(numberOfPeakMatches >= static_cast<int>(countOfSpectraOne))
            break;
        std::pop_heap(queue.begin(), queue.end(), CompareScores());
    }
    return totalScore;
}
