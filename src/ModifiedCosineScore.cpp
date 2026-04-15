//
// Created by gregj on 2/21/2024.
//

#include "ScoringMethods/GNPS/ModifiedCosineScore.h"
#include "NormalizeMs2/SquareRootNormalize.h"
#include <algorithm>
ModifiedCosineScore::ModifiedCosineScore(const Rcpp::List& parameter):tolerance(parameter["tolerance"]) {

}

std::vector<double> ModifiedCosineScore::ScoreRData(const std::vector<double>& listOneMZ, std::vector<double>& listOneInt,
    const std::vector<double>& listTwoMz, std::vector<double>& listTwoInt , const double tolerance, const double shift)
{
    SquareRootNormalize::Normalize(listOneInt);
    SquareRootNormalize::Normalize(listTwoInt);
    std::vector<ScoreValues> scoreMap = ConstructPeaks(listOneMZ, listTwoMz,  listOneInt,
        listTwoInt, tolerance,0);
    std::vector<ScoreValues> shiftMap;
    if(std::abs(shift) > tolerance)
    {
        shiftMap = ConstructPeaks(listOneMZ, listTwoMz, listOneInt, listTwoInt, tolerance,shift);
    }
    scoreMap.insert(scoreMap.end(), shiftMap.begin(), shiftMap.end());
    std::sort(scoreMap.begin(), scoreMap.end(), CompareScores());
    // std::make_heap(scoreMap.begin(), scoreMap.end(), CompareScores());
    size_t peakCount = 0;
    return {ScoreMatches(scoreMap, listOneInt.size(), peakCount),
    static_cast<double>(peakCount)};
}

double ModifiedCosineScore::CalculateScore(const Spectra &firstSpectra, const Spectra &secondSpectra) {

    std::vector<double> intensityListFirst = firstSpectra.intensity;
    std::vector<double> intensityListSecond = secondSpectra.intensity;
    return ScoreRData(firstSpectra.mz, intensityListFirst, secondSpectra.mz,
        intensityListSecond, tolerance,
        (firstSpectra.precursorMz - secondSpectra.precursorMz))[0];
}

std::vector<ScoreValues> ModifiedCosineScore::ConstructPeaks(const std::vector<double>& mzVector,
    const std::vector<double>& otherMzVector, const std::vector<double>& intensitiesOne,
    const std::vector<double>& intensitiesTwo,
    const double tol, const double shift)
{
    //Declare a vector and look at its size
    const double adjTolerance = tol + 0.000001;
    // TODO: Remove the unordered maps and sets from the code and generate the priorty queue sooner.
    // TODO: This will allow us to reduce the amount of hashmap collisions and set overhead.
    std::vector<ScoreValues> scoreValuesVector;
    scoreValuesVector.reserve(mzVector.size());

    // Iterate over
    for(size_t i = 0; i < mzVector.size(); i++)
    {
        const double highBounds = mzVector[i] + adjTolerance - shift;
        const double lowBounds = mzVector[i]  - adjTolerance - shift;
        auto lowerIters = std::lower_bound(otherMzVector.begin(), otherMzVector.end(), lowBounds);
        while(lowerIters != otherMzVector.end() && *lowerIters < highBounds)
        {
            const size_t index = lowerIters - otherMzVector.begin();
            scoreValuesVector.emplace_back(ScoreValues{i, index,
                intensitiesOne[i] * intensitiesTwo[index]});
            ++lowerIters;
        }

        //Should probably construct the Priorty queue here?
        //Although im sending over a good amount of memory
    }
    return scoreValuesVector;
}

std::vector<ScoreValues> ModifiedCosineScore::ConstructPriorityQueue(std::unordered_map<size_t, std::unordered_set<size_t>>& map,
    const std::vector<double>& intensitiesOne, const std::vector<double>& intensitiesTwo, const size_t sizeOfHeapArray)
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

// double ModifiedCosineScore::ScoreMatches(
//     std::vector<ScoreValues>& queue,const size_t countOfSpectraOne, size_t& numberOfPeakMatches) {
//     std::unordered_set<size_t> usedPeakOne;
//     std::unordered_set<size_t> usedPeakTwo;
//     double totalScore = 0;
//     while(!queue.empty()) {
//         ScoreValues value = queue.front();
//         if(usedPeakOne.find(value.indexOne) != usedPeakOne.end()||
//             usedPeakTwo.find(value.indexTwo) != usedPeakTwo.end()) {
//             std::pop_heap(queue.begin(), queue.end(), CompareScores());
//             queue.pop_back();
//             continue;
//         }
//         usedPeakOne.insert(value.indexOne);
//         usedPeakTwo.insert(value.indexTwo);
//         totalScore += value.score;
//         numberOfPeakMatches ++;
//         if(numberOfPeakMatches >= countOfSpectraOne)
//             break;
//         std::pop_heap(queue.begin(), queue.end(), CompareScores());
//     }
//     return totalScore;
// }

double ModifiedCosineScore::ScoreMatches(
    std::vector<ScoreValues>& queue,const size_t countOfSpectraOne, size_t& numberOfPeakMatches) {
    std::vector<bool> usedPeakOne(countOfSpectraOne, false);
    std::vector<bool> usedPeakTwo(countOfSpectraOne,  false);
    double totalScore = 0;

    for (const auto& value : queue) {
        if(usedPeakOne[value.indexOne]||
           usedPeakTwo[value.indexTwo]) {
            continue;
           }
        usedPeakOne[value.indexOne] = true;
        usedPeakTwo[value.indexTwo] = true;
        totalScore += value.score;
        numberOfPeakMatches ++;
        if(numberOfPeakMatches >= countOfSpectraOne)
            break;
    }

    return totalScore;
}
