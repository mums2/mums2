//
// Created by gregj on 2/21/2024.
//

#ifndef SCOREVALUES_H
#define SCOREVALUES_H
struct ScoreValues {
    int indexOne;
    int indexTwo;
    double score;
};
struct CompareScores {
    bool operator()(ScoreValues const& s1, ScoreValues const & s2) const {
        return s1.score < s2.score;
    }
};

#endif // SCOREVALUES_H
