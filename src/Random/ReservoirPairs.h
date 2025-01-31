//
// Created by gregj on 1/21/2025.
//

#ifndef RESERVOIRPAIRS_H
#define RESERVOIRPAIRS_H
struct ReservoirPairs {
    double key;
    size_t value;
};
struct CompareReservoir {
    bool operator()(const ReservoirPairs& first, const ReservoirPairs& other) const {
        return first.key > other.key;
    }
};
#endif //RESERVOIRPAIRS_H
