#ifndef SPARSEVALUE
#define SPARSEVALUE

struct SparseValue {
    SparseValue(const int i, const int j, const double distance):i(i), j(j), distance(distance) {};
    int i;
    int j;
    double distance;
};

#endif //SPARSEVALUE