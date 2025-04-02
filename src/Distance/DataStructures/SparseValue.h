#ifndef SPARSEVALUE
#define SPARSEVALUE

class SparseValue {
public:
    SparseValue(const int i, const int j, const double distance):i(i), j(j), distance(distance) {};

    // accessors
    double getDistance() const {return distance;}
    int getI() const {return i;}
    int getJ() const {return j;}

private:
    int i;
    int j;
    double distance;
};

#endif //SPARSEVALUE
// void Function()
// {
//     SparseMatrix matrix;
//     matrix.i.emplace_back(0);
//      matrix.j.emplace_back(0);
//       matrix.score.emplace_back(0);
//     std::vector<SparseMatrix> data;
//     data.emplace_back(i,j, score);
// }