//
// Created by gregj on 7/26/2025.
//

#ifndef CPPMATRIX_H
#define CPPMATRIX_H
#include <vector>
#include <Rcpp.h>
class CppMatrix {
public:
    explicit CppMatrix(const std::vector<std::vector<double>>& matrix) {
        const size_t rows = matrix.size();
        const size_t cols = matrix[0].size();
        rowSize = rows;
        colSize = cols;
        data = std::vector<double>(rows * cols);
        size_t counter = 0;
        for (const auto& row : matrix) {
            for (const auto& value : row) {
                data[counter++] = value;
            }
        }
    };
    CppMatrix() = default;
    ~CppMatrix() = default;
    CppMatrix operator+(const CppMatrix& other) const {
        const size_t size = data.size();
        if (other.data.size() != size)
            Rcpp::stop("CppMatrix::operator+, matrices are of different sizes");
        std::vector<double> result(size);
        for (size_t i = 0; i < size; i++) {
            result[i] = data[i] + other.data[i];
        }
        return CppMatrix(result, rowSize, colSize);
    }
    void operator+=(const CppMatrix& other) {
        const size_t size = data.size();
        if (other.data.size() != size)
            Rcpp::stop("CppMatrix::operator+=, matrices are of different sizes");
        for (size_t i = 0; i < size; i++) {
            data[i] += other.data[i];
        }
    }

    CppMatrix operator-(const CppMatrix& other) const {
        const size_t size = data.size();
        if (other.data.size() != size)
            Rcpp::stop("CppMatrix::operator-, matrices are of different sizes");
        std::vector<double> result(size);
        for (size_t i = 0; i < size; i++) {
            result[i] = data[i] - other.data[i];
        }
        return CppMatrix(result, rowSize, colSize);
    }

    void operator/=(const double other) {
        const size_t size = data.size();
        for (size_t i = 0; i < size; i++) {
            data[i] /= other;
        }
    }


    bool operator==(const CppMatrix& other) const {
        const size_t size = data.size();
        if (other.data.size() != size)
            Rcpp::stop("CppMatrix::operator==, matrices are of different sizes");
        for (size_t i = 0; i < size; i++) {
            if (data[i] != other.data[i]) return false;
        }
        return true;
    }
    std::vector<double> GetRow(const size_t rowIndex) const{
        std::vector<double> result(colSize);
        for (size_t i = 0; i < colSize; i++) {
            result[i] = data[rowIndex * colSize + i];
        }
        return result;
    }

    Rcpp::NumericMatrix ToRcppMatrix() const {
        Rcpp::NumericMatrix result(rowSize, colSize);
        for (size_t i = 0; i < rowSize; i++) {
            for (size_t j = 0; j < colSize; j++) {
                result(i, j) = data[i * colSize + j];
            }
        }
        return result;
    }
    size_t GetRowSize() const {return rowSize;}
    size_t GetColSize() const {return colSize;}

    explicit CppMatrix(const std::vector<double>& data, const size_t rows, const size_t cols) : data(data),
    rowSize(rows), colSize(cols){}

private:
    std::vector<double> data;
    size_t rowSize;
    size_t colSize;
};
#endif //CPPMATRIX_H
