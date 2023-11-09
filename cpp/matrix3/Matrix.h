#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <cassert>
#include <ctime>
using namespace std;

class Matrix
{
public:
    Matrix(int row = 0, int col = 0, double init = 0.0);
    Matrix(const Matrix &m);
    Matrix &operator=(const Matrix &m);
    Matrix(Matrix &&m);

public:
    inline int row() const { return m_Row; }
    inline int col() const { return m_Col; }
    inline shared_ptr<double[]> data() const { return m_Data; }
    inline int numel() const { return m_Col * m_Row; }
    inline bool squareMatrix() const { return m_Row == m_Col; }
    inline int length() const { return m_Row > m_Col ? m_Row : m_Col; }
    inline bool valid() const { return m_Data != nullptr; }

public:
    inline double &operator()(int row, int col) { return m_Data[m_Col * row + col]; }
    inline double operator()(int row, int col) const { return m_Data[m_Col * row + col]; }
    inline double &operator[](int number) { return m_Data[number]; }
    inline double operator[](int number) const { return m_Data[number]; }

public:
    static Matrix zeros(int row, int col) { return Matrix(row, col); }
    static Matrix ones(int row, int col) { return Matrix(row, col, 1.0); }
    static Matrix diag(int n, double init = 1.0)
    {
        Matrix result(n, n);
        for (int i = 0; i < n; ++i)
            result(i, i) = init;
        return result;
    }
    static Matrix randMatrix(int row, int col, double start, double end)
    {
        srand((unsigned)time(NULL));

        Matrix result(row, col);
        for (int i = 0; i < row * col; ++i)
            result[i] = std::rand() / double(RAND_MAX) * (end - start) + start;
        return result;
    }

public:
    double sum() const;
    Matrix T() const;
    Matrix inv() const;

    friend Matrix cat(const Matrix &m1, const Matrix &m2);
    friend Matrix diag(const Matrix &m);
    friend Matrix sqrt(const Matrix &m);
    friend Matrix sum(const Matrix &m, int dim);
    friend Matrix mul(const Matrix &m1, const Matrix &m2);
    friend Matrix operator+(const Matrix &m1, const Matrix &m2);
    friend Matrix operator-(const Matrix &m1, const Matrix &m2);
    friend Matrix operator*(double data, const Matrix &m);
    friend Matrix operator*(const Matrix &m1, const Matrix &m);
    friend ostream &operator<<(ostream &_out, const Matrix &m);

private:
    void swapRow(int iRow, int jRow);

private:
    int m_Row{0};
    int m_Col{0};
    shared_ptr<double[]> m_Data{nullptr};
};
