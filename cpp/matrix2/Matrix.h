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

    /*
    //矩阵行列式
    double getDet() const;

#define EPSINON 1.0e-14
static inline bool isEqual(const double &num1, const double &num2)
{
    return (num1 - num2 <= EPSINON && num1 - num2 >= -EPSINON);
}

double Matrix::getDet() const
{
    //if (!isSquareMatrix())
    //    throw "Error in getDet(): the matrix is not square.";

    Matrix copyMatrix(*this);

    // 记录初等行变换的次数
    int iter = 0;

    for (int i = 0; i < m_Row; ++i)
    {
        if (isEqual(copyMatrix(i, i), 0.0))
        {
            for (int j = i; j < m_Row; ++j)
            {
                if (!isEqual(copyMatrix(j, i), 0.0))
                {
                    copyMatrix.swapRow(i, j);
                    ++iter;
                }
            }
        }

        for (int k = i + 1; k < m_Row; ++k)
        {
            double scale = -1 * copyMatrix(k, i) / copyMatrix(i, i);
            double data = 0;
            for (int u = 0; u < m_Row; ++u)
            {
                data = copyMatrix(k, u) + copyMatrix(i, u) * scale;
                copyMatrix.set(k, u, data);
            }
        }
    }

    double det = 1;
    for (int i = 0; i < m_Row; ++i)
        det *= copyMatrix(i, i);

    if (1 == iter % 2)
        det *= -1;

    return det;
}
*/

    /*
        //LU分解求逆矩阵
        Matrix inv_LU() const;
        Matrix Matrix::inv_LU() const
    {
        double det = getDet();

        //if (det >= -EPSINON && det <= EPSINON)
        //    throw "Error in inv_LU(): the determinant is equal to 0.";

        int n = m_Row;
        // 建立l、l_inverse、u、u_inverse矩阵
        Matrix l(n, n);
        Matrix l_inverse(n, n);
        Matrix u(n, n);
        Matrix u_inverse(n, n);

        // 计算矩阵对角线
        for (int i = 0; i < n; ++i)
        {
            l.set(i, i, 1.0);
            l_inverse.set(i, i, 1.0);
        }

        for (int i = 0; i < n; ++i)
        {
            for (int j = i; j < n; ++j)
            {
                double s = 0.0;
                for (int k = 0; k < i; ++k)
                    s += l(i, k) * u(k, j); // 按行计算u值
                u.set(i, j, get(i, j) - s);
            }

            for (int j = i + 1; j < n; ++j)
            {
                double s = 0.0;
                for (int k = 0; k < i; ++k)
                    s += l(j, k) * u(k, i);
                l.set(j, i, (get(j, i) - s) / u(i, i)); // 按列计算l值
            }
        }

        for (int i = 1; i < n; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                double s = 0.0;
                for (int k = 0; k < i; ++k)
                    s += l(i, k) * l_inverse(k, j);
                l_inverse.set(i, j, -s);
            }
        }

        for (int i = 0; i < n; ++i)
            u_inverse.set(i, i, 1 / u(i, i));

        for (int i = 1; i < n; ++i)
        {
            for (int j = i - 1; j >= 0; --j)
            {
                double s = 0;
                for (int k = j + 1; k <= i; ++k)
                    s += u(j, k) * u_inverse(k, i);

                u_inverse.set(j, i, -s / u(j, j));
            }
        }

        Matrix result(u_inverse * l_inverse);

        for (int i = 0; i < result.m_Row; ++i)
        {
            for (int j = 0; j < result.m_Col; ++j)
            {
                if (_isnan(result.get(i, j)))
                    throw "Error in inv_LU(): the result is 1.#INF or -1.#IND or 1.#INF000 or -1.#INF000.";
            }
        }

        return result;
    }
        */
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
