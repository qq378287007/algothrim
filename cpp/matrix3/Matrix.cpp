#include "Matrix.h"

Matrix::Matrix(int row, int col, double init)
{
    if (row > 0 && col > 0)
    {
        m_Row = row;
        m_Col = col;
        m_Data.reset(new double[m_Row * m_Col]);
        for (int i = 0; i < m_Row * m_Col; ++i)
            m_Data[i] = init;
    }
}

Matrix::Matrix(const Matrix &m)
{
    if (!m.valid())
        return;

    m_Row = m.row();
    m_Col = m.col();
    m_Data.reset(new double[m_Row * m_Col]);
    for (int row = 0; row < m_Row; ++row)
        for (int col = 0; col < m_Col; ++col)
            operator()(row, col) = m(row, col);
}

Matrix &Matrix::operator=(const Matrix &m)
{
    if (this != &m)
    {
        if (m.valid())
        {
            m_Row = m.row();
            m_Col = m.col();
            m_Data.reset(new double[m_Row * m_Col]);
            for (int row = 0; row < m_Row; ++row)
                for (int col = 0; col < m_Col; ++col)
                    operator()(row, col) = m(row, col);
        }
        else
        {
            m_Row = 0;
            m_Col = 0;
            m_Data = nullptr;
        }
    }
    return *this;
}

Matrix::Matrix(Matrix &&m)
{
    m_Row = m.m_Row;
    m_Col = m.m_Col;
    m_Data = m.m_Data;
    m.m_Row = 0;
    m.m_Col = 0;
    m.m_Data = nullptr;
}

double Matrix::sum() const
{
    double data = 0.0;
    for (int i = 0; i < numel(); i++)
        data += operator[](i) * operator[](i);
    return data;
}

Matrix Matrix::T() const
{
    Matrix result(m_Col, m_Row);
    for (int row = 0; row < m_Col; ++row)
        for (int col = 0; col < m_Row; ++col)
            result(row, col) = operator()(col, row);
    return result;
}

// 高斯消元法求逆矩阵
Matrix Matrix::inv() const
{
    assert(m_Row == m_Col);

    Matrix I(*this);
    Matrix result = Matrix::diag(m_Row, 1.0);

    for (int row = 0; row < m_Row; row++)
    {
        // 选主元
        int pivotRow = row;
        for (int cur_row = row + 1; cur_row < m_Row; cur_row++)
        {
            if (std::abs(I(cur_row, row)) > std::abs(I(pivotRow, row)))
                pivotRow = cur_row;
        }

        I.swapRow(row, pivotRow);
        result.swapRow(row, pivotRow);

        // 每行第一个系数归一化
        double tmp = I(row, row);
        /*
        for (int col = 0; col < m_Col; col++) {
            I(row, col) /= tmp;
            result(row, col) /= tmp;
        }
        */
        for (int col = row; col < m_Col; col++) // I中col从row开始就行
            I(row, col) /= tmp;
        for (int col = 0; col < m_Col; col++)
            result(row, col) /= tmp;

        // 按行消元
        for (int cur_row = row + 1; cur_row < m_Row; cur_row++)
        {
            double tmp = I(cur_row, row);
            /*
            for (int col = 0; col < m_Col; col++)
            {
                I(cur_row, col) -= I(row, col) * tmp;
                result(cur_row, col) -= result(row, col) * tmp;
            }
            */
            for (int col = row; col < m_Col; col++) // I中col从row开始就行
                I(cur_row, col) -= I(row, col) * tmp;
            for (int col = 0; col < m_Col; col++)
                result(cur_row, col) -= result(row, col) * tmp;
        }
    }

    for (int row = m_Row - 1; row >= 1; row--)
    {
        for (int cur_row = row - 1; cur_row >= 0; cur_row--)
        {
            double tmp = I(cur_row, row);
            /*
            for (int col = m_Col-1; col >=0 ; col--)
            {
                I(cur_row, col) -= I(row, col) * tmp;
                result(cur_row, col) -= result(row, col) * tmp;
            }
            */
            I(cur_row, row) -= I(row, row) * tmp; // I中col取row就行
            for (int col = 0; col < m_Col; col++)
                result(cur_row, col) -= result(row, col) * tmp;
        }
    }

    return result;
}

// 横向拼接
Matrix cat(const Matrix &m1, const Matrix &m2)
{
    assert(m1.col() == m2.col());

    int m_Row1 = m1.m_Row;
    int m_Row2 = m2.m_Row;
    int m_Col = m1.m_Col;

    Matrix result(m_Row1 + m_Row2, m_Col);
    for (int row = 0; row < m_Row1; row++)
        for (int col = 0; col < m_Col; col++)
            result(row, col) = m1(row, col);
    for (int row = m_Row1; row < m_Row1 + m_Row2; row++)
        for (int col = 0; col < m_Col; col++)
            result(row, col) = m2(row - m_Row1, col);
    return result;
}

Matrix diag(const Matrix &m)
{
    int n = m.length();
    Matrix result(n, n);
    for (int i = 0; i < n; ++i)
        result(i, i) = m[i];
    return result;
}

Matrix sqrt(const Matrix &m)
{
    Matrix result(m);
    for (int i = 0; i < result.numel(); i++)
        result[i] = sqrt(m[i]);
    return result;
}

Matrix sum(const Matrix &m, int dim)
{
    int m_Row = m.row();
    int m_Col = m.col();
    if (dim == 1) // 按列求和
    {
        Matrix result(1, m_Col);
        for (int col = 0; col < m_Col; ++col)
        {
            double data = 0.0;
            for (int row = 0; row < m_Row; ++row)
                data += m(row, col);
            result[col] = data;
        }
        return result;
    }
    else // 按行求和
    {
        Matrix result(m_Row, 1);
        for (int row = 0; row < m_Row; ++row)
        {
            double data = 0.0;
            for (int col = 0; col < m_Col; ++col)
                data += m(row, col);
            result[row] = data;
        }
        return result;
    }
}

Matrix mul(const Matrix &m1, const Matrix &m2)
{
    assert(m1.row() == m2.row() && m1.col() == m2.col());

    Matrix result(m1.m_Row, m2.m_Col);
    for (int i = 0; i < m1.numel(); ++i)
        result[i] = m1[i] * m2[i];
    return result;
}

Matrix operator+(const Matrix &m1, const Matrix &m2)
{
    assert(m1.row() == m2.row() && m1.col() == m2.col());

    Matrix result(m1.m_Row, m2.m_Col);
    for (int i = 0; i < m1.numel(); ++i)
        result[i] = m1[i] + m2[i];
    return result;
}

Matrix operator-(const Matrix &m1, const Matrix &m2)
{
    assert(m1.row() == m2.row() && m1.col() == m2.col());

    Matrix result(m1.m_Row, m2.m_Col);
    for (int i = 0; i < m1.numel(); ++i)
        result[i] = m1[i] - m2[i];
    return result;
}

Matrix operator*(double data, const Matrix &m)
{
    Matrix result(m);
    for (int i = 0; i < result.numel(); i++)
        result[i] *= data;
    return result;
}

Matrix operator*(const Matrix &m1, const Matrix &m2)
{
    assert(m1.col() == m2.row());

    Matrix result(m1.m_Row, m2.m_Col);
    for (int row = 0; row < m1.m_Row; ++row)
    {
        for (int col = 0; col < m2.m_Col; ++col)
        {
            double data = 0.0;
            for (int k = 0; k < m1.m_Col; k++)
                data += m1(row, k) * m2(k, col);
            result(row, col) = data;
        }
    }
    return result;
}

ostream &operator<<(ostream &_out, const Matrix &m)
{
    _out << "[\n";
    for (int i = 0; i < m.row(); ++i)
    {
        _out << " ";
        for (int j = 0; j < m.col(); ++j)
            _out << m(i, j) << ", ";
        _out << "\n";
    }
    _out << "]\n";

    return _out;
}

void Matrix::swapRow(int iRow, int jRow)
{
    if (iRow == jRow)
        return;

    for (int col = 0; col < m_Col; ++col)
    {
        double temp = operator()(iRow, col);
        operator()(iRow, col) = operator()(jRow, col);
        operator()(jRow, col) = temp;
    }
}
