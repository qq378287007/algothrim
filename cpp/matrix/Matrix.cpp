#include "Matrix.h"

Matrix::Matrix(int row, int col, double init)
{
    initialize(row, col, init);
}

Matrix::Matrix(const vector<double> &data)
{
    int number = data.size();
    m_Row = number;
    m_Col = 1;

    m_Data.reset(new double[number]);
    for (int i = 0; i < number; ++i)
        m_Data[i] = data[i];
}

Matrix::Matrix(const Matrix &m)
{
    m_Row = m.row();
    m_Col = m.col();

    m_Data.reset(new double[m_Row * m_Col]);
    for (int row = 0; row < m_Row; ++row)
        for (int col = 0; col < m_Col; ++col)
            set(row, col, m.get(row, col));
}

Matrix &Matrix::operator=(const Matrix &m)
{
    if (m.m_Data == m_Data)
        return *this;

    if (m_Row != m.m_Row || m_Col != m.m_Col)
        throw "Error in operator=(): size of matrix is not equal.";

    for (int row = 0; row < m_Row; ++row)
        for (int col = 0; col < m_Col; ++col)
            set(row, col, m.get(row, col));
    return *this;
}

Matrix::Matrix(Matrix &&m)
{
    m_Row = m.m_Row;
    m_Col = m.m_Col;
    m_Data = m.m_Data;
    m.clear();
}

Matrix Matrix::mul(const Matrix &m) const
{
    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
    {
        for (int col = 0; col < m_Col; ++col)
        {
            double data = get(row, col) * m.get(row, col);
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix Matrix::div(const Matrix &m) const
{
    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
    {
        for (int col = 0; col < m_Col; ++col)
        {
            double data = get(row, col) / m.get(row, col);
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix Matrix::operator+(double data) const
{
    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
    {
        for (int col = 0; col < m_Col; ++col)
        {
            double m_Data = get(row, col) + data;
            result.set(row, col, m_Data);
        }
    }
    return result;
}

Matrix Matrix::operator-(double data) const
{
    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
    {
        for (int col = 0; col < m_Col; ++col)
        {
            double m_Data = get(row, col) - data;
            result.set(row, col, m_Data);
        }
    }
    return result;
}

Matrix Matrix::operator*(double scale) const
{
    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
    {
        for (int col = 0; col < m_Col; ++col)
        {
            double data = get(row, col) * scale;
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix Matrix::operator/(double scale) const
{

    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
    {
        for (int col = 0; col < m_Col; ++col)
        {
            double data = get(row, col) / scale;
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix operator+(double data, const Matrix &m)
{

    return m + data;
}

Matrix operator-(double data, const Matrix &m)
{

    return m - data;
}

Matrix operator*(double data, const Matrix &m)
{
    return m * data;
}

Matrix operator/(double data, const Matrix &m)
{
    return m / data;
}

Matrix operator+(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1.m_Row, m2.m_Col);
    for (int row = 0; row < m1.m_Row; ++row)
    {
        for (int col = 0; col < m2.m_Col; ++col)
        {
            double data = m1.get(row, col) + m2.get(row, col);
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix operator-(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1.m_Row, m2.m_Col);
    for (int row = 0; row < m1.m_Row; ++row)
    {
        for (int col = 0; col < m2.m_Col; ++col)
        {
            double data = m1.get(row, col) - m2.get(row, col);
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix operator*(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1.m_Row, m2.m_Col);
    for (int row = 0; row < m1.m_Row; ++row)
    {
        for (int col = 0; col < m2.m_Col; ++col)
        {
            double data = 0.0;
            for(int k=0; k<m1.m_Col; k++){
                data += m1.get(row, k) * m2.get(k, col);
            }
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix operator/(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1.m_Row, m2.m_Col);
    for (int row = 0; row < m1.m_Row; ++row)
    {
        for (int col = 0; col < m2.m_Col; ++col)
        {
            double data = m1.get(row, col) / m2.get(row, col);
            result.set(row, col, data);
        }
    }
    return result;
}

Matrix Matrix::T() const
{
    Matrix result(m_Col, m_Row);
    for (int row = 0; row < m_Col; ++row)
        for (int col = 0; col < m_Row; ++col)
            result.set(row, col, get(col, row));
    return result;
}

Matrix Matrix::sum(int dim) const
{

    if (dim == 1) // 按列求和
    {
        Matrix result(1, m_Col, 0.0);
        for (int col = 0; col < m_Col; ++col)
        {
            double data = 0.0;
            for (int row = 0; row < m_Row; ++row)
                data += get(row, col);
            result.set(0, col, data);
        }
        return result;
    }
    else // 按行求和
    {
        Matrix result(m_Row, 1, 0.0);
        for (int row = 0; row < m_Row; ++row)
        {
            double data = 0.0;
            for (int col = 0; col < m_Col; ++col)
                data += get(row, col);
            result.set(row, 0, data);
        }
        return result;
    }
}


Matrix Matrix::abs() const {

    Matrix result(m_Row, m_Col, 0.0);
    for (int row = 0; row < m_Row; ++row)
        for (int col = 0; col < m_Col; ++col)
            result.set(row, col, std::abs(get(row, col)));
    return result;

}

Matrix Matrix::sqrt() const {
    Matrix result(m_Row, m_Col);
    for (int row = 0; row < m_Row; ++row)
        for (int col = 0; col < m_Col; ++col)
            result.set(row, col, std::sqrt(get(row, col)));
    return result;
}

void Matrix::swapRow(int iRow, int jRow)
{
    if(iRow == jRow)
        return;

    double *pI = &m_Data[m_Col * iRow];
    double *pJ = &m_Data[m_Col * jRow];
    for (int col = 0; col < m_Col; ++col)
    {
        double temp = pI[col];
        pI[col] = pJ[col];
        pJ[col] = temp;
    }
}


#define EPSINON 1.0e-14
static inline bool isEqual(const double &num1, const double &num2)
{
    return (num1 - num2 <= EPSINON && num1 - num2 >= -EPSINON);
}

double Matrix::getDet() const
{
    /*
    if (!isSquareMatrix())
        throw "Error in getDet(): the matrix is not square.";
    */

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

double Matrix::dot(const Matrix &m) const
{
    Matrix result = mul(m).sum(1).sum(2);
    double r = *result.m_Data;

    return r;
}

Matrix Matrix::inv_LU() const
{
    double det = getDet();
    /*
    if (det >= -EPSINON && det <= EPSINON)
        throw "Error in inv_LU(): the determinant is equal to 0.";
*/
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

Matrix Matrix::matrixRand(int row, int col, double start, double end)
{
    Matrix result(row, col);

    auto pR = result.m_Data;
    for (int i = 0; i < row * col; ++i)
    {
        //	srand((unsigned)time(NULL));
        pR[i] = rand() / double(RAND_MAX) * (end - start) + start;
    }
    return result;
}

ostream &operator<<(ostream &_out, const Matrix &m)
{
    auto p = m.m_Data;

    _out << "[\n";
    for (int i = 0; i < m.m_Row; ++i)
    {
        _out << "  ";
        for (int j = 0; j < m.m_Col; ++j)
            _out << p[ i * m.m_Col + j] << ",\t ";
        _out << "\n";
    }
    _out << "]\n";

    return _out;
}
