#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
using namespace std;

class Matrix
{
public:
    explicit Matrix(int row=0, int col=0, double init = 0.0);
    explicit Matrix(const vector<double> &data);
    Matrix(const Matrix &m);
    Matrix &operator=(const Matrix &m);
    Matrix(Matrix &&m);

    static Matrix zeros(int row, int col) { return Matrix(row, col); }
    static Matrix ones(int row, int col) { return Matrix(row, col, 1.0); }
    static Matrix diag(int n, double init = 1.0)
    {
        Matrix result(n, n);
        for (int i = 0; i < n; ++i)
            result.set(i, i, init);
        return result;
    }

    //高斯消元法求逆矩阵
    static Matrix inv(const Matrix &m){
        int m_Row = m.m_Row;
        int m_Col = m.m_Col;

        Matrix result(m);
        Matrix I = diag(m_Row, 1.0);

        for(int row=0; row<m_Row; row++){
            //选主元
            int pivotRow = row;
            for(int j=row+1; j<m_Row;j++){
                if( std::abs( result(j, row) ) > std::abs( result(pivotRow, row) ) )
                    pivotRow = j;
            }
            result.swapRow(row, pivotRow);
            I.swapRow(row, pivotRow);

            //每行第一个系数归一化
            double tmp = result(row, row);
            for(int j=row; j<m_Col;j++){
                result.set(row, j, result.get(row, j) / tmp);
                I.set(row, j, I.get(row, j) / tmp);
            }

            //按行消元
            for(int cur_row=row+1; cur_row<m_Row;cur_row++){
                double tmp = result.get(cur_row, row);
                double I_tmp = I.get(cur_row, row);
                for(int col=row; col<m_Col;col++){
                    result.set(cur_row, col, result.get(cur_row, col) - result.get(row, col) * tmp);
                    I.set(cur_row, col, I.get(cur_row, col) - I.get(row, col) * I_tmp);
                }
            }
        }

        for(int row=m_Row-1; row>=1; row--){
            int col = row;
            for(int cur_row=row-1; cur_row>=0;cur_row--){
                I.set(cur_row, col, I.get(row, col) - I.get(cur_row, col) );
            }
        }
        return I;
    }

    inline bool valid(){return m_Data != nullptr;}
    inline bool initialize(int row, int col, double init=0.0)
    {
        if(row>0 && col>0){
            m_Row = row;
            m_Col = col;
            m_Data.reset(new double[m_Row * m_Col]);
            for (int i = 0; i < m_Row * m_Col; ++i)
                m_Data[i] = init;
        }
    }
    inline void clear(){m_Row=0;m_Col=0;m_Data=nullptr;}

    inline int row() const { return m_Row; }
    inline int col() const { return m_Col; }
    inline int number() const {return m_Col * m_Row;}
    inline bool isSquareMatrix() const { return m_Row == m_Col; }
    inline int length() const { return m_Row > m_Col ? m_Row : m_Col; }

    inline double get(int row, int col) const { return    m_Data[m_Col * row + col]; }
    inline void set(int row, int col, double init) { m_Data[m_Col * row + col] = init; }
    inline double operator()(int row, int col) const {return get(row, col);}
    inline double get(int number) const {return m_Data[number];}
    inline void set(int number, double init) { m_Data[number] = init; }

    Matrix mul(const Matrix &m) const;// 矩阵点乘
    Matrix div(const Matrix &m) const;// 矩阵点除

    // 矩阵运算
    Matrix operator+(double data) const;
    Matrix operator-(double data) const;
    Matrix operator*(double data) const;
    Matrix operator/(double data) const;
    friend Matrix operator+(double data, const Matrix &m);
    friend Matrix operator-(double data, const Matrix &m);
    friend Matrix operator*(double data, const Matrix &m);
    friend Matrix operator/(double data, const Matrix &m);
    friend Matrix operator+(const Matrix &m1, const Matrix &m2);
    friend Matrix operator-(const Matrix &m1, const Matrix &m2);
    friend Matrix operator*(const Matrix &m1, const Matrix &m2);
    friend Matrix operator/(const Matrix &m1, const Matrix &m2);

    Matrix T() const;// 矩阵转置
    Matrix diag()
    {
        int n = length();
        Matrix result(n, n);
        for (int i = 0; i < n; ++i)
            result.set(i, i, get(i));
        return result;
    }

    Matrix cat(const Matrix &m){
        int row = m_Row + m.m_Row;
        int col = m_Col;

        Matrix result(row, col);
        for(row=0; row<m_Row;row++){
            for(col=0;col<m_Col;col++){
                result.set(row,col,get(row,col));
            }
        }
        for(row=m_Row; row<m_Row + m.m_Row;row++){
            for(col=0;col<m_Col;col++){
                result.set(row, col, m.get(row-m_Row, col));
            }
        }
    }

    Matrix sum(int dim) const;
    Matrix abs() const;
    Matrix sqrt() const;

    //矩阵交换两行
    void swapRow(int iRow, int jRow);

    //矩阵行列式
    double getDet() const;

    //矩阵点乘后求和
    double dot(const Matrix &m) const;

    //LU分解求逆矩阵
    Matrix inv_LU() const;

    //随机矩阵
    static Matrix matrixRand(int row, int col, double start, double end);

    //输出矩阵
    friend ostream &operator<<(ostream &_out, const Matrix &m);

private:
    int m_Row{0};
    int m_Col{0};
    shared_ptr<double[]> m_Data{nullptr};
};
