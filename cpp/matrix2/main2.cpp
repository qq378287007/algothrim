#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include <istream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cassert>
using namespace std;

#include "libLua/lua.hpp"

static int l_log10(lua_State *L)
{
    double d = luaL_checknumber(L, 1);
    lua_pushnumber(L, log10(d));
    return 1;
}

class Matrix
{
public:
    Matrix(int row = 0, int col = 0, double init = 0.0)
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
    Matrix &operator=(const Matrix &m)
    {
        if (m_Data.get() != m.data().get() && m.valid())
        {
            m_Row = m.row();
            m_Col = m.col();
            m_Data.reset(new double[m_Row * m_Col]);
            for (int row = 0; row < m_Row; ++row)
                for (int col = 0; col < m_Col; ++col)
                    operator()(row, col) = m(row, col);
        }
        return *this;
    }

public:
    inline int row() const { return m_Row; }
    inline int col() const { return m_Col; }
    inline shared_ptr<double[]> data() const { return m_Data; }
    inline int numel() const { return m_Col * m_Row; }
    inline int length() const { return m_Row > m_Col ? m_Row : m_Col; }
    inline bool valid() const { return m_Data != nullptr; }

public:
    inline double &operator()(int row, int col) { return m_Data[m_Col * row + col]; }
    inline double operator()(int row, int col) const { return m_Data[m_Col * row + col]; }
    inline double &operator[](int number) { return m_Data[number]; }
    inline double operator[](int number) const { return m_Data[number]; }

public:
    static Matrix zeros(int row, int col) { return Matrix(row, col); }

public:
    double sum() const
    {
        double data = 0.0;
        for (int i = 0; i < m_Col * m_Row; i++)
            data += m_Data[i];
        return data;
    }

public:
    friend Matrix mul(const Matrix &m1, const Matrix &m2);
    friend Matrix operator+(const Matrix &m1, const Matrix &m2);
    friend Matrix operator-(const Matrix &m1, const Matrix &m2);
    friend Matrix operator*(double data, const Matrix &m);
    friend Matrix operator*(const Matrix &m1, const Matrix &m);

private:
    int m_Row{0};
    int m_Col{0};
    shared_ptr<double[]> m_Data{nullptr};
};

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
    // assert(m1.row() == m2.row() && m1.col() == m2.col());
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
    // assert(m1.col() == m2.row());

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

class ExpressValue
{
public:
    ExpressValue(const string &expr, int n)
    {
        m_LuaStater = luaL_newstate();
        luaL_openlibs(m_LuaStater);

        luaL_dostring(m_LuaStater, "sin = math.sin; cos = math.cos;"
                                   "tan = math.tan; abs = math.abs;"
                                   "acos = math.acos; asin = math.asin;"
                                   "atan = math.atan; exp = math.exp;"
                                   "log = math.log; "
                                   "sqrt = math.sqrt; int = math.floor;"
                                   "min = math.min; max = math.max;");

        lua_pushcfunction(m_LuaStater, l_log10);
        lua_setglobal(m_LuaStater, "log10");

        lua_checkstack(m_LuaStater, 100);

        setStr(expr);
        setValue("x", 1.0);
        setNum(n);
    }
    ~ExpressValue()
    {
        lua_close(m_LuaStater);
    }

    void setStr(const string &expr)
    {
        str = expr;
    }
    void setNum(int n)
    {
        for (int i = 0; i < n; i++)
            setValue("beta" + to_string(i), 1.0);
    }

    void setValue(const string &key, double value) const
    {
        const char *name = key.c_str();

        lua_getglobal(m_LuaStater, name);
        lua_pop(m_LuaStater, 1);
        lua_pushnumber(m_LuaStater, value);
        lua_setglobal(m_LuaStater, name);
    }

    void setValue(const vector<pair<string, double>> &vp) const
    {
        for (const pair<string, double> &p : vp)
            setValue(p.first, p.second);
    }

    double eval() const
    {
        luaL_dostring(m_LuaStater, ("return 0 + " + str).c_str());
        double value = lua_tonumber(m_LuaStater, -1);
        lua_pop(m_LuaStater, 1);
        return value;
    }

    Matrix F2(const Matrix &beta, const Matrix &X)
    {
        int m_Row = X.length();
        int m_Col = beta.length();

        Matrix Y(m_Row, m_Col);

        for (int i = 0; i < m_Col; i++)
            setValue("beta" + to_string(i), beta[i]);

        for (int i = 0; i < m_Row; i++)
        {
            setValue("x", X[i]);
            Y[i] = eval();
        }

        return Y;
    }

    Matrix statjacobian(const Matrix &X, const Matrix &beta)
    {
        const double eps = 2.220446049250313e-16;
        const double DerivStep = pow(eps, 1.0 / 3.0);
        const double nTheta = sqrt(beta.sum());
        const double DerivStep2 = DerivStep * (nTheta + (nTheta == 0.0));

        int m_Row = X.length();
        int m_Col = beta.length();

        Matrix J(m_Row, m_Col);

        for (int i = 0; i < m_Col; i++)
            setValue("beta" + to_string(i), beta[i]);

        Matrix thetaNew(m_Col);
        Matrix dy(m_Row);
        Matrix y0(m_Row);
        Matrix yplus(m_Row);
        for (int i = 0; i < m_Row; i++)
        {
            setValue("x", X[i]);
            y0[i] = eval();
        }

        for (int i = 0; i < m_Col; i++)
        {
            double delta = DerivStep * beta[i];
            if (delta == 0.0)
                delta = DerivStep2;
            thetaNew = beta;
            thetaNew[i] += delta;

            for (int j = 0; j < m_Col; j++)
                setValue("beta" + to_string(j), thetaNew[j]);
            for (int j = 0; j < m_Row; j++)
            {
                setValue("x", X[j]);
                yplus[j] = eval();
            }

            dy = 1.0 / delta * (yplus - y0);

            for (int row = 0; row < m_Row; row++)
                J(row, i) = dy[row];
        }

        return J;
    }

private:
    lua_State *m_LuaStater;
    string str;
};

Matrix myLMfit(const string &str, int n, const Matrix &X, const Matrix &Y, int maxiter = 200)
{
    const double betatol = 1.00E-08;
    const double rtol = 1.00E-08;
    const double eps = 2.220446049250313e-16;
    const double sqrteps = 1.490116119384766e-08;

    Matrix beta(n, 1, 1.0);
    ExpressValue ev(str, n);

    Matrix yfit = ev.F2(beta, X); // 求值
    Matrix r = Y - yfit;          // 求差
    double sse = r.sum();         // 误差平方和

    int m_Row = X.length();
    int m_Col = beta.length();

    Matrix step;
    Matrix JplusT;

    double lambda = 0.01;
    for (int i = 0; i < maxiter; i++)
    {
        Matrix betaold = beta;
        double sseold = sse;

        Matrix J = ev.statjacobian(X, beta);
        Matrix diagJtJ = sum(mul(J, J), 1);

        Matrix Jplus = cat(J, diag(sqrt(lambda * diagJtJ)));
        Matrix rplus = cat(r, Matrix::zeros(beta.length(), 1));
    }

    return beta;
}

int readCsv(const string &name, vector<double> &X, vector<double> &Y)
{
    ifstream csv_data(name, ios::in);
    if (!csv_data.is_open())
        return 1;

    vector<string> lines;
    string line;
    while (getline(csv_data, line))
        lines.emplace_back(line);

    int number = lines.size();
    X.resize(number);
    Y.resize(number);

    for (int i = 0; i < number; i++)
    {
        istringstream iss(lines[i]);
        string token;

        getline(iss, token, ',');
        X[i] = atof(token.c_str());

        getline(iss, token, ',');
        Y[i] = atof(token.c_str());
    }

    csv_data.close();
    return 0;
}

int main()
{
    const int n = 9;
    ExpressValue ev("beta0 + (beta1 * x + beta2) / ((x + beta3) * x + beta4) + (beta5 * x + beta6) / ((x + beta7) * x + beta8)", n);
    double value = ev.eval();
    cout << "value = " << value << endl;

    vector<double> X;
    vector<double> Y;
    string name{"ChinaSteel_35CS250H.csv"};
    readCsv(name, X, Y);

    /*
    for (double x : X)
    {
        ev.setValue("x", x);
        double value = ev.eval();
        cout << "value = " << value << endl;
    }
    */

    const int N = X.size();

    // 求值
    vector<double> beta(n, 1.0);
    for (int i = 0; i < n; i++)
        ev.setValue("beta" + to_string(i), beta[i]);

    vector<double> yfit(N);
    vector<double> r(N);
    double sse = 0.0;
    for (int i = 0; i < N; i++)
    {
        ev.setValue("x", X[i]);
        yfit[i] = ev.eval();
        r[i] = Y[i] - yfit[i];
        sse += r[i] * r[i];
    }

    vector<double> betaold(n);
    double sseold;
    int maxiter = 200;
    double lambda = 0.01;
    for (int i = 0; i < maxiter; i++)
    {
        betaold.assign(beta.cbegin(), beta.cend());
        sseold = sse;
    }

    return 0;
}

// g++ -o main2 main2.cpp lua.dll && main2.exe
