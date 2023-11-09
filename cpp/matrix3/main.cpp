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

#include "ExpressValue.h"

static int readCsv(const string &name, Matrix &X, Matrix &Y)
{
    ifstream csv_data(name, ios::in);
    if (!csv_data.is_open())
        return 1;

    string line;
    vector<string> lines;
    while (getline(csv_data, line))
        lines.emplace_back(line);

    int number = lines.size();
    X = move(Matrix::zeros(number, 1));
    Y = move(Matrix::zeros(number, 1));
    // cout << "number: " << number << endl;

    string token;
    for (int i = 0; i < number; i++)
    {
        istringstream iss(lines[i]);
        getline(iss, token, ',');
        X[i] = atof(token.c_str());

        getline(iss, token, ',');
        Y[i] = atof(token.c_str());
    }

    csv_data.close();
    return 0;
}

static Matrix myLMfit(const string &str, int n, const Matrix &X, const Matrix &Y, int maxiter = 200)
{
    const double betatol = 1.00e-08;
    const double rtol = 1.00e-08;
    const double eps = 2.220446049250313e-16;
    const double sqrteps = 1.490116119384766e-08;

    Matrix beta(n, 1, 1.0);
    ExpressValue ev(str, n);

    Matrix yfit = ev.F2(beta, X); // 求值
    Matrix r = Y - yfit;          // 求差
    double sse = r.sum();         // 误差平方和

    double lambda = 0.01;
    for (int i = 0; i < maxiter; i++)
    {
        Matrix betaold = beta;
        double sseold = sse;

        Matrix J = ev.statjacobian(X, beta);
        Matrix diagJtJ = sum(mul(J, J), 1);

        Matrix Jplus = cat(J, diag(sqrt(lambda * diagJtJ)));
        Matrix rplus = cat(r, Matrix::zeros(beta.length(), 1));

        Matrix JplusT = Jplus.T();
        Matrix step = (JplusT * Jplus).inv() * JplusT * rplus;
        // step = (JplusT * Jplus + lambda * diagJtJ).inv() * JplusT * rplus;

        beta = beta + step;
        yfit = ev.F2(beta, X); // 求值
        r = Y - yfit;          // 求差
        sse = r.sum();         // 误差平方和

        if (sse < sseold)
        {
            lambda = max(0.1 * lambda, eps);
        }
        else
        {
            while (sse > sseold)
            {
                lambda = 10 * lambda;
                if (lambda > 1e16)
                    return beta;

                Jplus = cat(J, diag(sqrt(lambda * diagJtJ)));
                JplusT = Jplus.T();
                step = (JplusT * Jplus).inv() * JplusT * rplus;
                // step = (JplusT * Jplus + lambda * diagJtJ).inv() * JplusT * rplus;

                beta = betaold + step;
                yfit = ev.F2(beta, X); // 求值
                r = Y - yfit;          // 求差
                sse = r.sum();         // 误差平方和
            }
        }
        if (sqrt(step.sum()) < betatol * (sqrteps + sqrt(beta.sum())) || abs(sse - sseold) <= rtol * sse)
            break;
    }

    // 测试拟合效果
    /*
        int m_Row = X.length();
        int m_Col = beta.length();
        for (int i = 0; i < m_Col; i++)
            ev.setValue("beta" + to_string(i), beta[i]);

        for (int i = 0; i < m_Row; i++)
        {
            ev.setValue("x", X[i]);
            double value = ev.eval();
            cout << "x=" << X[i] << " y=" << Y[i] << " yfit=" << value << endl;
        }
    */
    return beta;
}

int main()
{
    string expr{"beta0 + (beta1 * x + beta2) / ((x + beta3) * x + beta4) + (beta5 * x + beta6) / ((x + beta7) * x + beta8)"};
    int n = 9;
    // ExpressValue ev(expr, n);
    // cout << "value = " << ev.eval() << endl;

    Matrix X;
    Matrix Y;
    string name{"ChinaSteel_35CS250H.csv"};
    readCsv(name, X, Y);

    Matrix beta = myLMfit(expr, n, X, Y);
    // cout << beta;

    return 0;
}

// g++ -o main main.cpp Matrix.cpp ExpressValue.cpp lua.dll && main.exe
