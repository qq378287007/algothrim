#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cassert>
#include <functional>
using namespace std;

#include "Matrix.h"

void readCsv(const string &name, Matrix &X, Matrix &Y)
{
    ifstream csv_data(name, ios::in);
    if (!csv_data.is_open())
    {
        cout << "Error: opening file fail" << endl;
        exit(1);
    }

    string line;
    vector<string> lines;
    while (getline(csv_data, line))
        lines.emplace_back(line);

    int number = lines.size();
    X = move(Matrix::zeros(number, 1));
    Y = move(Matrix::zeros(number, 1));

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
}

using func = Matrix (*)(const Matrix &, const Matrix &);

Matrix model(const Matrix &beta, const Matrix &X)
{
    assert(beta.length() >= 9);

    int number = X.length();
    Matrix Y(number, 1);
    for (int i = 0; i < number; i++)
        Y[i] = beta[0] + (beta[1] * X[i] + beta[2]) / ((X[i] + beta[3]) * X[i] + beta[4]) + (beta[5] * X[i] + beta[6]) / ((X[i] + beta[7]) * X[i] + beta[8]);
    return Y;
}

Matrix statjacobian(const Matrix &X, const Matrix &beta, func F2 = model)
{
    const double eps = 2.220446049250313e-16;
    const double DerivStep = pow(eps, 1.0 / 3.0);
    const double nTheta = sqrt(beta.sum());
    const double DerivStep2 = DerivStep * (nTheta + (nTheta == 0.0));

    int m_Row = X.length();
    int m_Col = beta.length();
    Matrix J(m_Row, m_Col);
    auto F1 = bind(F2, placeholders::_1, X); // model（F2）第一个参数beta待定，第二个参数X绑定

    Matrix thetaNew;
    Matrix yplus;
    Matrix dy;

    Matrix y0 = F1(beta);
    for (int i = 0; i < m_Col; i++)
    {
        double delta = DerivStep * beta[i];
        if (delta == 0.0)
            delta = DerivStep2;
        thetaNew = beta;
        thetaNew[i] += delta;
        yplus = F1(thetaNew);
        dy = 1.0 / delta * (yplus - y0);

        for (int row = 0; row < m_Row; row++)
            J(row, i) = dy[row];
    }

    return J;
}

Matrix myLMfit(const Matrix &X, const Matrix &Y, func F2 = model, int maxiter = 200)
{
    const double betatol = 1.00E-08;
    const double rtol = 1.00E-08;
    const double eps = 2.220446049250313e-16;
    const double sqrteps = 1.490116119384766e-08;

    Matrix beta(9, 1, 1.0);

    Matrix yfit = F2(beta, X); // 求值
    Matrix r = Y - yfit;       // 求差
    double sse = r.sum();      // 误差平方和

    Matrix betaold;
    double sseold;
    int number = X.length();
    Matrix J(number, 9);
    Matrix diagJtJ;
    Matrix Jplus;
    Matrix rplus;
    Matrix step;
    Matrix JplusT;

    double lambda = 0.01;
    for (int i = 0; i < maxiter; i++)
    {
        betaold = beta;
        sseold = sse;

        J = statjacobian(X, beta, F2);

        diagJtJ = sum(mul(J, J), 1);

        Jplus = cat(J, diag(sqrt(lambda * diagJtJ)));
        rplus = cat(r, Matrix::zeros(beta.length(), 1));

        JplusT = Jplus.T();
        step = (JplusT * Jplus).inv() * JplusT * rplus;
        // step = (JplusT * Jplus + lambda * diagJtJ).inv() * JplusT * rplus;

        beta = beta + step;

        yfit = F2(beta, X); // 求值
        r = Y - yfit;       // 求差
        sse = r.sum();      // 误差平方和

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
                yfit = F2(beta, X); // 求值
                r = Y - yfit;       // 求差
                sse = r.sum();      // 误差平方和
            }
        }
        if (sqrt(step.sum()) < betatol * (sqrteps + sqrt(beta.sum())) || abs(sse - sseold) <= rtol * sse)
            break;
    }

    return beta;
}

int main()
{
    Matrix X;
    Matrix Y;
    string name{"ChinaSteel_35CS250H.csv"};
    readCsv(name, X, Y);
    // cout<<X;
    // cout<<Y;

    Matrix beta = myLMfit(X, Y);
    cout << beta;
    cout << "Over!\n";
    return 0;
}
