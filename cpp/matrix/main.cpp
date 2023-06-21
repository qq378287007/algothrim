#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
using namespace std;

#include "Matrix.h"

void readCsv(const string &name, vector<double> &X, vector<double> &Y)
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
}

using func = Matrix (*)(const Matrix &, const Matrix &);

Matrix model(const Matrix &beta, const Matrix &X)
{
    int number = X.row();
    Matrix Y(number, 1);

    double beta0 = beta.get(0, 0);
    double beta1 = beta.get(1, 0);
    double beta2 = beta.get(2, 0);
    double beta3 = beta.get(3, 0);
    double beta4 = beta.get(4, 0);
    double beta5 = beta.get(5, 0);
    double beta6 = beta.get(6, 0);
    double beta7 = beta.get(7, 0);
    double beta8 = beta.get(8, 0);

    for (int row = 0; row < number; ++row)
    {
        for (int col = 0; col < 1; ++col)
        {
            double x = X.get(row, col);
            double data = beta0 + (beta1 * x + beta2) / (x * x + beta3 * x + beta4) + (beta5 * x + beta6) / (x * x + beta7 * x + beta8);
            Y.set(row, col, data);
        }
    }
    return Y;
}

Matrix statjacobian(const Matrix &X, const Matrix &beta){
    int m_Row = X.row();
    int m_Col = beta.row();
    Matrix J(m_Row, m_Col);
    Matrix y0 = model(beta, X);

    double DerivStep = 0.0;
    double nTheta = 0.0;
    double DerivStep2 = 0.0;

    for(int col=0; col<m_Col; col++){
        double delta = DerivStep * beta.get(col, 0);
        if(delta == 0.0)
            delta = DerivStep2;
        Matrix thetaNew = beta;
        thetaNew.set(col,0,delta);
        Matrix yplus = model(thetaNew, X);
        Matrix dy = (yplus - y0) / delta;
        for(int row=0; row<m_Row; row++)
            J.set(row, col, dy.get(m_Col, 0));
    }

    return J;
}

Matrix myLMfit(const Matrix &X, const Matrix &Y)
{
    int number = X.row();
    Matrix beta(9, 1, 1.0);

    const int maxiter = 200;

    const double betatol = 1.00E-08;
    const double rtol = 1.00E-08;
    const double sqrteps = 1.490116119384766e-08;
    const double eps = 2.220446049250313e-16;

    double lambda = 0.01;
    Matrix yfit = model(beta, X); // 求值

    Matrix r = Y - yfit;          // 求差
    Matrix sse = r.T() * r;       // 误差平方和


    for (int i = 0; i < maxiter; i++)
    {
        Matrix betaold = beta;
        Matrix sseold = sse;

        Matrix J(number, 9);

        Matrix diagJtJ =  J.mul(J).sum(1);
        //Matrix Jplus = J.cat( (lambda * diagJtJ).sqrt().diag() );

        //cout <<r.col();
       //Matrix rplus = r.cat(Matrix::zeros(beta.length(), 1));

        yfit = model(beta, X); // 求值
        r = Y - yfit;          // 求差
        sse = r.T() * r;       // 误差平方和

        lambda = max(0.1 * lambda, eps);
    }

    return beta;
}

int main()
{
    string name{"ChinaSteel_35CS250H.csv"};
    vector<double> v_X;
    vector<double> v_Y;
    readCsv(name, v_X, v_Y);
    Matrix X(v_X);
    Matrix Y(v_Y);

    myLMfit(X, Y);
    cout << "Over!\n";

    return 0;
}
