#include <iostream>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

int main()
{
     MatrixXf A = MatrixXf::Random(3, 2);
     cout << "Here is the matrix A:\n"
          << A << endl;

     VectorXf b = VectorXf::Random(3);
     cout << "Here is the right hand side b:\n"
          << b << endl;

     // SVD
     cout << "The least-squares solution is:\n"
          << A.bdcSvd(ComputeThinU | ComputeThinV).solve(b) << endl;

     // QR
     cout << "The solution using the QR decomposition is:\n"
          << A.colPivHouseholderQr().solve(b) << endl;

     // 常规表达式
     // Ax = b is equivalent to solving the normal equation ATAx = ATb
     cout << "The solution using normal equations is:\n"
          << (A.transpose() * A).ldlt().solve(A.transpose() * b) << endl;
     return 0;
}