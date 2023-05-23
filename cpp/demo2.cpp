#include <iostream>
using namespace std;

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main()
{
    MatrixXd m = MatrixXd::Random(3, 3);
    m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;
    cout << "m = " << m << endl;
    
    VectorXd v(3);
    v << 1, 2, 3;

    cout << "m * v = " << m * v << endl;

    return 0;
}