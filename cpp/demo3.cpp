#include <iostream>
using namespace std;

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main()
{
    Matrix3d m = Matrix3d::Random();
    m = (m + Matrix3d::Constant(1.2)) * 50;
    std::cout << "m = " << m << std::endl;

    Vector3d v(1, 2, 3);
    cout << "m * v = " << m * v << endl;

    return 0;
}