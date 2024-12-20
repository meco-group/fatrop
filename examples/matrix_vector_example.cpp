#include "fatrop/context/generic.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <iostream>

using namespace fatrop;
using namespace std;
int main() {
    MatrixAllocated A(5, 5);
    MatrixAllocated B(5, 5);
    A = 1.;
    cout << " A: \n" << A << endl;
    A.diagonal() = VecScalar(A.m(), 2);
    cout << " A: \n" << A << endl;
    A.diagonal() = A.diagonal() + VecScalar(A.m(), 1);
    cout << " A: \n" << A << endl;
    cout << " A diagonal: \n" << A.diagonal() << endl;
    cout << " sum A diagonal: " << sum(A.diagonal()) << endl;
    A(2,2) = 10;
    cout << " A(1:, 1:) : \n" << A.block(4, 4, 1, 1) << endl;
    cout << " A(1:, 1:) diagonal: \n" << A.block(4, 4, 1, 1).diagonal() << endl;
    cout << " A(1:, 1:) row 1 : \n" << A.block(4, 4, 1, 1).row(1) << endl;
    cout << " A(1:, 1:) col 2 : \n" << A.block(4, 4, 1, 1).col(2) << endl;
    cout << "log(A(1:, 1:)) col 2 : \n" << log(A.block(4, 4, 1, 1).col(2)) << endl;
    return 0;
}
