#include "fatrop/context/generic.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <iostream>

using namespace fatrop;
using namespace std;
int main() {
    MatRealAllocated A(5, 5);
    MatRealAllocated B(5, 5);
    auto Ablock = A.block(4, 4, 1, 1).block(2, 2, 1, 1).diagonal();
    A = 1.;
    cout << " A: \n" << A << endl;
    A.diagonal() = VecRealScalar(A.m(), 2);
    cout << " A: \n" << A << endl;
    A.diagonal() = A.diagonal() + 1.;
    cout << " A: \n" << A << endl;
    cout << " 5*A diagonal: \n" << 5.*A.diagonal() << endl;
    cout << " - 1./A diagonal**2: \n" << -1. / (A.diagonal()*A.diagonal()) << endl;
    cout << " sum A diagonal: " << sum(A.diagonal()) << endl;
    A(2,2) = 10;
    cout << " A(1:, 1:) : \n" << A.block(4, 4, 1, 1) << endl;
    cout << " A(1:, 1:) diagonal: \n" << A.block(4, 4, 1, 1).diagonal() << endl;
    cout << " A(1:, 1:) row 1 : \n" << A.block(4, 4, 1, 1).row(1) << endl;
    cout << " A(1:, 1:) col 2 : \n" << A.block(4, 4, 1, 1).col(2) << endl;
    cout << "log(A(1:, 1:)) col 2 : \n" << log(A.block(4, 4, 1, 1).col(2)) << endl;
    cout << "Ablock: \n" << Ablock.block(2, 0) << endl;
    return 0;
}
