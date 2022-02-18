#include "sparse/FatropSparse.hpp"
using namespace fatrop;
namespace fatrop{
    fe_sp operator+(const fe_sp &fe1, const fe_sp &fe2)
    {
        fe_sp res = make_shared<FatropSum1>(fe1, fe2);
        return res;
    }
    fe_sp operator*(const Eig &mat, var_sp var)
    {
        fe_sp res = make_shared<MatrixVector>(mat, var);
        return res;
    };
}
