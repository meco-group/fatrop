#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hxx"

using namespace fatrop;

// For dense problems the "inequality" constraint-violation block sits at
// offset_g_eq_slack inside the constraint-violation vector.
template <>
const VecRealView IpIterate<DenseType>::constr_viol_ineq()
{
    return constr_viol().block(info_->number_of_g_eq_slack, info_->offset_g_eq_slack);
}

// explicit template instantiation
template class fatrop::IpIterate<DenseType>;
