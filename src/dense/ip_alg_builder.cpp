#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_solver_resto.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/pd_system_resto.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_alg_builder.hxx"

namespace fatrop
{
    template class IpAlgBuilder<DenseType>;
} // namespace fatrop
