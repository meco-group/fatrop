#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hxx"

using namespace fatrop;

template class fatrop::IpEqMultInitializer<DenseType>;
